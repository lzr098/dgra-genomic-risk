"""L2 unit tests for batch 6: variant_filter, review_gate, myvariant, proxy_routes."""
import asyncio
import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock, AsyncMock

from dgra_variant_filter import (
    filter_variants,
    FILTER_PRESETS,
    _get_impact,
    _get_consequence_terms,
    _has_splice_consequence,
    _is_synonymous_in_tissue_gene,
    _is_clinvar_conflicting,
    get_tissue_relevant_genes,
)
from gpa_review_gate import GPAReviewGate, Severity, Issue
from dgra_myvariant import (
    _build_variant_id,
    _parse_gnomad_from_myvariant,
    _parse_clinvar_from_myvariant,
    _parse_cadd_from_myvariant,
    _parse_dbsnp_from_myvariant,
    MyVariantResult,
    apply_myvariant_results,
)
from gpa_proxy_routes import (
    ProxyRoute,
    ProxyRouteMap,
    _build_candidate_proxies,
    _probe_single,
    probe_api_routes,
    build_route_map,
)


# =============================================================================
# dgra_variant_filter
# =============================================================================

class TestVariantFilterHelpers:
    def test_get_impact(self):
        assert _get_impact({"IMPACT": "HIGH"}) == "HIGH"
        assert _get_impact({"IMPACT": " high "}) == "HIGH"
        assert _get_impact({}) == ""

    def test_get_consequence_terms(self):
        assert _get_consequence_terms({"Consequence": "missense_variant"}) == {"missense_variant"}
        assert _get_consequence_terms({"Consequence": "错义变异"}) == {"missense_variant"}
        assert _get_consequence_terms({}) == set()

    def test_has_splice_consequence(self):
        assert _has_splice_consequence({"splice_region_variant"})
        assert not _has_splice_consequence({"missense_variant"})

    def test_is_synonymous_in_tissue_gene(self):
        assert _is_synonymous_in_tissue_gene({"Consequence": "synonymous_variant", "GENE": "TP53"}, {"TP53"})
        assert not _is_synonymous_in_tissue_gene({"Consequence": "missense_variant", "GENE": "TP53"}, {"TP53"})
        assert not _is_synonymous_in_tissue_gene({"Consequence": "synonymous_variant", "GENE": "TP53"}, None)

    def test_is_clinvar_conflicting(self):
        assert _is_clinvar_conflicting({"CLIN_SIG": "Pathogenic/Benign"})
        assert _is_clinvar_conflicting({"CLIN_SIG": "致病, 良性"})
        assert not _is_clinvar_conflicting({"CLIN_SIG": "Pathogenic"})
        assert not _is_clinvar_conflicting({"CLIN_SIG": ""})


class TestFilterVariants:
    def mk(self, impact, cons, gene="GENE", clin=""):
        return {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": gene,
                "IMPACT": impact, "Consequence": cons, "CLIN_SIG": clin}

    def test_strict_preset(self):
        variants = [
            self.mk("HIGH", "stop_gained"),
            self.mk("MODERATE", "missense_variant"),
            self.mk("LOW", "synonymous_variant"),
        ]
        filtered, stats = filter_variants(variants, preset="strict")
        assert stats["output_count"] == 2
        assert stats["excluded"] == 1

    def test_clinical_preset_splice(self):
        variants = [self.mk("LOW", "splice_region_variant")]
        filtered, stats = filter_variants(variants, preset="clinical")
        assert stats["output_count"] == 1
        assert stats["splice_retained"] == 1

    def test_clinical_preset_synonymous_tissue(self):
        variants = [self.mk("LOW", "synonymous_variant", gene="TP53")]
        filtered, stats = filter_variants(variants, preset="clinical", tissue_relevant_genes={"TP53"})
        assert stats["output_count"] == 1
        assert stats["synonymous_tissue_retained"] == 1

    def test_clinical_preset_clinvar_conflict(self):
        variants = [self.mk("LOW", "synonymous_variant", clin="Pathogenic/Benign")]
        filtered, stats = filter_variants(variants, preset="clinical")
        assert stats["output_count"] == 1
        assert stats["clinvar_conflicting_retained"] == 1

    def test_broad_preset(self):
        variants = [self.mk("LOW", "synonymous_variant")]
        filtered, stats = filter_variants(variants, preset="broad")
        assert stats["output_count"] == 1

    def test_unknown_preset_raises(self):
        with pytest.raises(ValueError):
            filter_variants([], preset="unknown")

    def test_filter_metadata_added(self):
        variants = [self.mk("HIGH", "stop_gained")]
        filtered, _ = filter_variants(variants, preset="clinical")
        assert "_filter_reason" in filtered[0]
        assert "_filter_preset" in filtered[0]


class TestGetTissueRelevantGenes:
    def test_missing_file(self, tmp_path):
        result = get_tissue_relevant_genes("general", tmp_path / "no.json")
        assert result == set()

    def test_valid_profile(self, tmp_path):
        data = {"profiles": {"general": {"special_gene_lists": {"list1": ["TP53", "BRCA1"]}}}}
        path = tmp_path / "tissue.json"
        path.write_text(json.dumps(data))
        result = get_tissue_relevant_genes("general", path)
        assert result == {"TP53", "BRCA1"}

    def test_unknown_profile(self, tmp_path):
        data = {"profiles": {}}
        path = tmp_path / "tissue.json"
        path.write_text(json.dumps(data))
        assert get_tissue_relevant_genes("general", path) == set()


# =============================================================================
# gpa_review_gate
# =============================================================================

class TestReviewGateBlockers:
    def _run(self, source: str, tmp_path):
        f = tmp_path / "test.py"
        f.write_text(source)
        gate = GPAReviewGate([f])
        gate._check_file(f)
        return gate.issues

    def test_bare_except(self, tmp_path):
        issues = self._run("try:\n    pass\nexcept:\n    pass\n", tmp_path)
        assert any(i.code == "BARE_EXCEPT" for i in issues)

    def test_broad_except(self, tmp_path):
        issues = self._run("try:\n    pass\nexcept Exception:\n    pass\n", tmp_path)
        assert any(i.code == "BROAD_EXCEPT" for i in issues)

    def test_mutable_default_list(self, tmp_path):
        issues = self._run("def foo(a=[]):\n    pass\n", tmp_path)
        assert any(i.code == "MUTABLE_DEFAULT" for i in issues)

    def test_mutable_default_dict(self, tmp_path):
        issues = self._run("def foo(a={}):\n    pass\n", tmp_path)
        assert any(i.code == "MUTABLE_DEFAULT" for i in issues)

    def test_open_no_encoding(self, tmp_path):
        issues = self._run("with open('x') as f:\n    pass\n", tmp_path)
        assert any(i.code == "OPEN_NO_ENCODING" for i in issues)

    def test_open_with_encoding_ok(self, tmp_path):
        issues = self._run("with open('x', encoding='utf-8') as f:\n    pass\n", tmp_path)
        assert not any(i.code == "OPEN_NO_ENCODING" for i in issues)

    def test_eval_exec(self, tmp_path):
        issues = self._run("eval('1+1')\n", tmp_path)
        assert any(i.code == "UNSAFE_EVAL" for i in issues)

    def test_is_literal(self, tmp_path):
        issues = self._run('x = "a" is "b"\n', tmp_path)
        assert any(i.code == "IS_LITERAL" for i in issues)

    def test_is_none_ok(self, tmp_path):
        issues = self._run("x is None\n", tmp_path)
        assert not any(i.code == "IS_LITERAL" for i in issues)

    def test_duplicate_set_item(self, tmp_path):
        issues = self._run("x = {1, 2, 1}\n", tmp_path)
        assert any(i.code == "DUPLICATE_SET_ITEM" for i in issues)

    def test_duplicate_dict_key(self, tmp_path):
        issues = self._run("x = {'a': 1, 'a': 2}\n", tmp_path)
        assert any(i.code == "DUPLICATE_DICT_KEY" for i in issues)

    def test_undefined_name(self, tmp_path):
        issues = self._run("print(undefined_thing)\n", tmp_path)
        assert any(i.code == "UNDEFINED_NAME" for i in issues)

    def test_defined_name_ok(self, tmp_path):
        issues = self._run("y = 1\nprint(y)\n", tmp_path)
        assert not any(i.code == "UNDEFINED_NAME" for i in issues)


class TestReviewGateCriticals:
    def _run(self, source: str, tmp_path):
        f = tmp_path / "test.py"
        f.write_text(source)
        gate = GPAReviewGate([f])
        gate._check_file(f)
        return gate.issues

    def test_logging_fstring(self, tmp_path):
        issues = self._run("logger.info(f'hello {x}')\n", tmp_path)
        assert any(i.code == "LOGGING_FSTRING" for i in issues)

    def test_unused_import(self, tmp_path):
        # pickle is NOT in the builtins whitelist
        issues = self._run("import pickle\n", tmp_path)
        assert any(i.code == "UNUSED_IMPORT" for i in issues)

    def test_used_import_ok(self, tmp_path):
        issues = self._run("import pickle\nprint(pickle)\n", tmp_path)
        assert not any(i.code == "UNUSED_IMPORT" for i in issues)


class TestReviewGateNits:
    def _run(self, source: str, tmp_path):
        f = tmp_path / "test.py"
        f.write_text(source)
        gate = GPAReviewGate([f])
        gate._check_file(f)
        return gate.issues

    def test_pointless_fstring(self, tmp_path):
        issues = self._run("x = f'hello'\n", tmp_path)
        assert any(i.code == "POINTLESS_FSTRING" for i in issues)


class TestReviewGateReport:
    def test_no_issues_returns_0(self, tmp_path):
        f = tmp_path / "empty.py"
        f.write_text("x = 1\n")
        gate = GPAReviewGate([f])
        assert gate.run() == 0

    def test_blockers_returns_1(self):
        gate = GPAReviewGate([Path("test.py")])
        gate.issues.append(Issue("test.py", 1, Severity.BLOCKER, "X", "msg"))
        assert gate._report() == 1

    def test_critical_returns_2(self):
        gate = GPAReviewGate([Path("test.py")])
        gate.issues.append(Issue("test.py", 1, Severity.CRITICAL, "X", "msg"))
        assert gate._report() == 2

    def test_nits_returns_3(self):
        gate = GPAReviewGate([Path("test.py")])
        gate.issues.append(Issue("test.py", 1, Severity.NIT, "X", "msg"))
        assert gate._report() == 3


# =============================================================================
# dgra_myvariant
# =============================================================================

class TestBuildVariantId:
    def test_basic(self):
        assert _build_variant_id("chr1", 12345, "A", "G") == "chr1:g.12345A>G"

    def test_no_chr_prefix(self):
        assert _build_variant_id("1", 12345, "A", "G") == "chr1:g.12345A>G"

    def test_chrX(self):
        assert _build_variant_id("chrX", 12345, "C", "T") == "chrX:g.12345C>T"


class TestParseGnomadFromMyvariant:
    def test_no_gnomad(self):
        assert _parse_gnomad_from_myvariant({}) == (None, None, None, None)

    def test_genome_data(self):
        raw = {"gnomad_genome": {"af": {"af": 0.01, "ac": 10, "an": 1000}}}
        af, ac, an, pops = _parse_gnomad_from_myvariant(raw)
        assert af == 0.01
        assert ac == 10
        assert an == 1000

    def test_population_frequencies(self):
        raw = {
            "gnomad_genome": {
                "ac": {"ac_eas": 5, "ac_nfe": 10},
                "an": {"an_eas": 100, "an_nfe": 200},
            }
        }
        _, _, _, pops = _parse_gnomad_from_myvariant(raw)
        assert pops["EAS"]["af"] == 0.05
        assert pops["NFE"]["af"] == 0.05


class TestParseClinvarFromMyvariant:
    def test_no_clinvar(self):
        assert _parse_clinvar_from_myvariant({}) == (None, None, None)

    def test_rcv_list(self):
        raw = {"clinvar": {"rcv": [{"clinical_significance": "Pathogenic", "review_status": "criteria provided", "accession": "RCV001"}]}}
        sig, review, acc = _parse_clinvar_from_myvariant(raw)
        assert sig == "Pathogenic"
        assert review == "criteria provided"
        assert acc == "RCV001"

    def test_direct_fields(self):
        raw = {"clinvar": {"clinical_significance": "Likely_pathogenic", "review_status": "single submitter"}}
        sig, review, acc = _parse_clinvar_from_myvariant(raw)
        assert sig == "Likely_pathogenic"


class TestParseCaddFromMyvariant:
    def test_no_cadd(self):
        assert _parse_cadd_from_myvariant({}) is None

    def test_valid(self):
        assert _parse_cadd_from_myvariant({"cadd": {"phred": 25.3}}) == 25.3


class TestParseDbsnpFromMyvariant:
    def test_no_dbsnp(self):
        assert _parse_dbsnp_from_myvariant({}) is None

    def test_int_rsid(self):
        assert _parse_dbsnp_from_myvariant({"dbsnp": {"rsid": 12345}}) == "rs12345"

    def test_str_rsid(self):
        assert _parse_dbsnp_from_myvariant({"dbsnp": {"rsid": "rs12345"}}) == "rs12345"


class TestApplyMyvariantResults:
    def test_fill_gnomad(self):
        from dgra_core import Variant
        v = Variant(chrom="chr1", pos=12345, ref="A", alt="G", gene="X",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="", gnomad_af=None)
        result = MyVariantResult(
            variant_id="chr1:g.12345A>G",
            queried=True, status="success",
            gnomad_af=0.01, clinvar_significance="Pathogenic", cadd_phred=20.0,
        )
        stats = apply_myvariant_results([v], {"chr1:g.12345A>G": result})
        assert v.gnomad_af == 0.01
        assert v.clinvar == "Pathogenic"
        assert v.vcf_info.get("cadd_phred") == 20.0
        assert stats["gnomad_filled"] == 1
        assert stats["clinvar_filled"] == 1
        assert stats["cadd_filled"] == 1

    def test_not_found(self):
        from dgra_core import Variant
        v = Variant(chrom="chr1", pos=12345, ref="A", alt="G", gene="X",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        result = MyVariantResult(variant_id="chr1:g.12345A>G", queried=True, status="not_found")
        stats = apply_myvariant_results([v], {"chr1:g.12345A>G": result})
        assert stats["not_found"] == 1

    def test_no_overwrite_existing(self):
        from dgra_core import Variant
        v = Variant(chrom="chr1", pos=12345, ref="A", alt="G", gene="X",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="Pathogenic", gnomad_af=0.5)
        result = MyVariantResult(
            variant_id="chr1:g.12345A>G",
            queried=True, status="success",
            gnomad_af=0.01, clinvar_significance="Benign",
        )
        stats = apply_myvariant_results([v], {"chr1:g.12345A>G": result})
        assert v.gnomad_af == 0.5  # not overwritten
        assert v.clinvar == "Pathogenic"  # not overwritten
        assert stats["gnomad_filled"] == 0
        assert stats["clinvar_filled"] == 0


# =============================================================================
# gpa_proxy_routes
# =============================================================================

class _MockResponse:
    def __init__(self, status, json_data=None, text="err"):
        self.status = status
        self._json = json_data or {}
        self._text = text

    async def json(self):
        return self._json

    async def text(self):
        return self._text

    async def __aenter__(self):
        return self

    async def __aexit__(self, *args):
        pass


class _MockSession:
    def __init__(self, resp=None, raise_on_get=None):
        self._resp = resp
        self._raise_on_get = raise_on_get

    async def __aenter__(self):
        return self

    async def __aexit__(self, *args):
        pass

    def get(self, *args, **kwargs):
        if self._raise_on_get:
            raise self._raise_on_get
        return self._resp

    async def close(self):
        pass


def _make_mock_session_for_probe(status=200, json_data=None, side_effect=None):
    """Return a _MockSession that acts as both ClientSession and session."""
    resp = _MockResponse(status, json_data)
    return _MockSession(resp=resp, raise_on_get=side_effect)


class TestProxyRouteMap:
    def test_get_proxy_pass(self):
        route = ProxyRoute(api_name="ensembl", best_proxy="http://p1", latency_ms=100, status="PASS")
        pmap = ProxyRouteMap(routes={"ensembl": route})
        assert pmap.get_proxy("ensembl") == "http://p1"

    def test_get_proxy_fail(self):
        route = ProxyRoute(api_name="ensembl", best_proxy="http://p1", latency_ms=100, status="FAIL")
        pmap = ProxyRouteMap(routes={"ensembl": route})
        assert pmap.get_proxy("ensembl") is None

    def test_get_fallback(self):
        route = ProxyRoute(
            api_name="ensembl", best_proxy="http://p1", latency_ms=100, status="PASS",
            all_results=[{"proxy": "http://p1", "status": "PASS"}, {"proxy": "http://p2", "status": "PASS"}]
        )
        pmap = ProxyRouteMap(routes={"ensembl": route})
        assert pmap.get_fallback("ensembl", exclude="http://p1") == "http://p2"

    def test_to_dict(self):
        route = ProxyRoute(api_name="ensembl", best_proxy=None, latency_ms=50, status="PASS")
        pmap = ProxyRouteMap(routes={"ensembl": route})
        d = pmap.to_dict()
        assert d["routes"]["ensembl"]["best_proxy"] is None

    def test_to_markdown(self):
        route = ProxyRoute(api_name="ensembl", best_proxy=None, latency_ms=50, status="PASS")
        pmap = ProxyRouteMap(routes={"ensembl": route})
        md = pmap.to_markdown()
        assert "ensembl" in md
        assert "直连" in md


class TestBuildCandidateProxies:
    def test_includes_direct(self):
        proxies = _build_candidate_proxies()
        assert None in proxies


class TestProbeSingle:
    @pytest.mark.asyncio
    async def test_aiohttp_missing(self):
        with patch("gpa_proxy_routes.aiohttp", None):
            result = await _probe_single("http://test", None)
            assert result["status"] == "FAIL"
            assert "aiohttp not installed" in result["error"]

    @pytest.mark.asyncio
    async def test_success(self):
        mock_client = _make_mock_session_for_probe(status=200, json_data={"ok": True})
        with patch("gpa_proxy_routes.aiohttp.ClientSession", return_value=mock_client):
            result = await _probe_single("http://test", None)
            assert result["status"] == "PASS"

    @pytest.mark.asyncio
    async def test_success_with_validator(self):
        mock_client = _make_mock_session_for_probe(status=200, json_data={"ok": True})
        with patch("gpa_proxy_routes.aiohttp.ClientSession", return_value=mock_client):
            result = await _probe_single("http://test", None, validator=lambda d: d.get("ok"))
            assert result["status"] == "PASS"

    @pytest.mark.asyncio
    async def test_validator_fail(self):
        mock_client = _make_mock_session_for_probe(status=200, json_data={"ok": False})
        with patch("gpa_proxy_routes.aiohttp.ClientSession", return_value=mock_client):
            result = await _probe_single("http://test", None, validator=lambda d: d.get("ok"))
            assert result["status"] == "WARN"

    @pytest.mark.asyncio
    async def test_spliceai_400_pass(self):
        mock_client = _make_mock_session_for_probe(status=400)
        with patch("gpa_proxy_routes.aiohttp.ClientSession", return_value=mock_client):
            result = await _probe_single("http://spliceai/api", None)
            assert result["status"] == "PASS"

    @pytest.mark.asyncio
    async def test_non_200_fail(self):
        mock_client = _make_mock_session_for_probe(status=500)
        with patch("gpa_proxy_routes.aiohttp.ClientSession", return_value=mock_client):
            result = await _probe_single("http://test", None)
            assert result["status"] == "FAIL"

    @pytest.mark.asyncio
    async def test_timeout(self):
        mock_client = _make_mock_session_for_probe(side_effect=asyncio.TimeoutError())
        with patch("gpa_proxy_routes.aiohttp.ClientSession", return_value=mock_client):
            result = await _probe_single("http://test", None)
            assert result["status"] == "FAIL"
            assert "Timeout" in result["error"]


class TestProbeApiRoutes:
    @pytest.mark.asyncio
    async def test_pass_result(self):
        with patch("gpa_proxy_routes._probe_single", return_value={"proxy": None, "status": "PASS", "latency_ms": 100}):
            route = await probe_api_routes("test", "http://test")
            assert route.status == "PASS"
            assert route.best_proxy is None

    @pytest.mark.asyncio
    async def test_warn_fallback(self):
        with patch("gpa_proxy_routes._probe_single", return_value={"proxy": None, "status": "WARN", "latency_ms": 100}):
            route = await probe_api_routes("test", "http://test")
            assert route.status == "WARN"

    @pytest.mark.asyncio
    async def test_all_fail(self):
        with patch("gpa_proxy_routes._probe_single", return_value={"proxy": None, "status": "FAIL", "latency_ms": 999}):
            route = await probe_api_routes("test", "http://test", proxies=[None])
            assert route.status == "FAIL"


class TestBuildRouteMap:
    @pytest.mark.asyncio
    async def test_build(self):
        checks = {"api1": ("http://test1", 5.0, None)}
        with patch("gpa_proxy_routes.probe_api_routes", return_value=ProxyRoute(api_name="api1", best_proxy=None, latency_ms=50, status="PASS")):
            pmap = await build_route_map(checks)
            assert "api1" in pmap.routes
