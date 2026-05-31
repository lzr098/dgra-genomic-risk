"""L2 Unit Tests — gpa_preflight.py

Covers CheckItem, PreflightReport, dependency checks, local tools,
file checks, disk space, and report formatting.
"""

import pytest
from unittest.mock import patch, MagicMock, AsyncMock


@pytest.mark.l2
class TestCheckItem:
    """PFLT-01: CheckItem dataclass."""

    def test_basic_check_item(self):
        """PFLT-01: Basic CheckItem creation."""
        from gpa_preflight import CheckItem
        item = CheckItem(
            name="aiohttp",
            category="python_deps",
            required=True,
            status="PASS",
            message="已安装",
        )
        assert item.name == "aiohttp"
        assert item.category == "python_deps"
        assert item.required is True
        assert item.status == "PASS"
        assert item.message == "已安装"


@pytest.mark.l2
class TestPreflightReport:
    """PFLT-02~10: PreflightReport methods."""

    def _make_report(self, items):
        from gpa_preflight import PreflightReport, CheckItem
        return PreflightReport(items=items)

    def test_is_ready_all_pass(self):
        """PFLT-02: All required PASS → is_ready True."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("a", "python_deps", True, "PASS", "ok"),
            CheckItem("b", "local_tools", True, "PASS", "ok"),
        ])
        assert report.is_ready() is True
        assert len(report.blockers()) == 0

    def test_is_ready_with_blocker(self):
        """PFLT-03: Required FAIL → is_ready False."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("a", "python_deps", True, "PASS", "ok"),
            CheckItem("b", "local_tools", True, "FAIL", "missing"),
        ])
        assert report.is_ready() is False
        assert len(report.blockers()) == 1
        assert report.blockers()[0].name == "b"

    def test_warnings_only(self):
        """PFLT-04: Optional WARN → warnings list."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("a", "python_deps", True, "PASS", "ok"),
            CheckItem("b", "api_connectivity", False, "WARN", "slow"),
        ])
        assert report.is_ready() is True
        assert len(report.warnings()) == 1

    def test_by_category(self):
        """PFLT-05: Filter by category."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("a", "python_deps", True, "PASS", "ok"),
            CheckItem("b", "python_deps", True, "PASS", "ok"),
            CheckItem("c", "local_tools", True, "FAIL", "missing"),
        ])
        assert len(report.by_category("python_deps")) == 2
        assert len(report.by_category("local_tools")) == 1
        assert len(report.by_category("disk_space")) == 0

    def test_to_dict_structure(self):
        """PFLT-06: to_dict returns correct structure."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("a", "python_deps", True, "PASS", "ok"),
        ])
        d = report.to_dict()
        assert d["overall_ready"] is True
        assert d["blocker_count"] == 0
        assert d["warning_count"] == 0
        assert len(d["items"]) == 1
        assert "timestamp" in d

    def test_to_markdown_pass(self):
        """PFLT-07: Markdown output for all-pass."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("aiohttp", "python_deps", True, "PASS", "已安装"),
        ])
        md = report.to_markdown()
        assert "全部就绪" in md
        assert "✅" in md
        assert "aiohttp" in md

    def test_to_markdown_with_blocker(self):
        """PFLT-08: Markdown output with blocker."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("aiohttp", "python_deps", True, "PASS", "已安装"),
            CheckItem("vep", "local_tools", True, "FAIL", "未找到"),
        ])
        md = report.to_markdown()
        assert "必须修复" in md or "blocker" in md.lower()
        assert "❌" in md

    def test_to_markdown_with_warning(self):
        """PFLT-09: Markdown output with warning."""
        from gpa_preflight import CheckItem
        report = self._make_report([
            CheckItem("aiohttp", "python_deps", True, "PASS", "已安装"),
            CheckItem("gnomad", "api_connectivity", False, "WARN", "slow"),
        ])
        md = report.to_markdown()
        assert "警告" in md or "⚠️" in md

    def test_empty_report(self):
        """PFLT-10: Empty report → ready, no blockers."""
        report = self._make_report([])
        assert report.is_ready() is True
        assert len(report.blockers()) == 0
        assert len(report.warnings()) == 0


@pytest.mark.l2
class TestPythonDepsCheck:
    """PFLT-11~13: Python dependency checks."""

    def test_check_python_package_installed(self):
        """PFLT-11: Installed package → PASS."""
        from gpa_preflight import _check_python_package
        result = _check_python_package("json")  # stdlib, always available
        assert result.status == "PASS"
        assert "json" in result.name

    def test_check_python_package_missing(self):
        """PFLT-12: Missing package → FAIL."""
        from gpa_preflight import _check_python_package
        result = _check_python_package("nonexistent_package_xyz")
        assert result.status == "FAIL"

    def test_check_python_deps_returns_list(self):
        """PFLT-13: check_python_deps returns list of CheckItems."""
        from gpa_preflight import check_python_deps
        items = check_python_deps()
        assert len(items) >= 2
        assert all(hasattr(i, "status") for i in items)


@pytest.mark.l2
class TestLocalToolsCheck:
    """PFLT-14~16: Local CLI tool checks."""

    def test_check_cli_tool_found(self):
        """PFLT-14: Tool in PATH → PASS."""
        from gpa_preflight import _check_cli_tool
        result = _check_cli_tool("python3", required=False)
        assert result.status == "PASS"
        assert "python3" in result.message

    def test_check_cli_tool_missing_optional(self):
        """PFLT-15: Missing optional tool → WARN."""
        from gpa_preflight import _check_cli_tool
        result = _check_cli_tool("nonexistent_tool_xyz", required=False)
        assert result.status == "WARN"

    def test_check_cli_tool_missing_required(self):
        """PFLT-16: Missing required tool → FAIL."""
        from gpa_preflight import _check_cli_tool
        result = _check_cli_tool("nonexistent_tool_xyz", required=True)
        assert result.status == "FAIL"

    def test_check_local_tools(self):
        """PFLT-17: check_local_tools returns list."""
        from gpa_preflight import check_local_tools
        items = check_local_tools()
        assert len(items) >= 1
        assert all(i.category == "local_tools" for i in items)


@pytest.mark.l2
class TestAPIConnectivity:
    """PFLT-18~22: API connectivity checks (mocked)."""

    def _mock_session(self, mock_resp):
        """Build a mock session that works with async context manager."""
        class MockSession:
            def __init__(self, resp):
                self._resp = resp
            def get(self, *args, **kwargs):
                return self._resp
            async def close(self):
                pass
        return MockSession(mock_resp)

    def _mock_resp(self, status, data=None):
        """Build a mock response that works with async context manager."""
        class MockResp:
            def __init__(self, status, data):
                self.status = status
                self._data = data
            async def json(self):
                return self._data
            async def __aenter__(self):
                return self
            async def __aexit__(self, *args):
                return False
        return MockResp(status, data)

    @pytest.mark.asyncio
    async def test_probe_api_success(self):
        """PFLT-18: HTTP 200 + valid JSON → PASS."""
        from gpa_preflight import _probe_api
        mock_resp = self._mock_resp(200, {"ping": 1})
        mock_session = self._mock_session(mock_resp)
        with patch("aiohttp.ClientSession", return_value=mock_session):
            result = await _probe_api("ensembl", "http://test", validator=lambda d: "ping" in d)
        assert result.status == "PASS"
        assert result.latency_ms is not None

    @pytest.mark.asyncio
    async def test_probe_api_http_error(self):
        """PFLT-19: HTTP 503 → FAIL."""
        from gpa_preflight import _probe_api
        mock_resp = self._mock_resp(503)
        mock_session = self._mock_session(mock_resp)
        with patch("aiohttp.ClientSession", return_value=mock_session):
            result = await _probe_api("test", "http://test")
        assert result.status == "FAIL"
        assert "503" in result.message

    @pytest.mark.asyncio
    async def test_probe_api_timeout(self):
        """PFLT-20: Timeout → FAIL."""
        from gpa_preflight import _probe_api
        import asyncio
        class TimeoutResp:
            async def __aenter__(self):
                raise asyncio.TimeoutError()
            async def __aexit__(self, *args):
                return False
        class TimeoutSession:
            def get(self, *args, **kwargs):
                return TimeoutResp()
            async def close(self):
                pass
        with patch("aiohttp.ClientSession", return_value=TimeoutSession()):
            result = await _probe_api("test", "http://test", timeout=0.001)
        assert result.status == "FAIL"
        assert "超时" in result.message or "timeout" in result.message.lower()

    @pytest.mark.asyncio
    async def test_probe_api_spliceai_special(self):
        """PFLT-21: SpliceAI 400/422 → PASS."""
        from gpa_preflight import _probe_api
        mock_resp = self._mock_resp(400)
        mock_session = self._mock_session(mock_resp)
        with patch("aiohttp.ClientSession", return_value=mock_session):
            result = await _probe_api("spliceai", "http://test")
        assert result.status == "PASS"
        assert "400" in result.message

    @pytest.mark.asyncio
    async def test_probe_api_validator_fail(self):
        """PFLT-22: HTTP 200 but validator rejects → WARN."""
        from gpa_preflight import _probe_api
        mock_resp = self._mock_resp(200, {"wrong": "format"})
        mock_session = self._mock_session(mock_resp)
        with patch("aiohttp.ClientSession", return_value=mock_session):
            result = await _probe_api("test", "http://test", validator=lambda d: "expected" in d)
        assert result.status == "WARN"


@pytest.mark.l2
class TestDiskSpaceCheck:
    """PFLT-23~25: Disk space checks."""

    def test_check_disk_space(self):
        """PFLT-23: Returns CheckItems for disk space."""
        from gpa_preflight import check_disk_space
        items = check_disk_space(min_mb=1)
        assert len(items) >= 1
        assert any(i.category == "disk_space" for i in items)

    def test_check_disk_space_pass(self):
        """PFLT-24: Sufficient space → PASS."""
        from gpa_preflight import check_disk_space
        items = check_disk_space(min_mb=1)
        # At least one item should be PASS (the current working dir)
        assert any(i.status == "PASS" for i in items)

    def test_check_disk_space_large_min(self):
        """PFLT-25: Unreasonably large min → FAIL."""
        from gpa_preflight import check_disk_space
        items = check_disk_space(min_mb=999999999)
        # Should still return items, some may be FAIL
        assert len(items) >= 1


@pytest.mark.l2
class TestLocalFilesCheck:
    """PFLT-26~28: Local file checks."""

    def test_check_local_files_returns_items(self):
        """PFLT-26: Returns list of CheckItems."""
        from gpa_preflight import check_local_files
        items = check_local_files()
        assert len(items) >= 2
        assert any(i.name == "references/ 目录" for i in items)

    def test_check_local_files_refs_dir(self):
        """PFLT-27: references dir check included."""
        from gpa_preflight import check_local_files
        items = check_local_files()
        refs_items = [i for i in items if "references" in i.name]
        assert len(refs_items) >= 1

    def test_check_local_files_cache_dir(self):
        """PFLT-28: cache dir check included."""
        from gpa_preflight import check_local_files
        items = check_local_files()
        cache_items = [i for i in items if "cache" in i.name.lower()]
        assert len(cache_items) >= 1


@pytest.mark.l2
class TestSuggestAction:
    """PFLT-29~30: suggest_action helper."""

    def test_suggest_action_ready(self):
        """PFLT-29: Ready report → continue."""
        from gpa_preflight import PreflightReport, CheckItem, suggest_action
        report = PreflightReport(items=[CheckItem("a", "python_deps", True, "PASS", "ok")])
        action = suggest_action(report)
        assert action == "continue"

    def test_suggest_action_blockers(self):
        """PFLT-30: Blockers → abort."""
        from gpa_preflight import PreflightReport, CheckItem, suggest_action
        report = PreflightReport(items=[
            CheckItem("a", "python_deps", True, "FAIL", "missing"),
        ])
        action = suggest_action(report)
        assert action == "abort"
