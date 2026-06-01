"""
SpliceAI 模块 v0.11.3 优化专项测试
验证基于 GPA_SpliceAI_Module_Diagnosis.md 报告的所有修改点。
"""
import sys, os, asyncio, json, logging
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class FakeVariant:
    def __init__(self, chrom, pos, ref, alt, consequence):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.consequence = consequence


# ========================================================================
# Test 1: 默认超时 45s（P0-1）
# ========================================================================
def test_default_timeout_45s():
    from dgra_splice_predictor import SpliceAIPredictor, _get_predictor
    p1 = SpliceAIPredictor()
    assert p1.timeout == 45, f"Expected default timeout=45, got {p1.timeout}"
    p2 = _get_predictor()
    assert p2.timeout == 45, f"Expected _get_predictor default timeout=45, got {p2.timeout}"
    print("  ✅ test_default_timeout_45s")


# ========================================================================
# Test 2: 自定义超时参数传递（P1-6 CLI）
# ========================================================================
def test_custom_timeout():
    from dgra_splice_predictor import SpliceAIPredictor
    p = SpliceAIPredictor(timeout=60)
    assert p.timeout == 60, f"Expected timeout=60, got {p.timeout}"
    print("  ✅ test_custom_timeout")


# ========================================================================
# Test 3: 退避重试时间 [5, 10, 15, 20]（P1-3）
# ========================================================================
def test_backoff_intervals():
    from dgra_splice_predictor import SpliceAIPredictor
    p = SpliceAIPredictor()
    # 通过检查源码中的 backoff 变量来验证
    import inspect
    src = inspect.getsource(p._query_with_retry)
    assert "[5, 10, 15, 20]" in src, "Backoff intervals not updated to [5, 10, 15, 20]"
    print("  ✅ test_backoff_intervals")


# ========================================================================
# Test 4: VEP fallback 使用 POST 格式（P0-2）
# ========================================================================
def test_vep_fallback_post_format():
    from dgra_splice_predictor import SpliceAIPredictor
    p = SpliceAIPredictor()
    import inspect
    src = inspect.getsource(p._query_vep_rest)
    # 检查是否使用 POST 而不是 GET
    assert 'session.post(' in src, "VEP fallback should use POST method"
    # 检查是否包含 variants 参数
    assert '"variants"' in src, "VEP POST should include 'variants' in payload"
    # 检查是否包含 ref/alt 在 variant_str 中
    assert "{ref}/{alt}" in src, "VEP POST should explicitly include ref/alt"
    print("  ✅ test_vep_fallback_post_format")


# ========================================================================
# Test 5: URL 支持环境变量覆盖（P2-7）
# ========================================================================
def test_url_env_override():
    import dgra_splice_predictor as dsp
    # 环境变量未设置时，使用默认值
    assert "spliceai-38" in dsp.SPLICEAI_BASE_URL_GRCh38
    # 设置环境变量后重新导入（模拟）
    os.environ["GPA_SPLICEAI_URL_38"] = "https://custom.example.com/spliceai/"
    os.environ["GPA_SPLICEAI_URL_37"] = "https://custom37.example.com/spliceai/"
    # 重新读取模块级变量（通过重新导入或 exec）
    import importlib
    importlib.reload(dsp)
    assert dsp.SPLICEAI_BASE_URL_GRCh38 == "https://custom.example.com/spliceai/"
    assert dsp.SPLICEAI_BASE_URL_GRCh37 == "https://custom37.example.com/spliceai/"
    # 清理
    del os.environ["GPA_SPLICEAI_URL_38"]
    del os.environ["GPA_SPLICEAI_URL_37"]
    importlib.reload(dsp)
    print("  ✅ test_url_env_override")


# ========================================================================
# Test 6: 空 scores 区分 not_in_db vs no_impact（P2-6）
# ========================================================================
def test_empty_scores_distinction():
    from dgra_splice_predictor import SpliceAIPredictor, SpliceAIResult
    p = SpliceAIPredictor()
    # 模拟 "不在数据库" 响应
    not_in_db_data = {"source": "not_in_db", "scores": []}
    r1 = p._parse_response(not_in_db_data, "spliceai_lookup")
    assert r1.source == "not_in_db", f"Expected source='not_in_db', got '{r1.source}'"
    # 模拟 "无剪接影响" 响应（source 不是 not_in_db）
    no_impact_data = {"source": "spliceai_lookup", "scores": []}
    r2 = p._parse_response(no_impact_data, "spliceai_lookup")
    assert r2.source == "spliceai_lookup", f"Expected source='spliceai_lookup', got '{r2.source}'"
    assert r2.delta_score == 0.0
    print("  ✅ test_empty_scores_distinction")


# ========================================================================
# Test 7: 进度日志和异常 traceback（P1-5, P2-8）
# ========================================================================
def test_progress_and_traceback_in_source():
    from dgra_splice_predictor import SpliceAIPredictor
    p = SpliceAIPredictor()
    import inspect
    src = inspect.getsource(p.batch_query)
    # 检查进度日志
    assert "SpliceAI batch progress" in src, "batch_query should log progress"
    # 检查 traceback 导入
    assert "import traceback" in inspect.getsource(sys.modules["dgra_splice_predictor"])
    # 检查 batch_query 中使用 traceback
    assert "traceback.format_exception" in src, "batch_query should format exception traceback"
    print("  ✅ test_progress_and_traceback_in_source")


# ========================================================================
# Test 8: query_spliceai_batch 移除双重 semaphore（P1-4）
# ========================================================================
def test_no_double_semaphore():
    from dgra_splice_predictor import query_spliceai_batch
    import inspect
    src = inspect.getsource(query_spliceai_batch)
    # 不应该有外层的 async with semaphore
    assert "async with semaphore" not in src, "query_spliceai_batch should not have outer async with semaphore"
    # 但 predictor.query() 内部应该有 semaphore
    from dgra_splice_predictor import SpliceAIPredictor
    p_src = inspect.getsource(SpliceAIPredictor._query_with_retry)
    assert "async with self._semaphore" in p_src, "predictor.query should still have internal semaphore"
    print("  ✅ test_no_double_semaphore")


# ========================================================================
# Test 9: CLI wrapper 参数链 — spliceai_timeout 传递到 GPAConfig
# ========================================================================
def test_cli_wrapper_timeout_propagation():
    from dgra_cli_wrapper import _run_gpa_direct
    import inspect
    sig = inspect.signature(_run_gpa_direct)
    assert "spliceai_timeout" in sig.parameters, "_run_gpa_direct should accept spliceai_timeout"
    assert sig.parameters["spliceai_timeout"].default == 45, "Default should be 45"
    print("  ✅ test_cli_wrapper_timeout_propagation")


# ========================================================================
# Test 10: dgra_core.py GPAConfig 包含 spliceai_timeout
# ========================================================================
def test_gpa_config_has_spliceai_timeout():
    from dgra_core import GPAConfig
    import inspect
    sig = inspect.signature(GPAConfig)
    assert "spliceai_timeout" in sig.parameters, "GPAConfig should have spliceai_timeout"
    assert sig.parameters["spliceai_timeout"].default == 45, "Default should be 45"
    # 实例化验证
    cfg = GPAConfig(spliceai_timeout=60)
    assert cfg.spliceai_timeout == 60
    print("  ✅ test_gpa_config_has_spliceai_timeout")


# ========================================================================
# Test 11: batch_runner 参数传递
# ========================================================================
def test_batch_runner_timeout_propagation():
    from dgra_batch_runner import run_batch, run_gpa_batched
    import inspect
    sig1 = inspect.signature(run_batch)
    assert "spliceai_timeout" in sig1.parameters, "run_batch should accept spliceai_timeout"
    sig2 = inspect.signature(run_gpa_batched)
    assert "spliceai_timeout" in sig2.parameters, "run_gpa_batched should accept spliceai_timeout"
    print("  ✅ test_batch_runner_timeout_propagation")


# ========================================================================
# Main
# ========================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("SpliceAI v0.11.3 Optimization Tests")
    print("=" * 60)

    tests = [
        test_default_timeout_45s,
        test_custom_timeout,
        test_backoff_intervals,
        test_vep_fallback_post_format,
        test_url_env_override,
        test_empty_scores_distinction,
        test_progress_and_traceback_in_source,
        test_no_double_semaphore,
        test_cli_wrapper_timeout_propagation,
        test_gpa_config_has_spliceai_timeout,
        test_batch_runner_timeout_propagation,
    ]

    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            print(f"  ❌ {t.__name__}: {e}")
            failed += 1

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    if failed > 0:
        sys.exit(1)
