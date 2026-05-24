"""
dgra_splice_predictor.py — SpliceAI 剪接预测集成（GPA v0.8.0）

对剪接相关变异自动查询 SpliceAI 分数，作为剪接功能影响的独立证据。
默认关闭，需 --spliceai 显式开启。

外部依赖：aiohttp
"""

import asyncio
import aiohttp
import logging
from typing import Optional, Dict, Any, List
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# SpliceAI lookup API (Broad Institute)
# v0.9.4: URL migrated from spliceailookup.broadinstitute.org/api/variant (404)
# to Google Cloud Run endpoints (verified 2026-05-24 via GitHub SpliceAI-lookup README + live test)
# GRCh38: https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant=chr-pos-ref-alt
# GRCh37: https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=37&variant=chr-pos-ref-alt
# Rate limit: interactive use only (~few req/min). For bulk, run local Docker instance.
SPLICEAI_BASE_URL_GRCh38 = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/"
SPLICEAI_BASE_URL_GRCh37 = "https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/"

# 阈值配置（分 consequence 类型）
SPLICEAI_THRESHOLDS = {
    "canonical": {  # splice_acceptor_variant / splice_donor_variant
        "strong": 0.5,
        "moderate": 0.2,
        "weak": 0.1,
        "none": 0.0,
    },
    "splice_region": {  # splice_region_variant / splice_polypyrimidine / splice_donor_5th_base
        "strong": 0.2,
        "moderate": 0.1,
        "weak": 0.05,
        "none": 0.0,
    },
}

# 触发查询的 consequence 类型（英文 SO terms）
SPLICE_QUERY_TERMS = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_polypyrimidine_tract_variant",
    "splice_donor_5th_base_variant",
    # synonymous_variant 单独判断（仅 ±50bp 内外显子边界时查询）
}

# 中文映射（与 gpa_i18n.py 同步）
_CN_SPLICE_TERMS = {
    "剪接受体位点变异", "剪接供体位点变异", "剪接区域变异",
    "剪接多嘧啶束变异", "剪接供体第5位碱基变异",
}


@dataclass
class SpliceAIResult:
    """SpliceAI 查询结果。"""
    delta_score: float = 0.0          # max of (AG, AL, DG, DL)
    delta_acceptor_gain: float = 0.0
    delta_acceptor_loss: float = 0.0
    delta_donor_gain: float = 0.0
    delta_donor_loss: float = 0.0
    predicted_impact: str = "none"  # strong/moderate/weak/none/unknown
    threshold_type: str = "canonical" # canonical / splice_region
    source: str = "unknown"           # spliceai_lookup / not_in_db / api_error
    raw_response: Optional[Dict] = None


class SpliceAIPredictor:
    """SpliceAI 预测器，支持异步批量查询、缓存、退避重试。"""

    def __init__(self, max_concurrency: int = 5, timeout: int = 30):
        self.max_concurrency = max_concurrency
        self.timeout = timeout
        self._cache: Dict[str, SpliceAIResult] = {}
        self._semaphore = asyncio.Semaphore(max_concurrency)

    def _cache_key(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """生成缓存键。"""
        return f"{chrom}:{pos}:{ref}:{alt}"

    @staticmethod
    def should_query(consequence_terms: List[str], is_near_exon_boundary: bool = False) -> bool:
        """判断该变异是否需要 SpliceAI 评估。

        Args:
            consequence_terms: 英文 SO terms 列表（已 normalize）
            is_near_exon_boundary: 同义变异专用，是否在 ±50bp 内外显子边界

        Returns:
            True 表示需要查询 SpliceAI
        """
        for term in consequence_terms:
            if term in SPLICE_QUERY_TERMS:
                return True
            # 中文术语检查（fallback）
            if term in _CN_SPLICE_TERMS:
                return True
            # 同义变异仅在靠近外显子边界时查询
            if term == "synonymous_variant" and is_near_exon_boundary:
                return True
        return False

    @staticmethod
    def is_canonical_splice(consequence_terms: List[str]) -> bool:
        """判断是否为 canonical splice（±1,2 bp 的 acceptor/donor）。"""
        canonical = {"splice_acceptor_variant", "splice_donor_variant",
                     "剪接受体位点变异", "剪接供体位点变异"}
        return any(t in canonical for t in consequence_terms)

    @staticmethod
    def determine_impact(delta_score: float, threshold_type: str) -> str:
        """根据 delta score 和阈值类型判断影响等级。"""
        th = SPLICEAI_THRESHOLDS.get(threshold_type, SPLICEAI_THRESHOLDS["splice_region"])
        if delta_score >= th["strong"]:
            return "strong"
        elif delta_score >= th["moderate"]:
            return "moderate"
        elif delta_score >= th["weak"]:
            return "weak"
        else:
            return "none"

    async def query(self, chrom: str, pos: int, ref: str, alt: str, genome: str = "GRCh38") -> SpliceAIResult:
        """查询 SpliceAI lookup，带缓存和退避重试。

        Args:
            chrom: 染色体（"1", "chr1" 均可，内部标准化）
            pos: 1-based 坐标
            ref: 参考碱基
            alt: 变异碱基
            genome: 基因组版本，默认 GRCh38

        Returns:
            SpliceAIResult 对象
        """
        # 标准化染色体名（去掉 chr 前缀）
        chrom_std = chrom.replace("chr", "") if chrom.startswith("chr") else chrom

        cache_key = self._cache_key(chrom_std, pos, ref, alt)
        if cache_key in self._cache:
            logger.debug(f"SpliceAI cache hit: {cache_key}")
            return self._cache[cache_key]

        result = await self._query_with_retry(chrom_std, pos, ref, alt, genome)
        self._cache[cache_key] = result
        return result

    async def _query_with_retry(self, chrom: str, pos: int, ref: str, alt: str, genome: str) -> SpliceAIResult:
        """带退避重试的 SpliceAI 查询。
        
        v0.9.4: 适配新 API 格式（Google Cloud Run 端点）。
        新格式：GET /spliceai/?hg=38&variant=chr-pos-ref-alt
        响应格式：{"scores": [{ "DS_AG": "0.13", ... }], "source": "..."}
        """
        # 选择对应基因组版本的 URL
        if "37" in genome:
            base_url = SPLICEAI_BASE_URL_GRCh37
        else:
            base_url = SPLICEAI_BASE_URL_GRCh38
        
        # 新 API 格式：variant 参数为 chr-pos-ref-alt
        variant_str = f"{chrom}-{pos}-{ref}-{alt}"
        params = {
            "hg": genome.replace("GRCh", ""),   # GRCh38 → "38"
            "variant": variant_str,
        }

        # 指数退避：2s → 4s → 8s（新 API 有限流，延迟更长）
        backoff = [2, 4, 8]

        async with self._semaphore:
            for attempt, delay in enumerate(backoff + [None]):
                try:
                    timeout = aiohttp.ClientTimeout(total=self.timeout)
                    async with aiohttp.ClientSession(timeout=timeout, trust_env=False) as session:
                        async with session.get(base_url, params=params) as resp:
                            if resp.status == 200:
                                data = await resp.json()
                                return self._parse_response(data, "spliceai_lookup")
                            elif resp.status == 404:
                                # 不在数据库（某些变异可能返回 404 或空 scores）
                                return SpliceAIResult(source="not_in_db")
                            elif resp.status == 429:
                                retry_after = int(resp.headers.get("Retry-After", delay or 5))
                                logger.warning(f"SpliceAI 429, retry after {retry_after}s (attempt {attempt+1})")
                                if delay is not None:
                                    await asyncio.sleep(retry_after)
                                    continue
                            else:
                                logger.warning(f"SpliceAI HTTP {resp.status}, attempt {attempt+1}")
                                if delay is not None:
                                    await asyncio.sleep(delay)
                                    continue
                except asyncio.TimeoutError:
                    logger.warning(f"SpliceAI timeout, attempt {attempt+1}")
                    if delay is not None:
                        await asyncio.sleep(delay)
                        continue
                except aiohttp.ClientError as e:
                    logger.warning(f"SpliceAI client error: {e}, attempt {attempt+1}")
                    if delay is not None:
                        await asyncio.sleep(delay)
                        continue

            # 所有重试失败
            logger.error(f"SpliceAI query failed after {len(backoff)} retries: {chrom}:{pos}:{ref}:{alt}")
            return SpliceAIResult(source="api_error")

    def _parse_response(self, data: Dict, source: str) -> SpliceAIResult:
        """解析 SpliceAI lookup API 响应。
        
        v0.9.4: 适配新 API 格式（Google Cloud Run 端点）。
        新格式: {"scores": [{"DS_AG": "0.13", "DS_AL": "0.00", ...}, ...], "source": "..."}
        旧格式（已废弃）: {"delta_scores": {"AG": 0.0, "AL": 0.0, ...}}
        """
        try:
            scores_list = data.get("scores", [])
            
            if not scores_list:
                # 无 scores = 无剪接影响
                return SpliceAIResult(
                    delta_score=0.0,
                    source=source,
                    raw_response=data,
                )

            # 取所有 transcript 中 delta score 最大值（保守策略）
            max_delta = 0.0
            max_ag = 0.0
            max_al = 0.0
            max_dg = 0.0
            max_dl = 0.0

            for score_entry in scores_list:
                # 新 API 返回字符串格式的浮点数
                ag = float(score_entry.get("DS_AG", 0.0) or 0.0)
                al = float(score_entry.get("DS_AL", 0.0) or 0.0)
                dg = float(score_entry.get("DS_DG", 0.0) or 0.0)
                dl = float(score_entry.get("DS_DL", 0.0) or 0.0)
                entry_max = max(ag, al, dg, dl)

                if entry_max > max_delta:
                    max_delta = entry_max
                    max_ag = ag
                    max_al = al
                    max_dg = dg
                    max_dl = dl

            return SpliceAIResult(
                delta_score=max_delta,
                delta_acceptor_gain=max_ag,
                delta_acceptor_loss=max_al,
                delta_donor_gain=max_dg,
                delta_donor_loss=max_dl,
                predicted_impact="unknown",  # 外部根据 threshold_type 确定
                source=source,
                raw_response=data,
            )
        except Exception as e:
            logger.error(f"SpliceAI response parse error: {e}, data={data}")
            return SpliceAIResult(source="api_error")

    async def batch_query(self, variants: List[Any]) -> Dict[str, SpliceAIResult]:
        """批量查询多个变异的 SpliceAI 分数。

        Args:
            variants: 列表，每个元素需有 chrom, pos, ref, alt, consequence 属性

        Returns:
            Dict[cache_key, SpliceAIResult]
        """
        tasks = []
        keys = []
        for v in variants:
            cache_key = self._cache_key(v.chrom, v.pos, v.ref, v.alt)
            if cache_key in self._cache:
                continue
            tasks.append(self.query(v.chrom, v.pos, v.ref, v.alt))
            keys.append(cache_key)

        if not tasks:
            return dict(self._cache)

        results = await asyncio.gather(*tasks, return_exceptions=True)
        for key, res in zip(keys, results):
            if isinstance(res, Exception):
                logger.error(f"SpliceAI batch query error for {key}: {res}")
                self._cache[key] = SpliceAIResult(source="api_error")
            else:
                self._cache[key] = res

        return dict(self._cache)


def should_query_spliceai(consequence_terms, is_near_exon_boundary: bool = False) -> bool:
    """同步判断是否需要 SpliceAI 查询。

    兼容传入字符串（单条 consequence）或列表（多条 consequence）。
    """
    if isinstance(consequence_terms, str):
        from gpa_i18n import normalize_consequence
        consequence_terms = normalize_consequence(consequence_terms)
        if not consequence_terms:
            return False
    return SpliceAIPredictor.should_query(consequence_terms, is_near_exon_boundary)


def _cache_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """模块级缓存键生成（供外部使用）。"""
    chrom_std = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
    return f"{chrom_std}:{pos}:{ref}:{alt}"


# 全局 predictor 实例（延迟初始化）
_spliceai_predictor: Optional[SpliceAIPredictor] = None


def _get_predictor(max_concurrency: int = 5, timeout: int = 30) -> SpliceAIPredictor:
    """获取/创建全局 SpliceAI predictor 实例。"""
    global _spliceai_predictor
    if _spliceai_predictor is None:
        _spliceai_predictor = SpliceAIPredictor(max_concurrency=max_concurrency, timeout=timeout)
    return _spliceai_predictor


def reset_spliceai_cache() -> None:
    """重置全局 SpliceAI 缓存。"""
    global _spliceai_predictor
    if _spliceai_predictor is not None:
        _spliceai_predictor._cache.clear()
        _spliceai_predictor = None  # 也销毁实例，下次重新创建


async def query_spliceai_batch(
    variants: List[Any],
    session: Any,
    semaphore: Any,
    max_concurrency: int = 5,
    timeout: int = 30,
) -> Dict[str, SpliceAIResult]:
    """批量查询 SpliceAI，供 dgra_core.py 调用。

    Args:
        variants: 需要查询的变异列表（需有 chrom, pos, ref, alt 属性）
        session: aiohttp ClientSession（由调用方管理生命周期）
        semaphore: asyncio.Semaphore（由调用方管理并发）
        max_concurrency: 最大并发数
        timeout: 请求超时（秒）

    Returns:
        Dict[cache_key, SpliceAIResult]
    """
    predictor = _get_predictor(max_concurrency=max_concurrency, timeout=timeout)
    results: Dict[str, SpliceAIResult] = {}

    for v in variants:
        cache_key = _cache_key(v.chrom, v.pos, v.ref, v.alt)
        if cache_key in predictor._cache:
            results[cache_key] = predictor._cache[cache_key]
            continue

        async with semaphore:
            try:
                result = await predictor.query(v.chrom, v.pos, v.ref, v.alt)
                results[cache_key] = result
                predictor._cache[cache_key] = result
            except Exception as e:
                logger.error(f"SpliceAI batch query error for {cache_key}: {e}")
                err_result = SpliceAIResult(source="api_error")
                results[cache_key] = err_result
                predictor._cache[cache_key] = err_result

    return results


# 保持向后兼容的便捷函数
