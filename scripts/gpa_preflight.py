#!/usr/bin/env python3
"""
GPA Preflight Health Check — v0.10.1

每次接到全新分析任务时，先执行一次可用性检查，确认所有依赖就绪后再
启动耗时较长的分析流程。避免分析到一半才发现 API 不通、工具缺失或
磁盘空间不足。

检查范围（六大类）：
  1. Python 依赖包    — aiohttp, vcfpy, yaml
  2. 本地命令行工具   — vep (可选), git (可选)
  3. 在线 API 连通性   — Ensembl, UniProt, GTEx, gnomAD, NCBI, HGNC,
                         MyVariant.info, SpliceAI (可选)
  4. 本地文件/目录     — references/, cache/, 配置文件
  5. 磁盘空间         — cache 目录所在分区
  6. 网络/代理环境     — 直连 vs 代理可用性

用法（在 pipeline 入口调用）：
    from gpa_preflight import run_preflight_check, interactive_prompt

    report = await run_preflight_check(config)
    if not report.is_ready():
        action = interactive_prompt(report)
        if action == "abort":
            return {"error": "Preflight check failed", "report": report.to_dict()}
        elif action == "offline":
            config.offline_mode = True
        # "continue" → 继续执行，缺失的可选功能会被跳过
"""

import asyncio
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from gpa_proxy_routes import ProxyRouteMap, build_route_map

# =============================================================================
# 0. 对现有模块的弱依赖（导入失败时不崩溃，降级为 None）
# =============================================================================
try:
    from dgra_config import DGRAGlobalConfig
except Exception:
    DGRAGlobalConfig = None  # type: ignore[misc,assignment]

try:
    import aiohttp
except Exception:
    aiohttp = None  # type: ignore[misc,assignment]

# =============================================================================
# 1. 数据模型
# =============================================================================

CHECK_CATEGORIES = [
    "python_deps",
    "local_tools",
    "api_connectivity",
    "local_files",
    "disk_space",
    "network_env",
]


@dataclass
class CheckItem:
    """单个检查项的结果。"""

    name: str
    category: str  # 必须属于 CHECK_CATEGORIES
    required: bool  # True = blocker if FAIL, False = warning only
    status: str  # PASS | FAIL | WARN | SKIP
    message: str
    latency_ms: Optional[int] = None
    suggestion: Optional[str] = None


@dataclass
class PreflightReport:
    """完整的预检报告。"""

    items: List[CheckItem] = field(default_factory=list)
    timestamp: float = field(default_factory=time.time)

    # -------------------------------------------------------------------------
    # 查询方法
    # -------------------------------------------------------------------------
    def is_ready(self) -> bool:
        """所有 required 项均为 PASS → 可直接进入分析。"""
        return all(i.status == "PASS" for i in self.items if i.required)

    def blockers(self) -> List[CheckItem]:
        """返回所有 required 且未 PASS 的项（必须解决或切换到离线模式）。"""
        return [i for i in self.items if i.required and i.status != "PASS"]

    def warnings(self) -> List[CheckItem]:
        """返回所有 optional 且未 PASS 的项（可跳过，但功能受限）。"""
        return [i for i in self.items if not i.required and i.status != "PASS"]

    def by_category(self, category: str) -> List[CheckItem]:
        return [i for i in self.items if i.category == category]

    # -------------------------------------------------------------------------
    # 输出方法
    # -------------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return {
            "timestamp": self.timestamp,
            "overall_ready": self.is_ready(),
            "blocker_count": len(self.blockers()),
            "warning_count": len(self.warnings()),
            "items": [
                {
                    "name": i.name,
                    "category": i.category,
                    "required": i.required,
                    "status": i.status,
                    "message": i.message,
                    "latency_ms": i.latency_ms,
                    "suggestion": i.suggestion,
                }
                for i in self.items
            ],
        }

    def to_markdown(self) -> str:
        """生成人类可读的 Markdown 报告（供 CLI 或 Agent 展示）。"""
        lines: List[str] = []
        lines.append("## GPA 前置可用性检查报告\n")

        # 总体结论
        if self.is_ready():
            lines.append("**总体状态**: 全部就绪，可直接开始分析。\n")
        elif self.blockers():
            lines.append(
                f"**总体状态**: 存在 {len(self.blockers())} 项必须修复的 blocker，"
                f"{len(self.warnings())} 项警告。\n"
            )
        else:
            lines.append(
                f"**总体状态**: 无 blocker，但存在 {len(self.warnings())} 项可选功能不可用。\n"
            )

        # 按分类输出
        for cat in CHECK_CATEGORIES:
            items = self.by_category(cat)
            if not items:
                continue
            cat_display = {
                "python_deps": "Python 依赖",
                "local_tools": "本地命令行工具",
                "api_connectivity": "在线 API 连通性",
                "local_files": "本地文件/目录",
                "disk_space": "磁盘空间",
                "network_env": "网络/代理环境",
            }.get(cat, cat)
            lines.append(f"\n### {cat_display}\n")
            for i in items:
                icon = {"PASS": "✅", "FAIL": "❌", "WARN": "⚠️", "SKIP": "⏭️"}.get(
                    i.status, "❓"
                )
                req = "[必须]" if i.required else "[可选]"
                lines.append(f"- {icon} {req} **{i.name}**: {i.message}")
                if i.latency_ms is not None and i.status == "PASS":
                    lines.append(f"  - 延迟: {i.latency_ms} ms")
                if i.suggestion:
                    lines.append(f"  - 建议: {i.suggestion}")
        return "\n".join(lines)


# =============================================================================
# 2. 检查实现 — Python 依赖
# =============================================================================

def _check_python_package(name: str, import_name: Optional[str] = None) -> CheckItem:
    """检查单个 Python 包是否可导入。"""
    import_name = import_name or name
    spec = importlib.util.find_spec(import_name)
    if spec is not None:
        return CheckItem(
            name=name,
            category="python_deps",
            required=True,
            status="PASS",
            message=f"已安装 (路径: {spec.origin})",
        )
    return CheckItem(
        name=name,
        category="python_deps",
        required=True,
        status="FAIL",
        message="未安装",
        suggestion=f"pip install {name}",
    )


def check_python_deps() -> List[CheckItem]:
    """检查所有必需的 Python 第三方包。"""
    return [
        _check_python_package("aiohttp"),
        _check_python_package("vcfpy"),
        _check_python_package("yaml", "yaml"),  # PyYAML
    ]


# =============================================================================
# 3. 检查实现 — 本地命令行工具
# =============================================================================

def _check_cli_tool(name: str, required: bool = False) -> CheckItem:
    """检查单个命令行工具是否在 PATH 中。"""
    path = shutil.which(name)
    if path:
        return CheckItem(
            name=name,
            category="local_tools",
            required=required,
            status="PASS",
            message=f"已找到: {path}",
        )
    return CheckItem(
        name=name,
        category="local_tools",
        required=required,
        status="FAIL" if required else "WARN",
        message="未在 PATH 中找到",
        suggestion=f"安装 {name} 并将其加入 PATH",
    )


def check_local_tools() -> List[CheckItem]:
    """检查本地命令行工具。"""
    return [
        # vep: 可选 — 只有使用 --annotator vep_local 时才需要
        _check_cli_tool("vep", required=False),
        # git: 可选 — 仅用于版本信息收集
        _check_cli_tool("git", required=False),
    ]


# =============================================================================
# 4. 检查实现 — 在线 API 连通性
# =============================================================================

_API_CHECKS: Dict[str, Tuple[str, int, Optional[Callable[[Any], bool]]]] = {
    # name → (probe_url, timeout_sec, response_validator)
    "ensembl": (
        "https://rest.ensembl.org/info/ping?content-type=application/json",
        5,
        lambda d: isinstance(d, dict) and "ping" in d,
    ),
    "uniprot": (
        "https://rest.uniprot.org/uniprotkb/search?query=insulin&size=1",
        5,
        lambda d: isinstance(d, dict) and "results" in d,
    ),
    "gtex": (
        "https://gtexportal.org/api/v2/reference/gene?geneId=ENSG00000141510",
        5,
        lambda d: isinstance(d, dict) and "data" in d,
    ),
    "gnomad": (
        "https://gnomad.broadinstitute.org/api",
        5,
        None,  # 只要 HTTP 200 就算通
    ),
    "ncbi_eutils": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=clinvar&retmode=json",
        5,
        lambda d: isinstance(d, dict) and "header" in d,
    ),
    "hgnc": (
        "https://rest.genenames.org/search/status/Approved",
        5,
        lambda d: isinstance(d, dict) and "response" in d,
    ),
    "myvariant": (
        "https://myvariant.info/v1/metadata",
        5,
        lambda d: isinstance(d, dict) and "stats" in d,
    ),
    # SpliceAI: Broad 托管的 Cloud Run 服务。 probe 一个假 variant，
    # 预期 400/422（REF mismatch），但只要能连上就视为服务在线。
    "spliceai": (
        "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant=chr1-12345-A-T",
        5,
        None,  # 任何 HTTP 响应（包括 400）都算在线
    ),
}


async def _probe_api(
    name: str,
    url: str,
    timeout: float = 5.0,
    proxy: Optional[str] = None,
    validator: Optional[Callable[[Any], bool]] = None,
) -> CheckItem:
    """对单个 API 端点发起轻量探测请求。"""
    if aiohttp is None:
        return CheckItem(
            name=name,
            category="api_connectivity",
            required=True,
            status="SKIP",
            message="aiohttp 未安装，跳过 API 检查",
            suggestion="pip install aiohttp",
        )

    t0 = time.perf_counter()
    session: Optional[Any] = None
    try:
        session = aiohttp.ClientSession(
            trust_env=False,
            timeout=aiohttp.ClientTimeout(total=timeout),
        )
        async with session.get(url, proxy=proxy, allow_redirects=True) as resp:
            latency_ms = int((time.perf_counter() - t0) * 1000)
            # SpliceAI 特殊处理：400/422 也算在线
            if name == "spliceai" and resp.status in (200, 400, 422):
                return CheckItem(
                    name=name,
                    category="api_connectivity",
                    required=False,  # SpliceAI 有 VEP REST fallback，标记为可选
                    status="PASS",
                    message=f"服务在线 (HTTP {resp.status})",
                    latency_ms=latency_ms,
                )
            if resp.status != 200:
                return CheckItem(
                    name=name,
                    category="api_connectivity",
                    required=True,
                    status="FAIL",
                    message=f"HTTP {resp.status}",
                    latency_ms=latency_ms,
                    suggestion="检查网络或代理配置",
                )
            # 尝试解析 JSON 并运行 validator
            try:
                data = await resp.json()
            except Exception:
                data = None
            if validator is not None and not validator(data):
                return CheckItem(
                    name=name,
                    category="api_connectivity",
                    required=True,
                    status="WARN",
                    message=f"HTTP 200 但响应格式异常",
                    latency_ms=latency_ms,
                    suggestion="API 可能正在维护，稍后重试",
                )
            return CheckItem(
                name=name,
                category="api_connectivity",
                required=True,
                status="PASS",
                message="连通正常",
                latency_ms=latency_ms,
            )
    except asyncio.TimeoutError:
        latency_ms = int((time.perf_counter() - t0) * 1000)
        return CheckItem(
            name=name,
            category="api_connectivity",
            required=True,
            status="FAIL",
            message=f"超时（>{timeout}s）",
            latency_ms=latency_ms,
            suggestion="检查网络连接或增大 timeout",
        )
    except Exception as e:
        latency_ms = int((time.perf_counter() - t0) * 1000)
        return CheckItem(
            name=name,
            category="api_connectivity",
            required=True,
            status="FAIL",
            message=f"连接异常: {type(e).__name__}: {str(e)[:60]}",
            latency_ms=latency_ms,
            suggestion="检查网络、代理或防火墙设置",
        )
    finally:
        if session is not None:
            await session.close()


async def check_api_connectivity(
    config: Optional[Any] = None,
) -> List[CheckItem]:
    """并发检查所有在线 API 的连通性。"""
    # 从 config 中提取代理设置（如果有的话）
    proxy: Optional[str] = None
    if config is not None and hasattr(config, "apis"):
        # 优先使用 ensembl 的 proxy，回退到全局 proxy
        proxy = getattr(config.apis.get("ensembl"), "proxy", None)
        if not proxy and hasattr(config, "proxy"):
            proxy = config.proxy
        if proxy == "__DIRECT__":
            proxy = None

    # 并发执行所有探测
    tasks = [
        _probe_api(name, url, timeout, proxy, validator)
        for name, (url, timeout, validator) in _API_CHECKS.items()
    ]
    return await asyncio.gather(*tasks)


# =============================================================================
# 5. 检查实现 — 本地文件/目录
# =============================================================================

def _resolve_path_from_config(
    config: Optional[Any], attr_name: str, default: Path
) -> Path:
    """从 config 安全地解析路径，失败时回退到默认值。"""
    if config is None:
        return default
    val = getattr(config, attr_name, None)
    if val is None:
        return default
    if isinstance(val, Path):
        return val
    try:
        return Path(val)
    except Exception:
        return default


def check_local_files(config: Optional[Any] = None) -> List[CheckItem]:
    """检查必需的本地文件和目录。"""
    script_dir = Path(__file__).resolve().parent
    default_refs = script_dir.parent / "references"
    default_cache = script_dir.parent / "cache"

    refs_dir = _resolve_path_from_config(
        config, "override_files", default_refs  # config 中没有直接的 refs_dir
    )
    # 如果 config 的 override_files 是 dict，回退到默认
    if isinstance(refs_dir, dict):
        refs_dir = default_refs

    cache_db = _resolve_path_from_config(
        config, "cache_db_path", default_cache / "dgra_cache.db"
    )

    items: List[CheckItem] = []

    # references 目录
    if refs_dir.exists() and refs_dir.is_dir():
        items.append(
            CheckItem(
                name="references/ 目录",
                category="local_files",
                required=True,
                status="PASS",
                message=f"存在: {refs_dir}",
            )
        )
    else:
        items.append(
            CheckItem(
                name="references/ 目录",
                category="local_files",
                required=True,
                status="FAIL",
                message=f"缺失: {refs_dir}",
                suggestion="创建 references/ 目录并放置必要的 JSON 配置文件",
            )
        )

    # tissue_context.json（Tier 分类必需）
    tissue_json = refs_dir / "tissue_context.json"
    if tissue_json.exists():
        items.append(
            CheckItem(
                name="tissue_context.json",
                category="local_files",
                required=True,
                status="PASS",
                message=f"存在: {tissue_json}",
            )
        )
    else:
        items.append(
            CheckItem(
                name="tissue_context.json",
                category="local_files",
                required=True,
                status="FAIL",
                message=f"缺失: {tissue_json}",
                suggestion="该文件为 Tier 分类的必需输入",
            )
        )

    # cache 目录
    cache_dir = cache_db.parent if hasattr(cache_db, "parent") else default_cache
    if cache_dir.exists() and cache_dir.is_dir():
        items.append(
            CheckItem(
                name="cache/ 目录",
                category="local_files",
                required=True,
                status="PASS",
                message=f"存在: {cache_dir}",
            )
        )
    else:
        # 自动创建（非 blocker）
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            items.append(
                CheckItem(
                    name="cache/ 目录",
                    category="local_files",
                    required=True,
                    status="PASS",
                    message=f"已自动创建: {cache_dir}",
                )
            )
        except Exception as e:
            items.append(
                CheckItem(
                    name="cache/ 目录",
                    category="local_files",
                    required=True,
                    status="FAIL",
                    message=f"无法创建: {e}",
                    suggestion="检查目录写入权限",
                )
            )

    # dgra_cache.db（可选 — 首次运行时不存在，会自动创建）
    if cache_db.exists():
        items.append(
            CheckItem(
                name="dgra_cache.db",
                category="local_files",
                required=False,
                status="PASS",
                message=f"存在: {cache_db}",
            )
        )
    else:
        items.append(
            CheckItem(
                name="dgra_cache.db",
                category="local_files",
                required=False,
                status="WARN",
                message="不存在（首次运行会自动创建）",
                suggestion="无需操作",
            )
        )

    # dgra.yaml（可选配置）
    yaml_path = refs_dir / "dgra.yaml"
    if yaml_path.exists():
        items.append(
            CheckItem(
                name="dgra.yaml",
                category="local_files",
                required=False,
                status="PASS",
                message=f"存在: {yaml_path}",
            )
        )
    else:
        items.append(
            CheckItem(
                name="dgra.yaml",
                category="local_files",
                required=False,
                status="WARN",
                message="不存在（使用默认配置）",
                suggestion="如需自定义参数，可创建 references/dgra.yaml",
            )
        )

    return items


# =============================================================================
# 6. 检查实现 — 磁盘空间
# =============================================================================

def check_disk_space(min_mb: int = 500) -> List[CheckItem]:
    """检查 cache 目录所在分区的可用空间。"""
    script_dir = Path(__file__).resolve().parent
    check_path = script_dir.parent / "cache"
    try:
        check_path.mkdir(parents=True, exist_ok=True)
    except Exception:
        check_path = Path.home()

    try:
        usage = shutil.disk_usage(check_path)
        free_mb = usage.free // (1024 * 1024)
        if free_mb >= min_mb:
            return [
                CheckItem(
                    name="磁盘可用空间",
                    category="disk_space",
                    required=True,
                    status="PASS",
                    message=f"{free_mb} MB 可用（阈值: {min_mb} MB）",
                )
            ]
        return [
            CheckItem(
                name="磁盘可用空间",
                category="disk_space",
                required=True,
                status="FAIL",
                message=f"仅 {free_mb} MB 可用（低于阈值 {min_mb} MB）",
                suggestion="清理 cache/ 目录或扩容磁盘",
            )
        ]
    except Exception as e:
        return [
            CheckItem(
                name="磁盘可用空间",
                category="disk_space",
                required=True,
                status="FAIL",
                message=f"无法获取磁盘信息: {e}",
                suggestion="检查文件系统权限",
            )
        ]


# =============================================================================
# 7. 检查实现 — 网络/代理环境
# =============================================================================

async def check_network_proxy() -> List[CheckItem]:
    """探测网络环境：直连是否可用，常用代理是否可用。"""
    if aiohttp is None:
        return [
            CheckItem(
                name="网络环境探测",
                category="network_env",
                required=False,
                status="SKIP",
                message="aiohttp 未安装，跳过",
            )
        ]

    test_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        "?db=clinvar&term=BRCA1%5BGene%5D&retmode=json&retmax=1"
    )
    items: List[CheckItem] = []

    # 1. 直连测试
    t0 = time.perf_counter()
    session: Optional[Any] = None
    try:
        session = aiohttp.ClientSession(
            trust_env=False,
            timeout=aiohttp.ClientTimeout(total=5),
        )
        async with session.get(test_url) as resp:
            latency_ms = int((time.perf_counter() - t0) * 1000)
            if resp.status == 200:
                data = await resp.json()
                count = data.get("esearchresult", {}).get("count")
                if count and int(count) > 0:
                    items.append(
                        CheckItem(
                            name="直连 NCBI",
                            category="network_env",
                            required=False,
                            status="PASS",
                            message=f"直连可用，延迟 {latency_ms} ms",
                            latency_ms=latency_ms,
                        )
                    )
                else:
                    items.append(
                        CheckItem(
                            name="直连 NCBI",
                            category="network_env",
                            required=False,
                            status="WARN",
                            message="HTTP 200 但响应异常",
                            latency_ms=latency_ms,
                        )
                    )
            else:
                items.append(
                    CheckItem(
                        name="直连 NCBI",
                        category="network_env",
                        required=False,
                        status="FAIL",
                        message=f"HTTP {resp.status}",
                        latency_ms=latency_ms,
                        suggestion="可能需要配置代理",
                    )
                )
    except asyncio.TimeoutError:
        latency_ms = int((time.perf_counter() - t0) * 1000)
        items.append(
            CheckItem(
                name="直连 NCBI",
                category="network_env",
                required=False,
                status="FAIL",
                message=f"超时（>5s）",
                latency_ms=latency_ms,
                suggestion="检查网络连接，或尝试配置代理",
            )
        )
    except Exception as e:
        latency_ms = int((time.perf_counter() - t0) * 1000)
        items.append(
            CheckItem(
                name="直连 NCBI",
                category="network_env",
                required=False,
                status="FAIL",
                message=f"异常: {type(e).__name__}: {str(e)[:60]}",
                latency_ms=latency_ms,
                suggestion="检查网络连接或代理配置",
            )
        )
    finally:
        if session is not None:
            await session.close()

    # 2. 常用代理测试（仅当直连失败时才测试，节省时间）
    if not any(i.status == "PASS" for i in items if i.name == "直连 NCBI"):
        common_proxies = [
            "http://127.0.0.1:7897",
            "http://127.0.0.1:7890",
            "http://127.0.0.1:7891",
            "http://127.0.0.1:1080",
            "http://127.0.0.1:10808",
            "http://127.0.0.1:10809",
        ]
        # 并发测试，最多 3 个并发
        semaphore = asyncio.Semaphore(3)

        async def _test_proxy(proxy: str) -> CheckItem:
            t0 = time.perf_counter()
            sess: Optional[Any] = None
            try:
                async with semaphore:
                    sess = aiohttp.ClientSession(
                        trust_env=False,
                        timeout=aiohttp.ClientTimeout(total=5),
                    )
                    async with sess.get(test_url, proxy=proxy) as resp:
                        latency_ms = int((time.perf_counter() - t0) * 1000)
                        if resp.status == 200:
                            return CheckItem(
                                name=f"代理 {proxy}",
                                category="network_env",
                                required=False,
                                status="PASS",
                                message=f"代理可用，延迟 {latency_ms} ms",
                                latency_ms=latency_ms,
                                suggestion=f"可在 dgra_config 中设置 proxy='{proxy}'",
                            )
                        return CheckItem(
                            name=f"代理 {proxy}",
                            category="network_env",
                            required=False,
                            status="FAIL",
                            message=f"HTTP {resp.status}",
                            latency_ms=latency_ms,
                        )
            except Exception as e:
                latency_ms = int((time.perf_counter() - t0) * 1000)
                return CheckItem(
                    name=f"代理 {proxy}",
                    category="network_env",
                    required=False,
                    status="FAIL",
                    message=f"异常: {type(e).__name__}",
                    latency_ms=latency_ms,
                )
            finally:
                if sess is not None:
                    await sess.close()

        proxy_tasks = [asyncio.create_task(_test_proxy(p)) for p in common_proxies]
        proxy_results = await asyncio.gather(*proxy_tasks)
        # 只保留第一个通过的代理，其余省略以减少噪音
        passed = [r for r in proxy_results if r.status == "PASS"]
        if passed:
            items.append(passed[0])
            if len(passed) > 1:
                items.append(
                    CheckItem(
                        name="其他可用代理",
                        category="network_env",
                        required=False,
                        status="PASS",
                        message=f"另有 {len(passed) - 1} 个代理可用",
                    )
                )
        else:
            items.append(
                CheckItem(
                    name="代理探测",
                    category="network_env",
                    required=False,
                    status="FAIL",
                    message="未找到可用代理",
                    suggestion="手动配置 HTTP_PROXY / HTTPS_PROXY 环境变量",
                )
            )

    return items


# =============================================================================
# 8. 主入口
# =============================================================================

async def run_preflight_check(
    config: Optional[Any] = None,
    check_categories: Optional[List[str]] = None,
    skip_api_if_offline: bool = True,
) -> Tuple[PreflightReport, ProxyRouteMap]:
    """执行完整的前置可用性检查。

    Args:
        config: DGRAGlobalConfig 实例。为 None 时使用默认路径。
        check_categories: 指定只检查哪些分类。None = 全部。
        skip_api_if_offline: 如果 config.offline_mode=True，跳过 API 检查。

    Returns:
        (PreflightReport, ProxyRouteMap) 元组。
        ProxyRouteMap 包含每个 API 的最佳代理路由（即使 API 检查被跳过也会返回空表）。
    """
    if check_categories is None:
        check_categories = CHECK_CATEGORIES[:]

    report = PreflightReport()
    route_map = ProxyRouteMap()

    # 1. Python 依赖
    if "python_deps" in check_categories:
        report.items.extend(check_python_deps())

    # 2. 本地工具
    if "local_tools" in check_categories:
        report.items.extend(check_local_tools())

    # 3. 本地文件
    if "local_files" in check_categories:
        report.items.extend(check_local_files(config))

    # 4. 磁盘空间
    if "disk_space" in check_categories:
        report.items.extend(check_disk_space())

    # 5. 网络/代理（在 API 之前，结果可用于 API 检查）
    if "network_env" in check_categories:
        report.items.extend(await check_network_proxy())

    # 6. API 连通性 + 代理路由表构建
    if "api_connectivity" in check_categories:
        is_offline = False
        if config is not None and hasattr(config, "offline_mode"):
            is_offline = bool(getattr(config, "offline_mode"))
        if is_offline and skip_api_if_offline:
            report.items.append(
                CheckItem(
                    name="API 连通性",
                    category="api_connectivity",
                    required=False,
                    status="SKIP",
                    message="offline_mode=True，跳过所有 API 检查",
                )
            )
        else:
            # v0.10.12: 使用 per-API 多代理探测，生成代理路由表
            print("[Preflight] 正在探测每个 API 的最佳代理路由...")
            route_map = await build_route_map(_API_CHECKS)
            # 将路由探测结果转换为 CheckItem 并入报告
            for name, route in route_map.routes.items():
                proxy_str = route.best_proxy or "直连"
                if route.status == "PASS":
                    report.items.append(
                        CheckItem(
                            name=name,
                            category="api_connectivity",
                            required=True,
                            status="PASS",
                            message=f"通过 {proxy_str} 连通 ({route.latency_ms} ms)",
                            latency_ms=route.latency_ms,
                        )
                    )
                else:
                    report.items.append(
                        CheckItem(
                            name=name,
                            category="api_connectivity",
                            required=True,
                            status="FAIL",
                            message="所有代理路由均失败",
                            suggestion="检查网络连接或代理配置",
                        )
                    )

    return report, route_map


# =============================================================================
# 9. 交互式提示
# =============================================================================

def interactive_prompt(report: PreflightReport) -> str:
    """根据预检报告生成用户交互提示文本，返回建议动作。

    返回值（供调用方判断）：
        - "continue": 无 blocker，或有 blocker 但用户选择继续（进入离线模式）
        - "offline": 用户明确选择切换到离线模式
        - "abort": 用户选择中止任务
    """
    blockers = report.blockers()
    warnings = report.warnings()

    if report.is_ready():
        return "continue"

    lines: List[str] = []
    lines.append("=" * 60)
    lines.append("⚠️  GPA 前置可用性检查未完全通过")
    lines.append("=" * 60)

    if blockers:
        lines.append(f"\n必须修复的 blocker（{len(blockers)} 项）：")
        for i in blockers:
            lines.append(f"  ❌ [{i.category}] {i.name}: {i.message}")
            if i.suggestion:
                lines.append(f"     → {i.suggestion}")

    if warnings:
        lines.append(f"\n可选功能缺失（{len(warnings)} 项）：")
        for i in warnings:
            lines.append(f"  ⚠️  [{i.category}] {i.name}: {i.message}")

    lines.append("\n可选操作：")
    if blockers:
        lines.append("  [1] 切换到离线模式 — 跳过所有 API 调用，仅使用本地缓存")
        lines.append("     （适合已有缓存或仅需本地 Tier 分类的场景）")
    lines.append("  [2] 忽略警告继续 — 缺失的功能会在运行时自动跳过")
    lines.append("  [3] 中止任务 — 修复环境后重试")

    return "\n".join(lines)


def suggest_action(report: PreflightReport) -> str:
    """根据报告自动建议最佳动作（无用户交互时由调用方使用）。

    Returns:
        "continue" | "abort"
    """
    CRITICAL_APIS = {"ensembl"}  # VEP annotation is essential; gnomAD has MyVariant.info fallback
    blockers = report.blockers()
    # 如果无 blocker，只有 warnings → 建议继续
    if not blockers:
        return "continue"
    # 如果 blocker 包含关键 API → 必须中止（离线模式会导致大量变异无法评估）
    if any(i.name in CRITICAL_APIS for i in blockers):
        return "abort"
    # 如果 blocker 包含 python_deps / local_files / disk_space → 必须中止
    if any(i.category in ("python_deps", "local_files", "disk_space") for i in blockers):
        return "abort"
    # 剩下的 blocker 都是非关键 API 连通性失败 → 继续（功能降级但不阻断核心分析）
    return "continue"


# =============================================================================
# 10. CLI 自测入口
# =============================================================================

async def _main() -> None:
    """命令行自测：python scripts/gpa_preflight.py"""
    print("[GPA Preflight] 正在执行可用性检查...\n")
    report = await run_preflight_check()
    print(report.to_markdown())
    print("\n" + "=" * 60)
    if report.is_ready():
        print("结论: 全部就绪 ✅")
    else:
        action = suggest_action(report)
        print(f"建议动作: {action}")
        print("\n" + interactive_prompt(report))


if __name__ == "__main__":
    asyncio.run(_main())
