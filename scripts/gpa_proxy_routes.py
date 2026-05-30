#!/usr/bin/env python3
"""
GPA Proxy Route Map — Per-API proxy routing with preflight probe.

每次任务启动时，并发探测每个外部 API 通过不同代理的连通性，
为每个 API 选择延迟最低的成功路由。运行时各 API 按自己的
最佳路由独立走代理，互不影响。

v0.10.12
"""

import asyncio
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

try:
    import aiohttp
except Exception:  # noqa: BROAD_EXCEPT — graceful fallback when aiohttp is missing
    aiohttp = None  # type: ignore


# =============================================================================
# 数据模型
# =============================================================================

@dataclass
class ProxyRoute:
    """单个 API 的最佳代理路由及所有候选结果."""

    api_name: str
    best_proxy: Optional[str]  # None = direct
    latency_ms: int
    status: str  # PASS | FAIL
    all_results: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class ProxyRouteMap:
    """完整的代理路由表，按 API 名索引."""

    routes: Dict[str, ProxyRoute] = field(default_factory=dict)
    detected_at: float = field(default_factory=time.time)

    def get_proxy(self, api_name: str) -> Optional[str]:
        """获取指定 API 的最佳代理（None = 直连）."""
        route = self.routes.get(api_name)
        if route and route.status == "PASS":
            return route.best_proxy
        return None

    def get_fallback(self, api_name: str, exclude: Optional[str] = None) -> Optional[str]:
        """获取备选代理（排除已失败的主代理）."""
        route = self.routes.get(api_name)
        if not route:
            return None
        for r in route.all_results:
            if r["status"] == "PASS" and r.get("proxy") != exclude:
                return r.get("proxy")
        return None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "detected_at": self.detected_at,
            "routes": {
                name: {
                    "best_proxy": route.best_proxy,
                    "latency_ms": route.latency_ms,
                    "status": route.status,
                }
                for name, route in self.routes.items()
            },
        }

    def to_markdown(self) -> str:
        lines = ["## 代理路由表\n"]
        for name, route in sorted(self.routes.items()):
            proxy_str = route.best_proxy or "直连"
            icon = "✅" if route.status == "PASS" else "❌"
            lines.append(f"- {icon} **{name}**: {proxy_str} ({route.latency_ms} ms)")
        return "\n".join(lines)


# =============================================================================
# 候选代理列表（按优先级排序）
# =============================================================================

# =============================================================================
# 候选代理列表（按优先级排序）
# =============================================================================

def _build_candidate_proxies() -> List[Optional[str]]:
    """Build proxy candidate list including system env proxies."""
    candidates: List[Optional[str]] = [None]  # direct first
    # Add common local proxy ports
    candidates.extend([
        "http://127.0.0.1:7897",
        "http://127.0.0.1:7890",
        "http://127.0.0.1:7891",
        "http://127.0.0.1:1080",
        "http://127.0.0.1:10808",
        "http://127.0.0.1:10809",
        "http://127.0.0.1:52402",  # sandbox proxy
    ])
    # v0.10.12-fix: Also probe system env proxies (HTTP_PROXY, HTTPS_PROXY)
    # because curl defaults to them, making "direct" tests actually proxied.
    env_proxies = []
    for key in ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy"):
        val = os.environ.get(key)
        if val and val not in env_proxies:
            env_proxies.append(val)
    for ep in env_proxies:
        if ep not in candidates:
            candidates.append(ep)
    return candidates


CANDIDATE_PROXIES: List[Optional[str]] = _build_candidate_proxies()


# =============================================================================
# 探测实现
# =============================================================================

async def _probe_single(
    url: str,
    proxy: Optional[str],
    timeout: float = 5.0,
    validator: Optional[Callable[[Any], bool]] = None,
) -> Dict[str, Any]:
    """对单个 URL + 代理发起探测请求."""
    if aiohttp is None:
        return {
            "proxy": proxy,
            "status": "FAIL",
            "latency_ms": 999999,
            "error": "aiohttp not installed",
        }

    t0 = time.perf_counter()
    session = None
    try:
        session = aiohttp.ClientSession(
            trust_env=False,
            timeout=aiohttp.ClientTimeout(total=timeout),
        )
        async with session.get(url, proxy=proxy, allow_redirects=True) as resp:
            latency_ms = int((time.perf_counter() - t0) * 1000)

            # SpliceAI 特殊处理：400/422 也算在线
            is_spliceai = "spliceai" in url.lower()
            if is_spliceai and resp.status in (200, 400, 422):
                return {
                    "proxy": proxy,
                    "status": "PASS",
                    "latency_ms": latency_ms,
                    "http_status": resp.status,
                }

            if resp.status != 200:
                return {
                    "proxy": proxy,
                    "status": "FAIL",
                    "latency_ms": latency_ms,
                    "http_status": resp.status,
                    "error": f"HTTP {resp.status}",
                }

            try:
                data = await resp.json()
            except Exception:  # noqa: BROAD_EXCEPT — JSON parse failure is non-fatal for probe
                data = None

            if validator is not None and not validator(data):
                return {
                    "proxy": proxy,
                    "status": "WARN",
                    "latency_ms": latency_ms,
                    "error": "Response validation failed",
                }

            return {"proxy": proxy, "status": "PASS", "latency_ms": latency_ms}
    except asyncio.TimeoutError:
        latency_ms = int((time.perf_counter() - t0) * 1000)
        return {
            "proxy": proxy,
            "status": "FAIL",
            "latency_ms": latency_ms,
            "error": "Timeout",
        }
    except Exception as e:  # noqa: BROAD_EXCEPT — probe must survive any network error
        latency_ms = int((time.perf_counter() - t0) * 1000)
        return {
            "proxy": proxy,
            "status": "FAIL",
            "latency_ms": latency_ms,
            "error": f"{type(e).__name__}: {str(e)[:60]}",
        }
    finally:
        if session is not None:
            await session.close()


async def probe_api_routes(
    name: str,
    url: str,
    timeout: float = 5.0,
    validator: Optional[Callable[[Any], bool]] = None,
    proxies: Optional[List[Optional[str]]] = None,
) -> ProxyRoute:
    """并发探测单个 API 的所有候选代理，返回最佳路由."""
    proxies = proxies or CANDIDATE_PROXIES

    # 并发探测所有候选代理
    tasks = [_probe_single(url, proxy, timeout, validator) for proxy in proxies]
    results = await asyncio.gather(*tasks)

    # 选择最佳路由（延迟最低的 PASS；若无 PASS，则选延迟最低的 WARN）
    passed = [r for r in results if r["status"] == "PASS"]
    if passed:
        best = min(passed, key=lambda r: r["latency_ms"])
        return ProxyRoute(
            api_name=name,
            best_proxy=best["proxy"],
            latency_ms=best["latency_ms"],
            status="PASS",
            all_results=list(results),
        )

    # v0.10.12-fix: WARN means network reachable (e.g., HTTP 400 from VEP probe)
    # Treat best WARN as usable route rather than failing entirely.
    warned = [r for r in results if r["status"] == "WARN"]
    if warned:
        best = min(warned, key=lambda r: r["latency_ms"])
        return ProxyRoute(
            api_name=name,
            best_proxy=best["proxy"],
            latency_ms=best["latency_ms"],
            status="WARN",
            all_results=list(results),
        )

    return ProxyRoute(
        api_name=name,
        best_proxy=None,
        latency_ms=999999,
        status="FAIL",
        all_results=list(results),
    )


async def build_route_map(
    api_checks: Dict[str, Tuple[str, float, Optional[Callable[[Any], bool]]]],
    proxies: Optional[List[Optional[str]]] = None,
) -> ProxyRouteMap:
    """为所有 API 构建完整的路由表."""
    routes: Dict[str, ProxyRoute] = {}

    # 并发探测所有 API
    tasks = [
        probe_api_routes(name, url, timeout, validator, proxies)
        for name, (url, timeout, validator) in api_checks.items()
    ]
    results = await asyncio.gather(*tasks)

    for route in results:
        routes[route.api_name] = route

    return ProxyRouteMap(routes=routes)


# =============================================================================
# CLI 自测入口
# =============================================================================

async def _main() -> None:
    """命令行自测: python scripts/gpa_proxy_routes.py"""
    from gpa_preflight import _API_CHECKS

    print("[ProxyRouteMap] 正在探测所有 API 的代理路由...\n")
    route_map = await build_route_map(_API_CHECKS)
    print(route_map.to_markdown())
    print("\n" + "=" * 60)
    fail_count = sum(1 for r in route_map.routes.values() if r.status != "PASS")
    if fail_count == 0:
        print("全部 API 均有可用路由 ✅")
    else:
        print(f"{fail_count} 个 API 无可用路由 ❌")


if __name__ == "__main__":
    asyncio.run(_main())
