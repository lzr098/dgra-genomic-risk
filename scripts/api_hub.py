#!/usr/bin/env python3
"""
api_hub.py - Unified async HTTP layer for GPA skills.

Centralizes aiohttp ClientSession lifecycle, proxy handling, adaptive rate
limiting, retry logic, and response caching. Other modules should obtain a
session or make requests through this hub instead of creating ad-hoc
ClientSession instances.

ponytail: deliberately minimal. Does not wrap business parsing; only HTTP,
proxy, rate-limit, retry, and cache concerns live here.
"""

from __future__ import annotations

import asyncio
import json
import time
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Dict, Optional

import aiohttp

from dgra_config import APIConfig, DGRAGlobalConfig


class AdaptiveRateLimiter:
    """
    Self-tuning rate limiter moved from dgra_api.py.

    Starts conservative, backs off on HTTP 429, and slowly recovers after
    successful requests.
    """

    def __init__(
        self,
        initial_rate: float = 2.0,
        min_rate: float = 0.1,
        max_rate: float = 5.0,
        success_threshold: int = 5,
        rate_boost: float = 1.2,
        rate_cut: float = 0.5,
    ):
        self.rate = initial_rate
        self.min_rate = min_rate
        self.max_rate = max_rate
        self.success_threshold = success_threshold
        self.rate_boost = rate_boost
        self.rate_cut = rate_cut
        self._success_streak = 0
        self._429_count = 0
        self._last_request_time = 0.0
        self._lock = asyncio.Lock()

    async def acquire(self) -> float:
        async with self._lock:
            now = time.time()
            min_interval = 1.0 / self.rate
            elapsed = now - self._last_request_time
            if elapsed < min_interval:
                await asyncio.sleep(min_interval - elapsed)
            self._last_request_time = time.time()
            return self.rate

    def report_success(self) -> None:
        self._success_streak += 1
        if self._success_streak >= self.success_threshold:
            self.rate = min(self.max_rate, self.rate * self.rate_boost)
            self._success_streak = 0

    def report_429(self) -> None:
        self._success_streak = 0
        self._429_count += 1
        self.rate = max(self.min_rate, self.rate * self.rate_cut)

    def report_error(self, is_429: bool = False) -> None:
        if is_429:
            self.report_429()
        else:
            self._success_streak = max(0, self._success_streak - 1)

    @property
    def stats(self) -> Dict[str, Any]:
        return {
            "current_rate": round(self.rate, 2),
            "success_streak": self._success_streak,
            "429_count": self._429_count,
        }


class APIHub:
    """
    Single shared HTTP hub for all GPA external API calls.

    Usage:
        async with APIHub(config, cache) as hub:
            result = await hub.request("ensembl", "/lookup/symbol/...")
            async with hub.session.get(url) as resp:
                ...
    """

    _COMMON_PROXIES = [
        "http://127.0.0.1:7897",
        "http://127.0.0.1:7890",
        "http://127.0.0.1:7891",
        "http://127.0.0.1:1080",
        "http://127.0.0.1:10808",
        "http://127.0.0.1:10809",
        "http://127.0.0.1:8000",
        "http://127.0.0.1:8080",
        "http://127.0.0.1:3128",
    ]

    def __init__(
        self,
        config: DGRAGlobalConfig,
        cache: Optional[Any] = None,
        *,
        detect_proxy: bool = True,
        audit_trail: Optional[Any] = None,
    ):
        self.config = config
        self.cache = cache
        self._detect_proxy_on_setup = detect_proxy
        self.audit_trail = audit_trail
        self._session: Optional[aiohttp.ClientSession] = None
        self._proxy_url: Optional[str] = None
        self._rate_limiters: Dict[str, AdaptiveRateLimiter] = {}

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------
    async def setup(self) -> None:
        if self._detect_proxy_on_setup:
            self._proxy_url = await self._detect_proxy()
        self._session = aiohttp.ClientSession(
            connector=aiohttp.TCPConnector(limit=50, limit_per_host=20),
            timeout=aiohttp.ClientTimeout(total=120),
            trust_env=False,  # ponytail: never trust env; use explicit proxy
        )

    async def close(self) -> None:
        if self._session:
            await self._session.close()
            self._session = None

    async def __aenter__(self) -> "APIHub":
        await self.setup()
        return self

    async def __aexit__(self, *exc: Any) -> None:
        await self.close()

    # ------------------------------------------------------------------
    # Proxy detection
    # ------------------------------------------------------------------
    @staticmethod
    async def _probe_endpoint(
        proxy: Optional[str],
        timeout: float = 3.0,
        session: Optional[aiohttp.ClientSession] = None,
    ) -> bool:
        test_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            "?db=clinvar&term=BRCA1%5BGene%5D&retmode=json&retmax=1"
        )
        close_session = session is None
        if session is None:
            session = aiohttp.ClientSession(
                trust_env=False,
                timeout=aiohttp.ClientTimeout(total=timeout),
            )
        try:
            async with session.get(test_url, proxy=proxy) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    count = data.get("esearchresult", {}).get("count")
                    return bool(count and int(count) > 0)
        except Exception:
            return False
        finally:
            if close_session:
                await session.close()
        return False

    async def _detect_proxy(self) -> Optional[str]:
        """Pick direct connection if possible, otherwise first working local proxy."""
        explicit = self.config.proxy
        if explicit and explicit != "__DIRECT__":
            return explicit

        print("[APIHub] Testing direct connection to NCBI...")
        if await self._probe_endpoint(proxy=None, timeout=3.0):
            print("[APIHub] Direct connection OK — no proxy needed")
            return None
        print("[APIHub] Direct connection FAILED")

        semaphore = asyncio.Semaphore(3)
        results: Dict[str, bool] = {}

        async def _test_one(proxy: str) -> None:
            async with semaphore:
                ok = await self._probe_endpoint(proxy=proxy, timeout=5.0)
                results[proxy] = ok
                status = "OK" if ok else "FAILED"
                print(f"[APIHub] Proxy {proxy} {status}")

        tasks = [asyncio.create_task(_test_one(p)) for p in self._COMMON_PROXIES]
        await asyncio.gather(*tasks, return_exceptions=True)

        for proxy in self._COMMON_PROXIES:
            if results.get(proxy):
                return proxy

        print("[APIHub] No working proxy found — using direct (expect failures)")
        return None

    def proxy_for(self, api_name: str) -> Optional[str]:
        """Return proxy for a specific API, respecting per-API and route-map overrides."""
        route_map = getattr(self.config, "_proxy_route_map", None)
        if route_map is not None:
            return route_map.get_proxy(api_name)
        cfg = self.config.apis.get(api_name)
        if cfg and cfg.proxy and cfg.proxy != "__DIRECT__":
            return cfg.proxy
        return self._proxy_url

    # ------------------------------------------------------------------
    # Rate limiting
    # ------------------------------------------------------------------
    def _get_rate_limiter(self, api_name: str) -> AdaptiveRateLimiter:
        if api_name not in self._rate_limiters:
            cfg = self.config.apis.get(api_name)
            initial_rate = getattr(cfg, "rate_limit_per_sec", 2.0) if cfg else 2.0
            self._rate_limiters[api_name] = AdaptiveRateLimiter(
                initial_rate=initial_rate,
                min_rate=0.05,
                max_rate=initial_rate,
            )
        return self._rate_limiters[api_name]

    async def rate_limit(self, api_name: str) -> None:
        limiter = self._get_rate_limiter(api_name)
        await limiter.acquire()

    # ------------------------------------------------------------------
    # Session access (for modules with custom batch logic)
    # ------------------------------------------------------------------
    @property
    def session(self) -> aiohttp.ClientSession:
        if self._session is None:
            raise RuntimeError("APIHub session not initialized; use async with APIHub(...)")
        return self._session

    # ------------------------------------------------------------------
    # Unified request with cache, rate limit, retry
    # ------------------------------------------------------------------
    async def request(
        self,
        api_name: str,
        endpoint: str,
        method: str = "GET",
        params: Optional[Dict[str, Any]] = None,
        json_body: Optional[Any] = None,
        headers: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        """
        Make a single API request through the hub.

        Returns a dict with keys: data, http_status, from_cache, confidence,
        and optionally error.
        """
        cfg = self.config.apis.get(api_name)
        if cfg is None:
            raise ValueError(f"No APIConfig for {api_name}")

        base_url = cfg.base_url.rstrip("/")
        url = f"{base_url}/{endpoint.lstrip('/')}"

        cache_key_params = {"url": url, **(params or {})}
        if json_body is not None:
            cache_key_params["_body"] = json.dumps(
                json_body, sort_keys=True, separators=(",", ":")
            )

        if self.cache:
            cached = self.cache.get(api_name, **cache_key_params)
            if cached:
                if self.audit_trail:
                    self.audit_trail.record_api_call(api_name, url, status=cached.get("http_status"), from_cache=True)
                return {
                    "data": cached["data"],
                    "http_status": cached["http_status"],
                    "from_cache": True,
                    "confidence": cached["confidence"],
                }

        if self.config.offline_mode:
            if self.audit_trail:
                self.audit_trail.record_api_call(api_name, url, status=None, from_cache=False, error="offline mode")
            return {
                "data": None,
                "http_status": None,
                "from_cache": False,
                "confidence": "low",
                "error": "Offline mode: no cached data available",
            }

        last_error: Optional[Exception] = None
        limiter = self._get_rate_limiter(api_name)
        req_start = time.time()

        for attempt in range(cfg.max_retries):
            try:
                await limiter.acquire()
                proxy = self.proxy_for(api_name)

                async with self._session.request(
                    method=method,
                    url=url,
                    params=params,
                    json=json_body,
                    headers=headers,
                    proxy=proxy,
                    timeout=aiohttp.ClientTimeout(total=cfg.timeout),
                ) as response:
                    status = response.status

                    duration_ms = round((time.time() - req_start) * 1000, 2)
                    if status == 200:
                        limiter.report_success()
                        try:
                            data = await response.json()
                        except Exception:
                            text = await response.text()
                            try:
                                data = json.loads(text)
                            except Exception:
                                err = f"HTTP 200 but response is not valid JSON (len={len(text)})"
                                if self.audit_trail:
                                    self.audit_trail.record_api_call(api_name, url, status=status, duration_ms=duration_ms, error=err)
                                return {
                                    "data": None,
                                    "http_status": status,
                                    "from_cache": False,
                                    "confidence": "low",
                                    "error": err,
                                }
                        if self.cache:
                            self.cache.set(
                                api_name=api_name,
                                response_data=data,
                                http_status=status,
                                confidence="medium",
                                **cache_key_params,
                            )
                        if self.audit_trail:
                            self.audit_trail.record_api_call(api_name, url, status=status, duration_ms=duration_ms)
                        return {
                            "data": data,
                            "http_status": status,
                            "from_cache": False,
                            "confidence": "medium",
                        }

                    if status == 404:
                        if self.cache:
                            self.cache.set(
                                api_name=api_name,
                                response_data={"error": "not_found", "status": 404},
                                http_status=404,
                                confidence="medium",
                                ttl_days=7,
                                **cache_key_params,
                            )
                        if self.audit_trail:
                            self.audit_trail.record_api_call(api_name, url, status=status, duration_ms=duration_ms)
                        return {
                            "data": None,
                            "http_status": 404,
                            "from_cache": False,
                            "confidence": "medium",
                            "error": "Not found",
                        }

                    if status == 429:
                        limiter.report_429()
                        retry_after = response.headers.get("Retry-After")
                        wait = int(retry_after) if retry_after else max(
                            cfg.retry_delay * (2 ** attempt),
                            max(1.0, 1.0 / max(limiter.rate, 0.05)),
                        )
                        last_error = Exception(f"Rate limited (429), retry after {wait}s")
                        if self.audit_trail:
                            self.audit_trail.record_api_call(api_name, url, status=status, duration_ms=duration_ms, error=str(last_error))
                        print(
                            f"[APIHub] {api_name}: 429, waiting {wait}s "
                            f"before retry {attempt + 1}/{cfg.max_retries}"
                        )
                        await asyncio.sleep(wait)
                        continue

                    if status in (502, 503, 504) or status >= 500:
                        last_error = Exception(f"Server error {status}")
                        if self.audit_trail:
                            self.audit_trail.record_api_call(api_name, url, status=status, duration_ms=duration_ms, error=str(last_error))
                        wait = cfg.retry_delay * (2 ** attempt)
                        print(
                            f"[APIHub] {api_name}: HTTP {status}, retrying in {wait}s "
                            f"({attempt + 1}/{cfg.max_retries})"
                        )
                        await asyncio.sleep(wait)
                        continue

                    text = await response.text()
                    err = f"HTTP {status}: {text[:200]}"
                    if self.audit_trail:
                        self.audit_trail.record_api_call(api_name, url, status=status, duration_ms=duration_ms, error=err)
                    return {
                        "data": None,
                        "http_status": status,
                        "from_cache": False,
                        "confidence": "low",
                        "error": err,
                    }

            except asyncio.TimeoutError:
                last_error = Exception(f"Timeout after {cfg.timeout}s")
                if self.audit_trail:
                    self.audit_trail.record_api_call(api_name, url, status=None, error=str(last_error))
                wait = cfg.retry_delay * (2 ** attempt)
                print(f"[APIHub] {api_name}: Timeout, retrying in {wait}s ({attempt + 1}/{cfg.max_retries})")
                await asyncio.sleep(wait)
                continue

            except aiohttp.ClientError as e:
                last_error = Exception(f"Connection error: {e}")
                if self.audit_trail:
                    self.audit_trail.record_api_call(api_name, url, status=None, error=str(last_error))
                wait = cfg.retry_delay * (2 ** attempt)
                print(f"[APIHub] {api_name}: Connection error ({e}), retrying in {wait}s ({attempt + 1}/{cfg.max_retries})")
                await asyncio.sleep(wait)
                continue

        duration_ms = round((time.time() - req_start) * 1000, 2)
        final_err = str(last_error) if last_error else "All retries failed"
        if self.audit_trail:
            self.audit_trail.record_api_call(api_name, url, status=None, duration_ms=duration_ms, error=final_err)
        return {
            "data": None,
            "http_status": None,
            "from_cache": False,
            "confidence": "low",
            "error": final_err,
        }


@asynccontextmanager
async def get_hub(
    config: Optional[DGRAGlobalConfig] = None,
    cache: Optional[Any] = None,
):
    """Convenience context manager for the common case."""
    if config is None:
        config = DGRAGlobalConfig.from_env()
    async with APIHub(config, cache) as hub:
        yield hub
