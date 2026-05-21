#!/usr/bin/env python3
"""
DGRA Configuration Manager
Phase 1 - v0.4 Architecture

Handles API keys, timeouts, retry policies, and offline mode.
"""

import os
import json
from dataclasses import dataclass, field
from typing import Dict, Optional, List
from pathlib import Path


@dataclass
class APIConfig:
    """API endpoint configuration with rate limiting and retry policy."""
    base_url: str
    timeout: float = 30.0  # seconds
    max_retries: int = 3
    retry_delay: float = 1.0  # initial backoff (exponential)
    rate_limit_per_sec: float = 10.0  # requests per second
    api_key: Optional[str] = None
    proxy: Optional[str] = None  # HTTP/HTTPS proxy URL
    

@dataclass
class DGRAGlobalConfig:
    """Global DGRA configuration including all API endpoints and runtime flags."""
    
    # Offline mode: skip all API calls, use cache + local overrides only
    offline_mode: bool = False
    
    # v0.5 P1-8: Gene list sync settings
    gene_sync_enabled: bool = True
    gene_sync_ttl_days: int = 7
    
    # v0.4.5: Somatic mode for tumor driver analysis (not germline donor screening)
    somatic_mode: bool = False
    
    # Cache settings
    cache_db_path: Path = field(default_factory=lambda: Path(__file__).parent.parent / "cache" / "dgra_cache.db")
    cache_ttl_days: int = 30  # Time-to-live for cached API responses
    
    # Analysis thresholds (same as v0.3)
    min_dp: int = 20
    min_gq: float = 90.0
    common_af_threshold: float = 0.01
    low_af_threshold: float = 0.001
    vaf_deviation_threshold: float = 0.20
    
    # Tissue context profile name
    tissue_profile: str = "hematopoietic"
    
    # v0.5 P1-7: Multi-organ assessment — evaluate multiple profiles simultaneously
    # Each profile runs independently; joint report takes max tier across profiles
    multi_organ_profiles: Optional[List[str]] = None
    
    # Confidence thresholds for manual review flagging
    high_confidence_min_apis: int = 3  # Need at least N APIs responding for HIGH
    
    # API endpoint configurations (no keys needed for public endpoints)
    apis: Dict[str, APIConfig] = field(default_factory=lambda: {
        "ensembl": APIConfig(
            base_url="https://rest.ensembl.org",
            timeout=20.0,
            max_retries=2,
            retry_delay=1.0,
            rate_limit_per_sec=10.0,
        ),
        "uniprot": APIConfig(
            base_url="https://rest.uniprot.org",
            timeout=25.0,
            max_retries=2,
            retry_delay=2.0,
            rate_limit_per_sec=5.0,
        ),
        "gtex": APIConfig(
            base_url="https://gtexportal.org/api/v2",
            timeout=20.0,
            max_retries=2,
            retry_delay=2.0,
            rate_limit_per_sec=3.0,
        ),
        "gnomad": APIConfig(
            base_url="https://gnomad.broadinstitute.org/api",
            timeout=15.0,
            max_retries=2,
            retry_delay=2.0,
            rate_limit_per_sec=2.0,
        ),
        "ncbi_eutils": APIConfig(
            base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
            timeout=15.0,
            max_retries=2,
            retry_delay=1.0,
            rate_limit_per_sec=3.0,  # NCBI: 3/sec without API key
        ),
        "clinvar_eutils": APIConfig(
            base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
            timeout=15.0,
            max_retries=2,
            retry_delay=1.0,
            rate_limit_per_sec=3.0,
        ),
        "hgnc": APIConfig(
            base_url="https://rest.genenames.org",
            timeout=15.0,
            max_retries=2,
            retry_delay=1.0,
            rate_limit_per_sec=5.0,
        ),
    })
    
    # Override files (local patches for API blind spots)
    override_files: Dict[str, Path] = field(default_factory=lambda: {
        "pseudogene": Path(__file__).parent.parent / "references" / "pseudogene_config.json",
        "api_corrections": Path(__file__).parent.parent / "references" / "api_corrections.json",
        "tissue_profiles": Path(__file__).parent.parent / "references" / "tissue_context.json",
    })
    
    @classmethod
    def from_env(cls) -> "DGRAGlobalConfig":
        """Load configuration from environment variables."""
        config = cls()
        
        # Offline mode
        if os.environ.get("DGRA_OFFLINE", "").lower() in ("1", "true", "yes"):
            config.offline_mode = True
        
        # Override cache path
        if cache_path := os.environ.get("DGRA_CACHE_PATH"):
            config.cache_db_path = Path(cache_path)
        
        # Override tissue profile
        if tissue := os.environ.get("DGRA_TISSUE"):
            config.tissue_profile = tissue
        
        # HTTP/HTTPS proxy support
        for api_name in config.apis:
            proxy = os.environ.get("HTTP_PROXY") or os.environ.get("http_proxy") or os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
            if proxy:
                config.apis[api_name].proxy = proxy
        
        # NCBI API key (increases rate limit from 3/sec to 10/sec)
        if ncbi_key := os.environ.get("NCBI_API_KEY"):
            config.apis["ncbi_eutils"].api_key = ncbi_key
            config.apis["clinvar_eutils"].api_key = ncbi_key
            config.apis["ncbi_eutils"].rate_limit_per_sec = 10.0
            config.apis["clinvar_eutils"].rate_limit_per_sec = 10.0
        
        # Custom API timeouts
        for api_name in config.apis:
            env_key = f"DGRA_{api_name.upper()}_TIMEOUT"
            if timeout := os.environ.get(env_key):
                config.apis[api_name].timeout = float(timeout)
        
        return config
    
    def get_override(self, name: str) -> Optional[dict]:
        """Load a local override JSON file if it exists."""
        path = self.override_files.get(name)
        if path and path.exists():
            with open(path, 'r') as f:
                return json.load(f)
        return None
