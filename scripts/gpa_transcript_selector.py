#!/usr/bin/env python3
"""
GPA Transcript Selector Module (v0.10.0)

Disease-aware transcript selection for genes with multiple isoforms.
1. Rule-based scoring (canonical, MANE, tissue expression, impact)
2. Ambiguity detection (top scores within <5 points)
3. LLM-assisted selection (only when ambiguous, reuses existing LLM pattern)

Outputs: primary transcript + alternative transcripts list.
"""

import asyncio
import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import aiohttp

try:
    from dgra_config import DGRAGlobalConfig
except Exception:
    DGRAGlobalConfig = None  # type: ignore[misc,assignment]

try:
    from api_hub import APIHub
except Exception:
    APIHub = None  # type: ignore[misc,assignment]

logger = logging.getLogger(__name__)

# Tissue profile → special gene lists mapping (loaded from tissue_context.json)
_TISSUE_CONTEXT_PATH = (
    Path(__file__).resolve().parent.parent / "references" / "tissue_context.json"
)


def _load_tissue_context() -> Dict[str, Any]:
    """Load tissue context profiles for gene relevance scoring."""
    try:
        with open(_TISSUE_CONTEXT_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


_tissue_context = _load_tissue_context()


@dataclass
class TranscriptSelectionResult:
    """Result of transcript selection for a single gene."""
    primary: Dict[str, Any]
    alternatives: List[Dict[str, Any]]
    is_ambiguous: bool
    method: str  # canonical / tissue_expression / llm_disease_match / ambiguous
    selection_reason: str


class TranscriptSelector:
    """
    Select the optimal transcript for a gene based on:
    - VEP annotations (canonical, MANE flags)
    - Tissue profile relevance
    - Disease description (LLM-assisted when ambiguous)
    """

    def __init__(
        self,
        tissue_profile: str = "general",
        disease_description: Optional[str] = None,
        llm_api_key: Optional[str] = None,
        llm_model: str = "gpt-4o-mini",
        ambiguity_threshold: int = 5,
    ):
        """
        Args:
            tissue_profile: tissue context (general, hematopoietic, cardiovascular, etc.)
            disease_description: optional clinical phenotype description
            llm_api_key: OpenAI API key for ambiguous cases
            llm_model: LLM model name
            ambiguity_threshold: score gap below which top transcripts are "ambiguous"
        """
        self.tissue_profile = tissue_profile
        self.disease_description = disease_description
        self.llm_api_key = llm_api_key or os.environ.get("OPENAI_API_KEY")
        self.llm_model = llm_model
        self.ambiguity_threshold = ambiguity_threshold

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def select(
        self,
        gene: str,
        transcripts: List[Dict[str, Any]],
    ) -> TranscriptSelectionResult:
        """
        Select primary transcript from a list of VEP transcript consequences.

        Args:
            gene: gene symbol
            transcripts: list of transcript dicts from VEP (each with canonical, mane_select, etc.)

        Returns:
            TranscriptSelectionResult with primary, alternatives, ambiguity flag
        """
        if not transcripts:
            return TranscriptSelectionResult(
                primary={},
                alternatives=[],
                is_ambiguous=False,
                method="none",
                selection_reason="No transcript consequences available",
            )

        if len(transcripts) == 1:
            tx = transcripts[0]
            return TranscriptSelectionResult(
                primary=tx,
                alternatives=[],
                is_ambiguous=False,
                method="canonical" if tx.get("canonical") else "single",
                selection_reason="Only one transcript available",
            )

        # Step 1: Rule-based scoring
        scored = []
        for tx in transcripts:
            score, reasons = self._score_transcript(tx, gene)
            scored.append({"tx": tx, "score": score, "reasons": reasons})

        # Sort by score descending
        scored.sort(key=lambda x: x["score"], reverse=True)

        # Step 2: Ambiguity detection
        top_score = scored[0]["score"]
        second_score = scored[1]["score"] if len(scored) > 1 else 0
        is_ambiguous = (top_score - second_score) < self.ambiguity_threshold

        # Step 3: LLM assist if ambiguous and disease description provided
        if is_ambiguous and self.disease_description and self.llm_api_key:
            # Gather top candidates (within ambiguity threshold)
            top_candidates = [
                s["tx"] for s in scored
                if top_score - s["score"] < self.ambiguity_threshold
            ][:3]  # Max 3 candidates for LLM
            
            # v0.9.3: Avoid asyncio.run() in async contexts
            try:
                loop = asyncio.get_running_loop()
                # We're in an async context — cannot use asyncio.run
                # Return without LLM selection; caller should use aselect()
                logging.warning(
                    "TranscriptSelector.select() called from async context. "
                    "LLM-assisted selection skipped. "
                    "Use 'await selector.aselect(gene, transcripts)' for full LLM support in async contexts."
                )
                llm_choice = None
            except RuntimeError:
                # No running loop — safe to use asyncio.run
                llm_choice = asyncio.run(
                    self._llm_assist_select(gene, top_candidates)
                )
                if llm_choice:
                    # Reorder: llm_choice becomes primary
                    primary = llm_choice
                    alternatives = [s["tx"] for s in scored if s["tx"] != primary]
                    return TranscriptSelectionResult(
                        primary=primary,
                        alternatives=alternatives,
                        is_ambiguous=True,
                        method="llm_disease_match",
                        selection_reason=f"LLM selected based on disease description '{self.disease_description[:50]}...' "
                                         f"(ambiguous: top scores {top_score} vs {second_score})",
                    )

        # Step 4: Use rule-based top choice
        primary = scored[0]["tx"]
        alternatives = [s["tx"] for s in scored[1:]]
        method = (
            "ambiguous" if is_ambiguous else
            ("tissue_expression" if any("tissue" in r for r in scored[0]["reasons"]) else "canonical")
        )
        return TranscriptSelectionResult(
            primary=primary,
            alternatives=alternatives,
            is_ambiguous=is_ambiguous,
            method=method,
            selection_reason="; ".join(scored[0]["reasons"]),
        )

    async def aselect(
        self,
        gene: str,
        transcripts: List[Dict[str, Any]],
    ) -> TranscriptSelectionResult:
        """Async version of select() with full LLM-assisted support.

        Use this when calling from an async context (e.g., within the dgra pipeline).
        Unlike select(), this method can await LLM-assisted selection without
        triggering RuntimeError from nested asyncio.run().
        """
        if not transcripts:
            return TranscriptSelectionResult(
                primary={},
                alternatives=[],
                is_ambiguous=False,
                method="none",
                selection_reason="No transcript consequences available",
            )

        if len(transcripts) == 1:
            tx = transcripts[0]
            return TranscriptSelectionResult(
                primary=tx,
                alternatives=[],
                is_ambiguous=False,
                method="canonical" if tx.get("canonical") else "single",
                selection_reason="Only one transcript available",
            )

        # Step 1: Rule-based scoring
        scored = []
        for tx in transcripts:
            score, reasons = self._score_transcript(tx, gene)
            scored.append({"tx": tx, "score": score, "reasons": reasons})

        scored.sort(key=lambda x: x["score"], reverse=True)

        # Step 2: Ambiguity detection
        top_score = scored[0]["score"]
        second_score = scored[1]["score"] if len(scored) > 1 else 0
        is_ambiguous = (top_score - second_score) < self.ambiguity_threshold

        # Step 3: LLM assist if ambiguous and disease description provided
        if is_ambiguous and self.disease_description and self.llm_api_key:
            top_candidates = [
                s["tx"] for s in scored
                if top_score - s["score"] < self.ambiguity_threshold
            ][:3]
            llm_choice = await self._llm_assist_select(gene, top_candidates)
            if llm_choice:
                primary = llm_choice
                alternatives = [s["tx"] for s in scored if s["tx"] != primary]
                return TranscriptSelectionResult(
                    primary=primary,
                    alternatives=alternatives,
                    is_ambiguous=True,
                    method="llm_disease_match",
                    selection_reason=f"LLM selected based on disease description '{self.disease_description[:50]}...' "
                                     f"(ambiguous: top scores {top_score} vs {second_score})",
                )

        # Step 4: Use rule-based top choice
        primary = scored[0]["tx"]
        alternatives = [s["tx"] for s in scored[1:]]
        method = (
            "ambiguous" if is_ambiguous else
            ("tissue_expression" if any("tissue" in r for r in scored[0]["reasons"]) else "canonical")
        )
        return TranscriptSelectionResult(
            primary=primary,
            alternatives=alternatives,
            is_ambiguous=is_ambiguous,
            method=method,
            selection_reason="; ".join(scored[0]["reasons"]),
        )

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------

    # v0.11.5: Biotype scoring penalties for non-protein-coding transcripts
    _BIOTYPE_PENALTIES = {
        # Strong penalty: these should almost never be selected for clinical interpretation
        "nonsense_mediated_decay": -20,
        "retained_intron": -20,
        "processed_pseudogene": -25,
        "transcribed_pseudogene": -25,
        "polymorphic_pseudogene": -25,
        "pseudogene": -25,
        # Moderate penalty: non-coding RNA types
        "lncRNA": -15,
        "misc_RNA": -15,
        "snRNA": -15,
        "snoRNA": -15,
        "rRNA": -15,
        "miRNA": -15,
        "scaRNA": -15,
        # Reward for reliable protein-coding
        "protein_coding": +5,
    }

    def _score_transcript(
        self,
        tx: Dict[str, Any],
        gene: str,
    ) -> Tuple[int, List[str]]:
        """
        Score a transcript based on rule-based criteria.
        v0.11.5: Added biotype filtering to prevent NMD/pseudogene selection.
        Returns (score, list_of_reasons).
        """
        score = 0
        reasons = []

        # 0. Biotype check (v0.11.5) — applied FIRST, before other scoring
        biotype = tx.get("biotype", "").lower()
        biotype_penalty = self._BIOTYPE_PENALTIES.get(biotype, 0)
        if biotype_penalty != 0:
            score += biotype_penalty
            if biotype_penalty > 0:
                reasons.append(f"protein_coding biotype (+{biotype_penalty})")
            else:
                reasons.append(f"{biotype} biotype ({biotype_penalty})")

        # 1. Canonical flag (+15) — raised from +10 to ensure canonical dominates
        #    over HIGH impact from NMD transcripts
        if tx.get("canonical"):
            score += 15
            reasons.append("canonical")

        # 2. MANE Select (+10)
        if tx.get("mane_select"):
            score += 10
            reasons.append("MANE Select")

        # 3. MANE Plus Clinical (+5)
        if tx.get("mane_plus_clinical"):
            score += 5
            reasons.append("MANE Plus Clinical")

        # 4. Tissue expression relevance (+5~15)
        tissue_bonus = self._tissue_expression_bonus(gene)
        if tissue_bonus > 0:
            score += tissue_bonus
            reasons.append(f"tissue relevance (+{tissue_bonus})")

        # 5. Impact severity (+2~5) — lowered from +3~10
        #    HIGH impact on NMD transcript should NOT outrank canonical protein_coding
        impact = tx.get("impact", "").upper()
        if impact == "HIGH":
            score += 5
            reasons.append("HIGH impact")
        elif impact == "MODERATE":
            score += 3
            reasons.append("MODERATE impact")
        elif impact == "LOW":
            score += 1
            reasons.append("LOW impact")

        # 6. Protein domain involvement (+2~6) — lowered from +3~8
        domains = tx.get("protein_domains", [])
        if domains:
            score += min(len(domains) * 2, 6)
            reasons.append(f"{len(domains)} protein domains")

        return score, reasons

    def _tissue_expression_bonus(self, gene: str) -> int:
        """
        Check if gene is in tissue profile's special gene lists.
        Returns bonus score (0~15).
        """
        profiles = _tissue_context.get("profiles", {})
        profile = profiles.get(self.tissue_profile)
        if not profile:
            return 0

        special_lists = profile.get("special_gene_lists", {})
        # Check if gene appears in any special list for this tissue
        for list_name, genes in special_lists.items():
            if gene in genes:
                # Higher bonus for core functional lists
                if list_name in ("coagulation", "fa_dna_repair", "cardiac_safety"):
                    return 15
                return 10
        return 0

    # ------------------------------------------------------------------
    # LLM-assisted selection (only for ambiguous cases)
    # ------------------------------------------------------------------

    async def _llm_assist_select(
        self,
        gene: str,
        candidates: List[Dict[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        """
        Ask LLM to pick the best transcript given disease description.
        Only called when top candidates are ambiguous (score gap < threshold).

        Returns the chosen transcript dict, or None if LLM fails.
        """
        if not self.llm_api_key:
            return None

        # Build candidate descriptions
        candidate_desc = []
        for i, tx in enumerate(candidates, 1):
            desc = (
                f"{i}. {tx.get('transcript_id', 'N/A')} — "
                f"consequence: {', '.join(tx.get('consequence_terms', []))}, "
                f"impact: {tx.get('impact', 'N/A')}, "
                f"canonical: {bool(tx.get('canonical'))}, "
                f"MANE: {bool(tx.get('mane_select'))}, "
                f"protein domains: {len(tx.get('protein_domains', []))}"
            )
            candidate_desc.append(desc)

        prompt = (
            f"You are a clinical geneticist. A patient has the following phenotype: "
            f"'{self.disease_description}'.\n\n"
            f"Gene: {gene}\n"
            f"Candidate transcripts (isoforms):\n"
            f"{'\n'.join(candidate_desc)}\n\n"
            f"Based on the disease description, which transcript is most likely to be "
            f"clinically relevant? Consider tissue expression, protein domains, and "
            f"the specific consequence on each isoform.\n\n"
            f"Reply with ONLY the number (1, 2, or 3) of the best transcript. "
            f"If uncertain, reply '1' (default to first/candidate)."
        )

        try:
            cfg = DGRAGlobalConfig.from_env() if DGRAGlobalConfig is not None else None
            async with APIHub(cfg, None, detect_proxy=False) as hub:
                async with hub.session.post(
                    "https://api.openai.com/v1/chat/completions",
                    headers={
                        "Authorization": f"Bearer {self.llm_api_key}",
                        "Content-Type": "application/json",
                    },
                    json={
                        "model": self.llm_model,
                        "messages": [{"role": "user", "content": prompt}],
                        "temperature": 0.1,
                    },
                ) as resp:
                    if resp.status != 200:
                        logger.warning(f"LLM transcript selection API error: {resp.status}")
                        return None
                    data = await resp.json()
                    content = data["choices"][0]["message"]["content"].strip()
                    # Extract number
                    match = __import__("re").search(r"\d+", content)
                    if match:
                        idx = int(match.group()) - 1
                        if 0 <= idx < len(candidates):
                            return candidates[idx]
                    # Default to first candidate
                    return candidates[0]
        except Exception as e:
            logger.warning(f"LLM transcript selection failed: {e}")
            return None

    # ------------------------------------------------------------------
    # Batch selection helper
    # ------------------------------------------------------------------

    def select_batch(
        self,
        gene_transcripts: Dict[str, List[Dict[str, Any]]],
    ) -> Dict[str, TranscriptSelectionResult]:
        """
        Select transcripts for multiple genes at once.

        Args:
            gene_transcripts: {gene_symbol: [transcript_dict, ...]}

        Returns:
            {gene_symbol: TranscriptSelectionResult}
        """
        results = {}
        for gene, transcripts in gene_transcripts.items():
            results[gene] = self.select(gene, transcripts)
        return results
