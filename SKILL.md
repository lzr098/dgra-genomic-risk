---
name: dgra-genomic-risk
description: >
  Donor Genomic Risk Assessment (DGRA) v0.4. Analyzes donor VCF variants with
  three-tier risk classification. API-first with Ensembl/UniProt/GTEx live queries
  and 30-day SQLite cache. Offline archive mode: online analyses auto-save API
  results to local JSON, subsequent offline runs load archived data for identical
  results. Tissue-context adaptive: hematopoietic, cardiovascular, hepatic, renal,
  neurological. Use when evaluating donor genomic variants for transplant or
  intervention contexts.
---

# DGRA: Donor Genomic Risk Assessment v0.4

## When to Use

- Analyzing donor VCF/variant call files before hematopoietic stem cell transplantation
- Evaluating donor variants for any tissue context (cardiovascular, hepatic, renal, neurological)
- Determining if donor genetic variants affect collection safety (PBSC vs BM)
- Cross-checking patient somatic driver mutations against donor germline
- Offline analysis in clinical environments without internet access

## When NOT to Use

- General genetic counseling without tissue/intervention context
- Somatic tumor-only variant interpretation
- Population genetics or ancestry analysis

## Core Principles

### 1. Tissue-Context Adaptive

The same variant gets different tiers depending on clinical context.
Set context with --tissue. Default: hematopoietic.

### 2. Three-Tier Classification

| Tier | Name | Action |
|------|------|--------|
| 1 | Action Required | Pre-transplant intervention or exclusion |
| 2 | Inform and Monitor | Informed consent, post-intervention monitoring |
| 3 | No Concern | Document and dismiss |

### 3. API-First with Offline Archive

Primary sources (live APIs, cached 30 days): Ensembl REST, UniProt REST, GTEx Portal.

Offline archive (automatic):
- Every online analysis saves per-gene API results to references/offline_data/{gene}.json
- Subsequent offline runs load archived data = identical results to online
- If no archive exists, falls back to special_gene_lists + conservative rules

This means offline mode is NOT a downgrade -- it is "last-known-good state replay".

### 4. Confidence Annotation

HIGH = all APIs responded. MEDIUM = one API failed. LOW = multiple APIs failed or offline without archive.

## Input Format

Annotated variant table (CSV/TSV) with columns:
CHROM, POS, REF, ALT, GENE, Feature, EXON, IMPACT, Consequence, HGVSp, HGVSc, CLIN_SIG, GT, DP, GQ, VAF, gnomAD_AF (optional)

## CLI Usage

Normal mode (queries live APIs, caches + archives results):
  python3 scripts/dgra_core.py --input donor_variants.tsv --tissue hematopoietic --output report.md --json results.json

Offline mode (loads archived data if available):
  python3 scripts/dgra_core.py --input donor_variants.tsv --tissue hematopoietic --offline --output offline_report.md

Available tissue profiles: hematopoietic (default), cardiovascular, hepatic, renal, neurological.
Add profiles by editing references/tissue_context.json -- no code changes needed.

## Output

1. Markdown report -- clinical decision-oriented with tier sections
2. Structured JSON -- machine-readable for downstream integration

## Offline Archive Mechanism

Location: references/offline_data/{gene}.json

Saved automatically after every online analysis:
- ensembl: canonical transcript, biotype, description, coordinates
- uniprot: domains, GO terms, sequence length
- gtex: tissue-specific expression (RPKM)
- tissue_profile: which profile was used
- saved_at: ISO timestamp

Loaded automatically in offline mode:
- If archive exists for a gene, uses it (same as online result)
- If archive missing, falls back to special_gene_lists only
- Archive survives server restarts (file-based, not in-memory)

## Limitations

1. gnomAD coverage gaps: KIR cluster, X-linked genes, highly polymorphic regions
2. Protein domains: UniProt annotation coverage varies by gene
3. Phase confirmation: Multi-hit cis/trans requires family/trio or long-read sequencing
4. GTEx proxy: Tissue-level RPKM may not reflect cell-type-specific expression
5. Population bias: gnomAD primarily European ancestry

## File Structure

dgra-genomic-risk/
  SKILL.md              -- this file
  DESIGN_v0.4.md        -- architecture documentation
  config.json           -- skill metadata
  requirements.txt      -- aiohttp dependency
  scripts/
    dgra_core.py        -- main pipeline (async, v0.4)
    dgra_config.py      -- configuration
    dgra_cache.py       -- SQLite cache with TTL
    dgra_api.py         -- async API clients
  references/
    tissue_context.json -- tissue profiles
    offline_data/       -- per-gene API archive (auto-created)
  cache/
    dgra_cache.db       -- SQLite API response cache (auto-created)

## Version

DGRA v0.4 -- 2026-05-19

Key updates from v0.3:
- API-first architecture with live Ensembl/UniProt/GTEx queries and 30-day SQLite cache
- Offline archive: online analyses auto-save per-gene JSON, offline loads identical data
- Async batch queries with rate limiting and retry
- special_gene_lists as irreplaceable clinical rules
- Confidence annotation per result (HIGH/MEDIUM/LOW)
