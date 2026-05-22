# DGRA 更新日志

## [v0.5.2] - 2026-05-22

### 核心功能 — VEP Canonical Reannotation (Transcript Discrepancy Fix)

**问题**：ANNOVAR/VEP/SnpEff 选择的 "首选转录本" 与 Ensembl canonical 不一致，导致非编码转录本（NR_/XM_）被标注为 `splice_donor_variant`/`HIGH`，但 canonical 蛋白编码转录本（NM_/ENST_）下同一变异实为 `upstream_gene_variant`/`MODIFIER`。`HIGH` 被错误送入 `classify_variant_tier()`，产生假阳性 Tier 1。

**解决方案**：
- **Step 1.5 新增 `batch_query_vep_region()`**：收集 Step 1 中 `TRANSCRIPT_DISCREPANCY` 变异，用 Ensembl VEP API 以 `canonical=1&domains=1&protein=1&hgvs=1&mane_select=1` 重新注释
- **解析优先级**：canonical → MANE Select → protein_coding
- **覆盖字段**：`consequence`, `impact`, `hgvsc`, `hgvsp`, `transcript`
- **`transcript_warning.vep_reannotation`**：记录 original vs canonical 对比，Markdown 报告中加 ⚠️ 标注
- **失败/离线降级**：`quality_confidence="LOW"`, `tier_confidence="LOW"`, `vep_reannotation_failed=True`
- **Domain mapping 顺序修正**：Step 1.5 在 Step 4 之前，修正后的 HGVSp 自动流入 UniProt 功能域映射

**典型案例 — CRIP2 chr14:105473030**：
| 阶段 | 转录本 | 后果 | 影响 | 最终分级 |
|:---|:---|:---|:---|:---|
| 原始（ANNOVAR） | `NR_073082` | `splice_donor_variant` | **HIGH** | 可能 Tier 1 |
| VEP reannotation | `NM_001312` | `upstream_gene_variant` | **MODIFIER** | **Tier 3** |

**验证**：`test_vep_reannotation_e2e.py` 端到端测试通过（9 步 Pipeline 验证）。

**相关文件**：
- `scripts/dgra_api.py` — 新增 `query_ensembl_vep_region()`, `_parse_vep_batch_response()`, `batch_query_vep_region()` (+216 行)
- `scripts/dgra_core.py` — Step 1.5 插入 `run_dgra_pipeline()` (+104 行)，`transcript_warning` fallback confidence 降级，`_format_vep_reannotation_note()` 报告增强 (+~50 行)
- `scripts/test_vep_reannotation.py` — 3 个单元测试（canonical/MANE/protein_coding fallback）
- `scripts/test_vep_reannotation_e2e.py` — CRIP2 端到端测试（9 步验证）

### 核心逻辑修正（v0.5.2 同时包含）
- **Multi-hit 不再升级变异**：只标记 multi-hit 基因，各变异独立分级（Tier 1: 301→4 突变）
- **ClinVar 中文注释支持**：`_clinvar_pathogenic` 同时匹配 "Pathogenic" 和 "致病"
- **新增 Priority 1c**：ClinVar 致病 + HIGH + 组织相关无路可走 → Tier 1（CD36 正确分级）
- **Transcript discrepancy 降级**：NR_/XM_ 非编码转录本标注 HIGH，若 canonical 为 ENST 蛋白编码 → 降级为 MODERATE
- **统计格式**："X 基因 / Y 突变" 双维度
- **报告位点格式**：强制包含 `CHROM:POS:REF>ALT`

**提交信息**：`485b851` v0.5.2 - VEP canonical reannotation + transcript discrepancy fix

---

## [v0.4.4] - 2026-05-20

### 重大改进 - Multi-hit 致病性证据过滤

**问题**：原始 multi-hit 规则将同一基因内所有变异都 elevation 到 Tier 1，导致大量 false positive。556 个变异的供者分析中，334 个被标为 Tier 1，绝大多数是正常多态性。

**解决方案**：
- 新增 `_variant_has_pathogenic_evidence()` 函数，定义致病性证据的 3 条标准（或的关系）：
  1. **影响蛋白质功能域** + 目标组织表达（GTEx TPM ≥ 1.0）
  2. **ClinVar 致病性/可能致病性** 或 **IMPACT=HIGH** 或 **gnomAD AF < 0.001**
  3. **剪切位点变化**
- **ClinVar benign 排除**：即使落在功能域内，ClinVar benign 的变异也不视为致病证据
- **精准 elevation**：只 elevation 自身满足致病性条件的变异，不是整个基因的全部变异
- **HLA 排除保留**：继续排除 HLA 基因的天然多态性

**效果**：
| 版本 | Tier 1 | Multi-hit 基因 |
|------|--------|----------------|
| 原始 | 334 | 86 |
| v0.4.3 (HLA 排除) | 238 | 43 |
| **v0.4.4** | **28** | **14** |

### 报告格式改进

- Tier 1 按基因分组，展示详细变异表格：染色体位置、转录本、变异名称、功能域、合子型、ClinVar、原因
- 每个变异增加详细说明：影响程度、后果、功能域位置、组织相关性
- 方法学附录改为中文

### 提交信息
`87cf5c0` v0.4.4 - Multi-hit pathogenic evidence filtering + report format overhaul

---

## [v0.4.3] - 2026-05-20

### 改进 - HLA 基因从 multi-hit elevation 中排除

**问题**：HLA-A/B/C 在 WES 中各有 30-38 个变异，触发 multi-hit 后被全部 elevation 到 Tier 1。但 HLA 是人类基因组多态性最高的区域，这些变异是正常免疫多样性，不代表致病性。

**解决方案**：
- 定义 HLA 基因排除集合（HLA-A/B/C/DRB1/DQA1/DQB1/DPA1/DPB1/E/F/G 等 + MICA/MICB/TAP1/TAP2）
- HLA 基因的多态性变异不再触发 elevation
- HLA 仍在 multi-hit 列表中报告，但 tier 不会被强制降至 1

**效果**：Tier 1 从 334 → 238（减少 96 个 false positive）

### 提交信息
`bcc4e4c` v0.4.3 - Exclude HLA genes from multi-hit elevation

---

## [v0.4.2] - 2026-05-20

### 修复 - GTEx API v2 迁移 + 并发优化

**问题**：GTEx API v1 (`/rest/v1`) 返回空结果，API 覆盖率 0%。

**解决方案**：
- 迁移到 GTEx API v2 (`/api/v2/expression/medianGeneExpression`)
- 支持 versioned gencodeId 解析（两步查找 + 缓存）
- 并发优化：`asyncio.Semaphore(20)` + 30-gene chunks + 0.5s 批间间隔
- 连接池扩容：TCPConnector `limit=50, limit_per_host=20`
- 超时延长：Ensembl 15→30s, UniProt 20→45s
- 单位更新：GTEx 从 RPKM 改为 TPM

**效果**：GTEx 覆盖率从 0% → **304/309 (98.4%)**

### 提交信息
`9a5a0d8` v0.4.2 - GTEx API v2 fix + UniProt batch concurrency optimization

---

## [v0.4.1] - 2026-05-19

### 关键修复

| 修复 | 说明 |
|------|------|
| **UniProt fragment 选择** | `size=1` → `size=5`，选择 reviewed + 最长的 canonical entry，避免选择 60aa 的 fragment |
| **Tissue 默认值移除** | `--tissue` 参数设为 `required=True`，禁止默认 tissue |
| **Cache JSON 解码** | 增加 `.decode('utf-8')` 处理 bytes 类型的 JSON 响应 |
| **硬编码日期** | 替换为 `datetime.now().isoformat()` |
| **UniProt feature 大小写** | 统一 `.upper()` 处理 |

### 提交信息
`39dedbe` v0.4.1 - critical fixes

---

## [v0.4.0] - 2026-05-19

### 首次发布 - API-first 架构

**重大重构**：从离线静态字典迁移到 API-first 架构

- **Phase 1**：替换 `MANE_SELECT` / `PROTEIN_DOMAINS` / `tissue_gene_lists` 为实时 API 查询
- **Ensembl REST API**：基因注释、转录本校正、canonical 转录本选择
- **UniProt REST API**：蛋白质功能域映射
- **GTEx API**：组织特异性表达
- **async 并发**：批量 API 查询 + SQLite 缓存 30 天
- **离线模式**：`--offline` 跳过 API，使用缓存 + 本地回退

### 提交信息
`edc3f75` DGRA v0.4.0 - API-first with offline archive

---

## 版本概览

| 版本 | 日期 | 主题 |
|------|------|------|
| v0.4.4 | 2026-05-20 | Multi-hit 致病性证据过滤 |
| v0.4.3 | 2026-05-20 | HLA 排除 |
| v0.4.2 | 2026-05-20 | GTEx v2 + 并发优化 |
| v0.4.1 | 2026-05-19 | 关键修复 |
| v0.4.0 | 2026-05-19 | API-first 架构 |
