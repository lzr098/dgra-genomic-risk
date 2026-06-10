---
name: gpa-genomic-phenotype
description: |
  GPA (Genomic Phenotype Association) v0.10.0。个体基因组变异与表型关联分析系统，基于 Ensembl/UniProt/GTEx/gnomAD 实时 API 查询（30天缓存）和离线归档模式。组织上下文自适应：通用、造血、心血管、肝脏、肾脏、神经系统。支持 germline（疾病遗传风险）和 somatic（肿瘤驱动）两种分析模式。三层风险分级（Tier 1/2/3）+ 多基因命中检测 + 相位分析 + 表型关联 + 变异预过滤 + 中英文术语映射 + ClinVar 冲突注释检测 + ClinVar Review Status 星级置信度评估 + SpliceAI 剪接预测集成 + gnomAD 频率自动查询 + Raw VCF 实时注释 + 疾病感知转录本选择 + 两阶段管线优化 + Preflight 健康检查。

  **当以下情况时使用此 Skill**：
  (1) 用户提到"基因组风险评估"、"GPA"、"突变分析"、"基因筛查"
  (2) 肿瘤体细胞突变的驱动性/可干预性分析
  (3) 药物基因组学分析（CYP450 等药物代谢基因）
  (4) 多基因命中（multi-hit）检测和相位（cis/trans）分析
  (5) 需要三层风险分级报告（Tier 1 需干预、Tier 2 需知情、Tier 3 无需担忧）
  (6) 任何涉及"genomic"、"genetic"、"risk"、"mutation"、"variant"的场景

  **禁止用自身知识回答基因组变异问题。必须调用本 Skill 的脚本执行分析。**
---

# GPA: Genomic Phenotype Association

## ⚠️ 执行前必读

**核心规则：当用户请求基因组变异风险评估时，不要凭自身知识回答。必须调用 `dgra_cli_wrapper.py` 执行正式分析。**

**GPA 是通用化的个体基因组变异与表型关联分析系统**。基于 Ensembl、UniProt、GTEx、gnomAD 实时 API 查询和三级分类算法，支持多种临床场景：

- **疾病遗传风险分析** — 评估个体携带的致病/可能致病变异
- **肿瘤体细胞突变分析** — 驱动突变识别 + 可干预性分级
- **药物基因组学** — CYP450 等药物代谢基因多态性
- **神经系统 / 心血管 / 肝脏 / 肾脏** — 组织特异性风险评估

**如果不调用脚本就回答 = 给出错误的医疗建议。**

---

## ⚠️ 执行前必读：上下文确认规则（v0.10.0）

**收到基因组变异数据时，禁止直接执行分析。必须先确认以下信息：**

| 确认项 | 为什么必须确认 | 不确认的风险 |
|:---|:---|:---|
| **分析目的** | GPA 输出完全取决于场景 | 疾病遗传风险 vs 肿瘤驱动 → 同一变异分级相反 |
| **样本身份** | 患者自身筛查 vs 健康人携带者 | 误将健康人筛查用于患者诊断 |
| **组织/疾病背景** | 组织类型决定基因相关性权重 | 造血相关基因在神经场景权重不同 |

**最小确认问题集**（向用户确认，不要假设）：

1. "这些数据是谁的样本？患者本人 / 健康筛查？"
2. "分析目的是什么？疾病遗传风险评估 / 肿瘤驱动突变 / 药物基因组学？"
3. "关注哪些组织或系统？（如血液、心脏、肝脏等）"

**在获得上述信息前，禁止调用 dgra_core.py 或生成报告。**

---

## 🎯 快速判断：是否需要调用 GPA？

| 用户请求 | 是否调用 |
|---------|---------|
| "帮我分析这些基因变异" | ✅ 调用 |
| "这个突变有什么风险" | ✅ 调用 |
| "GPA 分析一下" | ✅ 调用 |
| "肿瘤突变驱动性分析" | ✅ 调用 |
| "一般性的基因突变知识" | ❌ 不需要，直接回答 |
| "TP53 突变是什么意思"（无具体样本） | ❌ 不需要 |

---

## 🚀 调用方式

### 自动分批（v0.7.1 新增，默认启用）

当变异数 > 500 时，wrapper **自动分批**处理，每批500个变异，每批有独立的5分钟超时。

**不需要任何额外参数**——只要数据量大，自动生效：

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --input-file variants.tsv \
  --tissue neurological
```

### 控制分批行为

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--batch-size` | 每批变异数 | 500 |
| `--timeout` | 每批超时（秒） | 300 |
| `--no-auto-batch` | 禁用自动分批 | 否 |

**示例：肌病基因子集（2638个变异），分6批处理**
```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --input-file myopathy_subset.tsv \
  --tissue neurological \
  --batch-size 500 \
  --output-json result.json
```

### 方式一：直接调用 wrapper（推荐）

用 `exec` 运行 `dgra_cli_wrapper.py`，传入 variant JSON 数组：

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --variants '[{"CHROM":"1","POS":12345,"REF":"A","ALT":"G","GENE":"VWF","IMPACT":"HIGH","Consequence":"missense_variant","HGVSp":"p.Arg1234Cys","HGVSc":"c.3700C>T","CLIN_SIG":"Pathogenic","GT":"0/1","DP":30,"GQ":99,"VAF":0.5}]' \
  --tissue general
```

**Somatic / 肿瘤模式**：添加 `--somatic` 标志，GPA 会按肿瘤驱动逻辑分级：

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --variants '[...]' \
  --tissue general \
  --somatic
```

Somatic 模式下：
- TSG 截断突变 + 造血相关 = **Tier 1**（核心驱动）
- 癌基因热点突变（如 IDH1 R132） = **Tier 1**
- OncoKB Oncogenic / Likely Oncogenic = **Tier 1**
- VAF > 0.5 的变异会被标记为可能的 germline 混入

**预过滤（v0.7.1 新增）**：变异数过大时，用 `--filter-preset` 先过滤再分级，减少噪音：

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --variants '[...]' \
  --tissue general \
  --filter-preset clinical
```

| Preset | 保留规则 | 用途 |
|--------|---------|------|
| `strict` | 仅 HIGH / MODERATE | 高置信度，少噪音 |
| `clinical` | HIGH / MODERATE + 剪接区 LOW + 组织相关基因同义 + ClinVar 冲突 | 推荐默认 |
| `broad` | HIGH / MODERATE / LOW | 保守，保留更多 |

**自动分批（v0.7.1）**：当变异数 > 500 时，wrapper 自动分批处理，每批500个：

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --input-file variants.tsv \
  --tissue neurological \
  --batch-size 500
```

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_core.py \
  --input /path/to/variants.tsv \
  --tissue general \
  --output /tmp/gpa_report.md \
  --json /tmp/gpa_results.json \
  --somatic \
  --spliceai                              # v0.8.0: 启用 SpliceAI 剪接预测
```

### SpliceAI 剪接预测（v0.8.0，默认关闭）

SpliceAI 仅对 canonical splice（acceptor/donor）和 splice_region 变异查询 Broad Institute lookup API，作为 VEP HIGH 剪接过调用的独立验证证据。

**必须显式开启：**

```bash
python3 ~/.openclaw/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --input-file variants.tsv \
  --tissue hematopoietic \
  --spliceai                              # 开启 SpliceAI
  --spliceai-concurrency 5                # 可选：调整并发（默认 5）
```

**效果：**
- delta=0（无剪接变化）→ 降级：HIGH 剪接过调用 → Tier 下调
- delta≥0.5（强剪接变化）→ 升级：MODERATE 剪接区 → Tier 上调
- API 失败 / 不在数据库 → graceful fallback，不阻断分析

---

### Raw VCF 注释 + 疾病感知转录本选择（v0.9.0）

GPA 支持原始未注释 VCF 输入。当检测到 raw VCF（无 `CSQ`/`ANN` INFO 字段）时，自动调用 Ensembl VEP REST API 实时注释，并通过 disease-aware 转录本选择器从多个候选转录本中挑选最相关的一条。

**必须显式提供疾病描述才能触发 disease-aware 选择**：

```bash
python3 ~/.openclaw/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --input-file raw_variants.vcf \
  --tissue hematopoietic \
  --disease-description "acute myeloid leukemia"  # 触发 disease-aware 转录本选择
```

**输入类型自动检测**：

| 输入特征 | 检测类型 | 处理方式 |
|---------|---------|---------|
| VCF 无 CSQ/ANN | `RAW_VCF` | VEP REST API 注释 → 转录本选择 → 分析 |
| VCF 有 CSQ | `ANNOTATED_VCF` | 直接解析，走原有 pipeline |
| TSV/CSV/Excel 有注释列 | `ANNOTATED_TABLE` | 直接解析，走原有 pipeline |

**转录本选择四层评分**：

1. **tissue_expression_bonus**：目标组织高表达基因额外加分
2. **consequence_bonus**：HIGH > MODERATE > LOW
3. **canonical_bonus**：canonical / MANE Select / MANE Plus Clinical 优先
4. **location_bonus**：外显子 > UTR > 内含子

**歧义处理**：
- 顶部分数差距 ≤ 10 分 → `is_ambiguous=True`
- 提供 `disease_description` + `llm_api_key` → LLM 辅助选择最相关转录本
- 未提供 → fallback 到 rule-based 最高分

**报告输出**：
- 当存在歧义或 LLM 介入时，Markdown 报告新增 **转录本选择评估章节**
- 列出 primary transcript、selection method、alternatives
- 歧义案例标 ⚠️

---

### 两阶段管线优化（v0.10.1+）

针对大型 VCF（>5,000 变异）的 API 调用优化。Phase 1 用本地规则快速 triage，Phase 2 仅对 Tier 1/2 候选变异调用外部 API。

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --input-file large_wes.vcf.gz \
  --tissue general \
  --two-phase
```

**Phase 1**（本地，<30 秒）：VEP 注释 + 本地基因列表 → 过滤 >95% 常见 SNP/低影响变异
**Phase 2**（API 仅查候选）：gnomAD、SpliceAI、表型 LLM 仅对候选变异执行

典型 germline VCF 的 API 调用量减少 **50-200x**。

---

### Preflight 健康检查（v0.10.1+）

每次接到全新分析任务时，可先执行可用性检查，确认依赖就绪后再启动耗时分析：

```bash
python3 ~/.workbuddy/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue hematopoietic \
  --preflight
```

检查范围：Python 依赖包、在线 API 连通性（8 个 API）、本地文件/目录、磁盘空间、网络/代理环境。未就绪项自动标记，可选跳过或切换离线模式。

---

## 📋 输入数据构造指南

### Variant JSON 格式

每个 variant 是一个 dict，**必填字段**：

| 字段 | 含义 | 示例 |
|-----|------|------|
| CHROM | 染色体 | "1", "X" |
| POS | 位置 | 12345 |
| REF | 参考碱基 | "A" |
| ALT | 突变碱基 | "G" |
| GENE | 基因符号 | "VWF" |

**建议填写的字段**（影响分级精度）：

| 字段 | 含义 | 示例 |
|-----|------|------|
| HGVSp | 蛋白变化 | "p.Arg1234Cys" |
| HGVSc | cDNA 变化 | "c.3700C>T" |
| IMPACT | 影响等级 | "HIGH", "MODERATE", "LOW" |
| Consequence | 突变类型 | "missense_variant", "frameshift_variant" |
| CLIN_SIG | ClinVar 状态 | "Pathogenic", "Likely_pathogenic", "Benign" |
| GT | 基因型 | "0/1" (杂合), "1/1" (纯合) |
| DP | 测序深度 | 30 |
| GQ | 基因质量 | 99 |
| VAF | 变异丰度 | 0.5 |
| gnomAD_AF | gnomAD 频率 | 0.0001 |
| **classification** | **OncoKB 致癌性** | **"Oncogenic", "Likely Oncogenic", "VUS"** |
| **is_tsg** | **是否为抑癌基因** | **"Yes" / "No"** |
| **is_oncogene** | **是否为癌基因** | **"Yes" / "No"** |

**如果用户提供了 OncoKB 标注的 CSV/表格数据，务必提取 `classification`、`is_tsg`、`is_oncogene` 字段传入 GPA，这对 somatic 模式分级至关重要。**

---

## 📊 输出解析

Wrapper 返回 JSON 结构：

```json
{
  "success": true,
  "results": {
    "meta": {...},
    "summary": {
      "tier1_count": 0,
      "tier2_count": 1,
      "tier3_count": 2,
      "multi_hit_genes": ["VWF"]
    },
    "tier1_variants": [],
    "tier2_variants": [...],
    "tier3_variants": [...],
    "multi_hit_details": [...],
    "report_markdown": "# GPA 报告..."
  },
  "report_md": "# GPA Report...",
  "stdout": "GPA Report Generated..."
}
```

### 给用户呈现的关键信息

**默认行为：分析完成后，直接读取报告文件内容并完整展示给用户。** 不要只给统计数字，不要让用户自己去看文件。

呈现内容优先级：
1. **完整 Markdown 报告** — 直接贴出报告全文（Tier 1/2/3 详情、多基因命中）
2. **若报告过长** — 先展示 Tier 1 和关键发现，再询问是否需要完整报告
3. **保存到指定路径** — 如果用户需要文件，用 `--output /path/to/file.md` 保存，不要用 /tmp

### 输出解析

Wrapper 返回 JSON 结构：

```json
{
  "success": true,
  "results": {
    "meta": {...},
    "summary": {
      "tier1_gene_count": 0,
      "tier1_variant_count": 0,
      "tier2_gene_count": 1,
      "tier2_variant_count": 2,
      "tier3_gene_count": 2,
      "tier3_variant_count": 5,
      "multi_hit_genes": ["VWF"]
    },
    "tier1_variants": [],
    "tier2_variants": [...],
    "tier3_variants": [...],
    "multi_hit_details": [...],
    "report_markdown": "# GPA 报告..."
  },
  "report_md": "# GPA Report...",
  "stdout": "GPA Report Generated..."
}
```

**呈现格式要求：**

每个变异必须展示以下字段（方便用户直接查询）：

| 字段 | 示例 | 用途 |
|------|------|------|
| **基因** | VWF | 基因符号 |
| **位点** | `chr12:6126538:G>A` | **CHROM:POS:REF:ALT，可直接在 IGV/UCSC/ClinVar 查询** |
| **转录本变化** | c.3931C>T | cDNA 水平 |
| **蛋白变化** | p.Gln1311Ter | 氨基酸水平 |
| **影响** | HIGH | VEP IMPACT |
| **类型** | stop_gained | Consequence（中英文自动映射，见下方） |
| **合子性** | 0/1 | GT |
| **ClinVar** | 致病 | 中文/英文 |
| **Tier** | 1 | GPA 分级 |
| **QC 标记** | `CLINVAR_CONFLICTING` | v0.7.1 新增：冲突注释标记 |

**Consequence 中英文映射（v0.7.1）**：
GPA 内部通过 `gpa_i18n.py` 自动标准化中英文 consequence 术语。输入可以是中文 VEP 注释（如 `错义变异`、`剪接供体变异`）或英文（`missense_variant`、`splice_donor_variant`），系统会自动映射到统一的标准术语并推断 IMPACT。

| 中文示例 | 映射后英文 | 推断 IMPACT |
|---------|-----------|------------|
| 错义变异 | missense_variant | MODERATE |
| 无义变异 | stop_gained | HIGH |
| 剪接供体变异 | splice_donor_variant | HIGH |
| 同义变异 | synonymous_variant | LOW |
| 移码变异 | frameshift_variant | HIGH |
| 框内缺失 | inframe_deletion | MODERATE |

**ClinVar 冲突注释（v0.7.1）**：
当 ClinVar 同时包含正反两种评级（如 `"良性, 致病"`、`"VUS, Pathogenic"`）时：
- 不触发 Tier 升级（weight=0）
- 标记 `CLINVAR_CONFLICTING` qc_flag
- 仍保留进入下游分析与报告
- 标准复合评级如 `"Pathogenic/Likely_pathogenic"` **不算冲突**

**ClinVar Review Status 星级（v0.7.2）**：
GPA 读取 `CLNREVSTAT` 字段，将 ClinVar 提交者星级纳入证据权重计算：

| CLNREVSTAT 文本 | 星级 | 置信度权重 |
|---|---|---|
| `practice_guideline` | ★★★★ | 0.95 |
| `reviewed_by_expert_panel` | ★★★☆ | 0.80 |
| `criteria_provided,_multiple_submitters,_no_conflicts` | ★★☆☆ | 0.55 |
| `single_submitter` | ★☆☆☆ | 0.40 |
| 缺失 / `no_assertion` / `conflicting` | — | 0.30 |

**效果**：
- Pathogenic 证据 weight = 基础值 × 星级权重（1.0 × 0.30~0.95）
- Benign 证据 weight = -0.5 × 星级权重（排除信号随星级增强）
- 单一提交者的 Pathogenic 不再自动 Tier 1（weight=0.40，需更多证据支持）
- 实践指南认可的 Pathogenic 仍高置信度（weight=0.95）

**冲突注释优先**：如果 ClinVar 评级冲突（v0.7.1 逻辑），星级评估被跳过，weight=0。
| **原因** | ClinVar pathogenic... | 为什么是这个 tier |

**若 chrom 为空或转录本选择有警告，必须标注并说明影响。**

| 场景 | 组织类型 |
|------|---------|
| 通用疾病遗传风险评估（默认） | `general` |
| 血液/肿瘤血液 | `hematopoietic` |
| 心血管 | `cardiovascular` |
| 肝脏 | `hepatic` |
| 肾脏 | `renal` |
| 神经系统 | `neurological` |

---

## 🔄 完整分析流程（含 Rescue）

### 标准流程图

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        GPA 完整分析流程 v0.10+                          │
└─────────────────────────────────────────────────────────────────────────┘

STEP 1: dgra-prefilter 预过滤
    └── 输入: raw VCF
    └── 输出: 过滤后 VCF (保留编码区/ncRNA/调控元件/ClinVar安全网)

STEP 2: VEP 注释
    └── 输入: 预过滤 VCF
    └── 输出: VEP注释 VCF (含CSQ字段)

STEP 3: 候选提取 + QC
    └── 输入: VEP注释 VCF
    └── 过滤: QUAL≥30, DP≥10, GQ≥30, GT≠0/0, MODERATE/HIGH, AF<0.01
    └── 输出: 候选变异 TSV

STEP 4: GPA 分级分析
    └── 输入: 候选变异 TSV
    └── 模式: germline (默认) 或 somatic (--somatic)
    └── 输出: Tier 1/2/3 分级结果 JSON

STEP 5: 报告 GPA 结果
    └── 呈现: Tier统计 + 关键变异详情 + 多基因命中

STEP 6: [条件触发] Phenotype Rescue Search
    └── 触发条件: Tier 1 为空, 或 Tier 1/2 不匹配表型
    └── 方法: phenotype-vcf-rescue skill (LLM Prompt + 数据库验证)
    └── 输出: Rescue候选变异 TSV (独立于GPA分级)

STEP 7: 综合报告
    └── 并列呈现: GPA结果 + Rescue结果
    └── 标注差异: "GPA Tier 3 / Rescue Priority 1"
    └── 给出建议: 家系验证 / 功能实验 / 复核影像
```

### 各步骤详细说明

| 步骤 | 工具/脚本 | 输入 | 输出 | 耗时 |
|:---|:---|:---|:---|:---|
| 1. 预过滤 | `dgra-prefilter` CLI | `patient.genotyper.vcf.gz` | `prefiltered.vcf.gz` | ~2min |
| 2. VEP注释 | VEP Docker (offline) | `prefiltered.vcf.gz` | `patient.vep.vcf.gz` | ~5-10min |
| 3. 候选提取 | `bcftools` + Python | `patient.vep.vcf.gz` | `candidates.tsv` | ~30s |
| 4. GPA分级 | `dgra_cli_wrapper.py` | `candidates.tsv` | `gpa_results.json` | ~5-30min |
| 5. 报告 | 自动解析 JSON | `gpa_results.json` | Markdown 报告 | 即时 |
| **6. Rescue** | `gpa_phenotype_rescue.py` | `patient.vep.vcf.gz` + 基因集 | `rescue.tsv` | ~1min |
| 7. 综合 | 大模型复核 + 文献 | GPA结果 + Rescue结果 | 最终诊断报告 | ~10-20min |

### Rescue 触发条件

**自动触发以下任一条件：**

1. GPA Tier 1 数量为 **0**
2. GPA Tier 1/2 变异**与患者表型不匹配**（如 PKD1L1 致多囊肾但患者无肾病）
3. 用户明确要求"深入分析"或"人工复核"

**不触发条件：**
- Tier 1 存在且表型匹配 → 直接报告，无需 Rescue

---

## 🔍 Step 6: Phenotype Rescue Search (详细)

当 Rescue 触发时，调用 `phenotype-vcf-rescue` skill：

### 6a. LLM Prompt 生成基因集

```
Patient: [年龄] [性别]
Phenotypes: [列出所有临床表型]
Distinctive features: [任何特殊组合，如"多指+大枕大池"]

Task: Identify molecular pathways and list associated genes.
Return 30-80 gene symbols, one per line, with brief rationale.
```

### 6b. 数据库验证

```bash
python3 scripts/gpa_gene_set_builder.py \
  --phenotypes "[keyword1],[keyword2]" \
  --omim-db ~/.workbuddy/data/omim/omim.db \
  --output genes_db.txt \
  --max-genes 80
```

### 6c. 合并基因集

```bash
# 合并 LLM 基因集 + 数据库基因集
cat genes_llm.txt genes_db.txt | sort -u > genes_final.txt
```

### 6d. VCF Rescue 搜索

```bash
python3 scripts/gpa_phenotype_rescue.py \
  --vcf patient.vep.vcf.gz \
  --gene-list genes_final.txt \
  --output patient.rescue.tsv \
  --patient-sex [male|female] \
  --min-impact MODERATE \
  --max-af 0.01
```

### 6e. 大模型复核 + 文献验证

对 Priority 1/2 的 Rescue 候选：
- 检查 GPA 为何将其降级（ClinVar标签？IMPACT？GTEx缺失？）
- 搜索文献确认基因-疾病关联
- 评估表型匹配度和遗传模式一致性
- 给出"suspected pathogenic"或"需排除"的结论

### 6f. 综合报告格式

```
=== Automated Tiering (GPA) ===
Tier 1: [n] | Tier 2: [n] | Tier 3: [n]

=== Phenotype Rescue Search ===
Priority 1 (high suspicion):
  - Gene | Variant | GT | AF | Rescue Reason | Phenotype Match

Priority 2 (moderate suspicion):
  - Gene | Variant | GT | AF | Rescue Reason | Phenotype Match

=== Key Discrepancies (GPA vs Rescue) ===
  - Gene X: GPA Tier 3 → Rescue Priority 1 (reason: ...)

=== Final Recommendations ===
1. [Top candidate] — [action]
2. [Second candidate] — [action]
3. [If inconclusive] — [next steps: CNV, non-coding, trio analysis]
```

---

## 🔧 离线模式

当网络不可用或 API 超时频繁时，添加 `--offline` 参数：

```bash
python3 .../dgra_cli_wrapper.py --variants '[...]' --tissue general --offline
```

离线模式使用本地缓存（`references/offline_data/` 下的基因 JSON），对于已有归档的基因结果与在线模式一致。未归档的基因 fallback 到保守规则。

---

## ❌ 常见错误

| 错误 | 原因 | 解决 |
|-----|------|------|
| `Invalid tissue 'xxx'` | 组织类型不对 | 用 general / hematopoietic / cardiovascular / hepatic / renal / neurological |
| `variants list is empty` | 输入为空 | 检查 JSON 是否解析正确 |
| `Failed to write TSV` | 输入字段缺失 | 确保必填字段 CHROM/POS/REF/ALT/GENE 存在 |
| `dgra_core.py exited with code 1` | 核心脚本执行失败 | 看 stderr 输出排查 |
| `Offline mode: no cached data` | 离线模式但基因未归档 | 先在线运行一次建立归档，或换在线模式 |

---

## 📁 文件结构

```
dgra-genomic-risk/
  SKILL.md                  # 本文件
  config.json               # 元数据配置
  scripts/
    dgra_cli_wrapper.py     # ⭐ 推荐入口：agent 调用此 wrapper
    dgra_core.py            # 核心分析引擎入口（向后兼容，~200行）
    gpa_pipeline.py         # Pipeline 主流程（v0.10.0 God Module 拆分）
    gpa_tier_classifier.py  # 三级分级 + 证据链（v0.10.0）
    gpa_report.py           # Markdown/JSON 报告生成（v0.10.0）
    gpa_phaser.py           # 相位分析（v0.10.0）
    gpa_multi_hit.py        # 多基因命中检测（v0.10.0）
    gpa_qc.py               # QC 检查（v0.10.0）
    gpa_two_phase.py        # 两阶段管线优化（v0.10.1）
    gpa_preflight.py        # Preflight 健康检查（v0.10.1）
    gpa_workflow.py         # Workflow-as-Code 定义（v0.11.0）
    dgra_api.py             # API 查询层
    dgra_cache.py           # SQLite 缓存
    dgra_config.py          # 配置管理
    dgra_build_state.py     # 构建状态持久化（v0.6.1）
    gpa_phenotype_match.py  # LLM 表型语义匹配（v0.7.0 Phase 2）
    gpa_i18n.py             # 中英文 consequence 术语映射（v0.7.1 Phase 2）
    dgra_variant_filter.py  # 预过滤模块：strict/clinical/broad（v0.7.1 Phase 3）
    gpa_vcf_annotator.py    # Raw VCF → VEP REST API 实时注释（v0.9.0）
    gpa_transcript_selector.py # 疾病感知转录本选择器（v0.9.0）
  references/
    tissue_context.json     # 组织上下文配置
    dgra.yaml               # 运行时参数配置
    offline_data/           # 离线归档（自动创建）
  cache/
    dgra_cache.db           # API 响应缓存

**Rescue 模块**（GPA Tier 1 为空或表型不匹配时自动触发）：
- `scripts/gpa_gene_set_builder.py` — 动态基因集构建（OMIM + HPO）
- `scripts/gpa_phenotype_rescue.py` — VCF 候选变异救援搜索
```

---

## 🩺 临床使用注意

### Tier 分级含义（按场景动态调整）

**Germline / 疾病遗传风险场景：**
- **Tier 1** = 必须干预的致病突变（如纯合截断、已知致病突变）
- **Tier 2** = 需知情同意并持续监测（如杂合携带者、药物代谢多态性）
- **Tier 3** = 记录归档，不影响当前决策

**Somatic / 肿瘤驱动场景：**
- **Tier 1** = 核心驱动突变（TSG 功能丧失、癌基因热点突变、OncoKB Oncogenic）
- **Tier 2** = 可能驱动突变（Likely Oncogenic、亚克隆突变、药物代谢相关）
- **Tier 3** = 乘客突变 / 无功能影响 / 胚系多态混入

**多基因命中**需确认相位：cis（同一条染色体）风险更高，trans（两条染色体）通常为复合杂合

**GPA 是辅助决策工具，最终临床决策需结合完整临床评估。**

---

## 🧪 测试方案

GPA 采用 **pytest** 框架，测试架构为 **L0~L6 分层 + E2E**，覆盖全部 ~30 个模块。详细测试计划见 `tests/TEST_PLAN.md`。

### 快速开始

```bash
# 安装开发依赖
cd tests && pip install -r requirements-dev.txt

# 运行全部测试
pytest

# 运行指定分层
pytest -m l2          # 单元测试
pytest -m l3          # 集成测试
pytest -m "l2 or l3"  # 单元+集成

# 运行指定优先级
pytest -m p0          # 仅 P0（关键路径）

# 运行纯 mock（无网络，最快）
pytest -m "mock and not recording"

# 录制-回放模式（外部 API）
pytest -m recording

# 重新录制 API 响应（需要网络）
pytest -m recording --record-mode=refresh

# 覆盖率报告
pytest --cov=scripts --cov-report=html
```

### 测试分层

| 分层 | 说明 | 用例数 |
|------|------|--------|
| **L0** | 契约测试：输入/输出 Schema 验证 | ~20 |
| **L1** | 静态测试：模块导入、数据结构、循环依赖 | ~20 |
| **L2** | 单元测试：单模块独立测试（纯 mock） | ~140 |
| **L3** | 集成测试：模块间交互、Pipeline 端到端 | ~20 |
| **L4** | 性能测试：基准、内存、缓存吞吐 | ~30 |
| **L5** | 边界测试：极端输入、畸形数据、编码、并发 | ~30 |
| **L6** | 回归测试：版本间向后兼容性 | ~16 |
| **E2E** | 端到端：完整临床场景 | ~16 |

### 覆盖率目标

- **核心模块**（`gpa_tier_classifier.py`, `gpa_pipeline.py`, `gpa_report.py`, `dgra_api.py`, `gpa_vcf_annotator.py`, `dgra_cli_wrapper.py`）：**≥ 80%**
- **其他模块**：**≥ 60%**

### 录制-回放机制

外部 API（Ensembl/UniProt/GTEx/gnomAD 等）采用录制-回放模式：
- **首次运行**：真实调用 API，响应保存到 `tests/recording/`
- **后续运行**：从录制文件加载，秒级完成
- **API 变更检测**：`--record-mode=refresh` 重新录制，diff 发现 schema 变更

### 关键文件

| 文件 | 说明 |
|------|------|
| `tests/TEST_PLAN.md` | 完整测试计划文档 |
| `tests/pytest.ini` | pytest 配置（markers、覆盖率、超时） |
| `tests/requirements-dev.txt` | 开发依赖（pytest、pytest-asyncio、pytest-cov） |
| `tests/conftest.py` | pytest fixtures + 录制-回放基础设施 + Mock 工具 |
| `tests/recording/` | API 录制响应存储目录 |
