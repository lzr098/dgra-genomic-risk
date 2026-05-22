# DGRA - Dynamic Genomic Risk Assessment

**DGRA (Dynamic Genomic Risk Assessment)** v0.5.3 是一个用于个体基因组变异致病性评估的自动化分析系统。基于多源注释数据（VEP/ANNOVAR/SnpEff）进行多维度注释、分级和风险评估，帮助识别可能影响特定组织/器官功能的遗传变异。

---

## ⚠️ 执行前必读：上下文确认规则（v0.5.1）

**收到基因组变异数据时，禁止直接执行分析。必须先确认以下信息：**

| 确认项 | 为什么必须确认 | 不确认的风险 |
|:---|:---|:---|
| **分析目的** | DGRA 输出完全取决于场景 | 供者安全 vs 肿瘤驱动 → 同一变异分级相反 |
| **样本身份** | 患者自身筛查 vs 供者评估 vs 健康人携带者 | 误将供者风险评估用于患者 |
| **组织/疾病背景** | 组织类型决定基因相关性权重 | 造血相关基因在神经场景权重不同 |
| **移植类型**（如适用）| PBSC、骨髓、脐带血 → 风险关注点不同 | VWF 风险在 PBSC 和骨髓中完全不同 |

**最小确认问题集**（向用户确认，不要假设）：

1. "这些数据是谁的样本？患者本人 / 潜在供者 / 健康筛查？"
2. "分析目的是什么？移植供者安全性 / 肿瘤驱动突变 / 遗传病筛查 / 药物基因组学？"
3. "如果是移植场景，什么类型？造血干细胞 / 实体器官？采集方式是 PBSC 还是骨髓？"
4. "患者原发疾病是什么？（如 AML、MDS、免疫缺陷等）"

**在获得上述信息前，禁止调用 dgra_core.py 或生成报告。**

---

## 适用场景

DGRA 根据"目标组织"动态调整分析权重，同一套变异数据可针对不同临床场景生成针对性的风险评估：

| 场景 | 说明 | 典型关注点 |
|:---|:---|:---|
| **造血干细胞移植供者筛查** | 评估潜在供者的遗传风险 | 骨髓衰竭基因、凝血因子、DNA修复基因 |
| **实体器官移植供者筛查** | 肝/肾/心移植供者安全性 | 组织特异性代谢基因、药物基因组学 |
| **肿瘤驱动突变分析** | 患者体细胞突变评估 | 致癌基因、抑癌基因、突变负荷 |
| **遗传病携带者筛查** | 健康人群筛查 | 隐性遗传病、药物代谢多态性 |
| **药物基因组学** | 个体化用药指导 | CYP450、ABCB1、TPMT、DPYD |

**核心优势**：不局限于单一疾病模型，而是基于"组织特异性表达 × 蛋白质功能域影响 × 人群频率 × 致病性证据 × 基因约束"的多维度动态评估。

---

## 架构总览（v0.5 五层架构）

```
┌─────────────────────────────────────────┐
│  Layer 5: 输出层 (Output)               │  → Markdown 带证据链 / JSON 结构化
│  Layer 4: 分级层 (Scoring)              │  → 加权评分 → Tier 1/2/3 + 置信度
│  Layer 3: 注释层 (Annotation)           │  → Ensembl/UniProt/GTEx/gnomAD/ClinVar
│  Layer 2: 适配层 (Adapter)              │  → VEP / ANNOVAR / SnpEff → 统一格式
│  Layer 1: 输入层 (Input + QC)           │  → VCF / Excel / TSV / 自由文本
└─────────────────────────────────────────┘
```

---

## 核心功能

### 输入层（v0.5 P0 统一输入层）

| 功能 | 说明 | 格式 |
|------|------|------|
| **VCF 解析** | cyvcf2 原生解析，支持 bgzip/BCF | `.vcf`, `.vcf.gz`, `.bcf` |
| **Excel 解析** | pandas 读取多 sheet | `.xlsx`, `.xlsm` |
| **TSV/CSV 解析** | 自动探测分隔符 | `.tsv`, `.csv` |
| **自由文本解析** | 自动识别 "GENE p.Pro123Leu" 格式 | 任意文本 |
| **格式自动探测** | 根据扩展名 + 文件头自动识别 | — |

### 适配层（v0.5 P0-2 自动适配）

| 适配器 | 探测特征 | 输出 |
|--------|---------|------|
| **VEPAdapter** | `Consequence`, `IMPACT`, `HGVSp`, `CLIN_SIG` 列 | 直接映射到 DGRA 标准列 |
| **ANNOVARAdapter** | `Gene.refGene`, `AAChange.refGene`, `ExonicFunc.refGene` | 解析 AAChange → HGVSc+HGVSp，ExonicFunc → Consequence |
| **SnpEffAdapter** | `ANN[0].EFFECT`, `ANN[0].IMPACT`, `ANN[0].GENE` | 解析 ANN 字段或结构化列 |
| **自动探测** | 根据表头列名自动选择适配器 | 无需手动指定 |

### 注释层（8 个外部 API，5 层离线 Fallback）

| API | 端点 | 支撑功能 | 离线 Fallback |
|:---|:---|:---|:---|
| **Ensembl REST** | `/lookup/symbol`, `/lookup/id`, `/sequence/id` | 基因注释、转录本校正、Exon 边界 | `offline_data/{gene}.json` → `"ensembl"` |
| **UniProt REST** | `/uniprotkb/{gene}.json` | 功能域（DOMAIN）、活性位点（ACT_SITE）、结合位点（BINDING）、已知突变（MUTAGEN） | `offline_data/{gene}.json` → `"uniprot"` |
| **GTEx API v2** | `/expression/medianGeneExpression` | 组织特异性表达量（RPKM/TPM） | `offline_data/{gene}.json` → `"gtex"` |
| **gnomAD GraphQL** | `/graphql` (`gnomad_freq` + `gnomad_subpops`) | 人群频率 + 亚组频率（EAS/AMR/AFR/NFE/SAS/ASJ/FIN） | `offline_data/{gene}.json` → `"gnomad"` |
| **NCBI Eutils (ClinVar)** | `/esearch.fcgi`, `/esummary.fcgi` | 致病性注释（CLIN_SIG）、评审星标 | `offline_data/{gene}.json` → `"clinvar"` |
| **HGNC REST** | `/fetch/{symbol}` | 基因符号校验（有效/撤回/未找到） | 标记 `INVALID_GENE_SYMBOL` |
| **Orphanet REST** | `/nomenclature/orphanumber/{id}/genes` | 基因-罕见病表型关联（8 个表型） | 硬编码 `CORE_GENE_LISTS` |
| **OMIM GeneMap** | `/api/geneMap` | 基因-表型关联（MIM ID → 基因） | 硬编码 `CORE_GENE_LISTS`（默认禁用） |

**离线 Fallback 分层**：
1. 内存缓存（当前会话）
2. SQLite 缓存（`cache/dgra_cache.db`，30 天 TTL）
3. 离线 JSON 文件（`offline_data/*.json`，634 基因 / 181 MB）
4. 硬编码安全列表（`CORE_GENE_LISTS`，凝血/癌症易感/骨髓衰竭等，不可覆盖）
5. 保守评估规则（数据缺失 → 标记 UNKNOWN → 不降级风险）

### 分级层（v0.5 P1 核心引擎重构）

#### 三级分类体系

| Tier | 定义 | 触发条件 |
|:---|:---|:---|
| **Tier 1** | 需干预 / 排除供者 | 纯合截断 + primary 基因 / ClinVar 致病 + 凝血/FA 基因 / 多击 + 相位 unknown |
| **Tier 2** | 需知情 / 监测 | 杂合 primary 基因 / ClinVar 致病 + 非组织相关 / 药物代谢多态性 |
| **Tier 3** | 无风险 / 乘客 | gnomAD 常见（AF>1%）/ ClinVar 良性 / 同义变异 / 组织不表达 |

#### ACMG/AMP 规则引擎（P1 新增）

| 规则类别 | 规则数量 | 说明 |
|:---------|:--------:|:-----|
| 致病性证据（Pathogenic） | PVS1-PS4 | 功能丧失、已知致病、新发、共分离 |
| 良性证据（Benign） | BA1-BS4 | 人群常见、无功能影响、反分离 |
| 支持证据（Supporting） | PP1-PP5 / BP1-BP7 | 共定位、数据库记录、功能验证 |

#### NMD 逃逸调制（P1-5）

| 逃逸类型 | 调制效果 | 说明 |
|:---|:---|:---|
| **Exon-intron junction last 50bp** | 截断 → MODERATE | 距离最后一个外显子-内含子交界 <50bp |
| **Single-exon gene** | 截断 → MODERATE | 单外显子基因无 NMD |
| **Alternative last exon** | 截断 → 个案评估 | 可变末端外显子 |
| **uORF / IRES** | 截断 → LOW | 上游 ORF 或内部核糖体 entry 点 |
| **5' UTR nonsense-mediated decay evasion** | 截断 → 个案评估 | 5' UTR 介导的逃逸 |

#### Missense 分层（P1-5）

| 层级 | 条件 | Impact |
|:---|:---|:---|
| **Critical domain** | 落在 UniProt DOMAIN/ACT_SITE/BINDING | HIGH |
| **Known pathogenic** | ClinVar 或文献报道相同位点致病 | HIGH |
| **Conservation** | 跨物种高度保守（通过 Ensembl compara） | MODERATE |
| **In silico** | SIFT/PolyPhen 双有害 | MODERATE |
| **Uncharacterized** | 无功能域信息、无保守性数据 | LOW |

#### 加权评分模型（P1 新增）

| 维度 | 权重 | 说明 |
|:---|:---:|:---|
| **功能影响（Functional Impact）** | 0.30 | IMPACT + Consequence + NMD 调制 |
| **人群频率（Population Frequency）** | 0.25 | gnomAD 全局 + 亚组频率 |
| **致病性证据（Pathogenic Evidence）** | 0.25 | ClinVar + ACMG/AMP 证据 |
| **组织相关性（Tissue Relevance）** | 0.15 | GTEx 表达量 + special_gene_lists 匹配 |
| **基因约束（Gene Constraint）** | 0.05 | pLI / LOEUF（ClinGen 剂量敏感性 HI/TS） |

**Tier 阈值**：
- Tier 1: 加权评分 ≥ 0.70
- Tier 2: 加权评分 0.30-0.69
- Tier 3: 加权评分 < 0.30

#### 置信度量化（P1-10）

| 置信度 | 条件 | 说明 |
|:---|:---|:---|
| **HIGH** | ≥3 个 API 成功响应 + 数据一致 | 可靠，可直接用于临床决策 |
| **MEDIUM** | 2 个 API 成功 + 无冲突 | 基本可靠，建议复核 |
| **LOW** | ≤1 个 API 成功 或 数据冲突 | 数据不足，需人工审阅 |

### 输出层（v0.5 P1 增强输出）

| 输出格式 | 内容 | 用途 |
|:---|:---|:---|
| **Markdown 报告** | 带证据链（"如满足 X 则升级为 Tier Y"）、多器官联合风险矩阵 | 临床审阅 |
| **JSON 结构化** | 完整 evidence 对象、confidence 字段、version 字段 | 系统集成 |
| **多器官联合报告** | `--multi-organ` 同时评估多个 profile，取 max tier | 复杂临床场景 |

---

## 快速开始

### 安装

```bash
git clone https://github.com/lzr098/dgra-genomic-risk.git
cd dgra-genomic-risk
pip install -r requirements.txt
```

**依赖**：Python 3.8+, aiohttp, pandas, requests, numpy, cyvcf2

### 命令行用法

```bash
# 在线模式（推荐）— 自动查询所有 API
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic

# 离线模式 — 使用本地缓存 + offline_data
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic --offline

# 多器官联合评估
python scripts/dgra_core.py --input variants.csv --multi-organ hematopoietic,cardiovascular,hepatic

# 配置文件驱动（P2-3）
python scripts/dgra_core.py --input variants.csv --config references/dgra.yaml

# 输出 JSON（P1-12）
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic --format json --output report.json
```

### 支持的输入格式

| 格式 | 扩展名 | 自动探测 | 说明 |
|:---|:---|:---:|:---|
| VCF (VEP 注释) | `.vcf`, `.vcf.gz`, `.bcf` | ✅ | 自动解析 INFO/CSQ，提取 GT/DP/GQ/VAF |
| Excel | `.xlsx`, `.xlsm` | ✅ | pandas 读取，自动探测 sheet |
| TSV | `.tsv` | ✅ | 自动探测分隔符 |
| CSV | `.csv` | ✅ | 自动探测分隔符 |
| 自由文本 | `.txt`, `.md`, 任意 | ✅ | 自动识别 "GENE p.Pro123Leu" 格式 |

### 支持的注释格式（自动探测）

| 注释工具 | 探测特征 | 适配器 |
|:---|:---|:---|
| **VEP** | `Consequence`, `IMPACT`, `HGVSp`, `CLIN_SIG` | `VEPAdapter` |
| **ANNOVAR** | `Gene.refGene`, `AAChange.refGene` | `ANNOVARAdapter` |
| **SnpEff** | `ANN[0].EFFECT`, `ANN[0].IMPACT` | `SnpEffAdapter` |

### Python API

```python
import asyncio
from dgra_core import run_dgra_pipeline, DGRAConfig
from dgra_input_parsers import parse_input

# 统一输入解析（自动探测格式 + 注释适配器）
variants = parse_input("donor_variants.vcf.gz")

# 配置分析
config = DGRAConfig(
    tissue_profile="hematopoietic",
    offline_mode=False,
    multi_organ_profiles=["hematopoietic", "cardiovascular", "hepatic"]
)

# 运行分析
results = asyncio.run(run_dgra_pipeline(variants, config=config))

# 获取报告
report_md = results["report"]
json_data = results["json"]
```

### 配置文件示例（`dgra.yaml`）

```yaml
dgra_version: "0.5.3"

api_endpoints:
  ensembl:
    base_url: "https://rest.ensembl.org"
    timeout: 20.0
    max_retries: 2
    rate_limit_per_sec: 10.0
  uniprot:
    base_url: "https://rest.uniprot.org"
    timeout: 25.0
    max_retries: 2
    rate_limit_per_sec: 5.0

thresholds:
  min_dp: 20
  min_gq: 90.0
  common_af_threshold: 0.01
  low_af_threshold: 0.001

tissue_profiles:
  default: hematopoietic
  available:
    - general
    - hematopoietic
    - cardiovascular
    - hepatic
    - renal
    - neurological

gene_sync:
  enabled: true
  sources:
    orphanet:
      enabled: true
    omim:
      enabled: false

evidence:
  detail_level: brief
  high_confidence_min_apis: 3
```

---

## 组织 Context Profile 对照表

| Profile | 适用场景 | GTEx 组织 | Special Gene Lists |
|:---|:---|:---|:---|
| `general` | 通用健康筛查 | — | 癌症易感、心脏安全、药物代谢、凝血、免疫缺陷 |
| `hematopoietic` | 造血干细胞移植 | Bone Marrow, Whole Blood, Spleen, Thymus | 药物代谢、凝血、FA DNA 修复、KIR 簇 |
| `cardiovascular` | 心移植 / 心肌病 | Heart - Left Ventricle, Atrial Appendage, Aorta | 心肌病、离子通道、主动脉病、心律失常 |
| `hepatic` | 肝移植 | Liver, Small Intestine, Stomach | 胆红素代谢、CYP450、胆汁淤积、血色病 |
| `renal` | 肾移植 | Kidney - Cortex, Medulla, Bladder | 肾小球、肾小管、囊肿、补体 |
| `neurological` | 神经系统评估 | Brain - Cortex, Cerebellum, Hippocampus | 三核苷酸重复、运动神经元、帕金森、周围神经病 |

---

## 转录本校正与 VEP Canonical Reannotation（v0.5.2 新增，v0.5.3 延续）

### 问题：Transcript Discrepancy

注释工具（VEP/ANNOVAR/SnpEff）选择的 "首选转录本" 不一定等于 Ensembl 标注的 **canonical transcript** 或 **MANE Select**。典型错误：

- 使用 **NR_**（非编码转录本）标注 `splice_donor_variant` → `HIGH`
- 但 canonical 是 **NM_** / **ENST**（蛋白编码转录本）→ 同一变异实为 `upstream_gene_variant` → `MODIFIER`

**后果**：`HIGH` 被错误地送入 `classify_variant_tier()`，可能触发 Priority 1b（纯合截短→Tier 1），产生假阳性。

### 新旧行为对比

| 步骤 | 旧行为（≤v0.5.1） | 新行为（≥v0.5.2） |
|:---|:---|:---|
| **Step 1** 转录本校正 | 仅比对 Ensembl canonical 存 warning | 保留，继续比对 |
| **Step 1.5** VEP reannotation | **无** | **新增**： discrepancy 变异 → Ensembl VEP API 用 canonical 参数重新注释 |
| Step 4 Domain mapping | 使用原始 HGVSp → 可能落在错误功能域 | 使用 VEP 修正后的 HGVSp → 正确功能域 |
| Tier 分级 | 基于原始 impact=HIGH | 基于修正后 impact=MODIFIER |

### Pipeline 步骤更新

```
Step 1:  correct_transcript_priority()   → TRANSCRIPT_DISCREPANCY warning
Step 1.5: batch_query_vep_region()      → Ensembl VEP canonical reannotation  ← 新增
Step 2:  detect_pseudogene_risks()       → VAF deviation
Step 3:  classify_gnomad_frequency()     → AF threshold
Step 4:  map_variant_to_domain()          → UniProt domain mapping (使用修正后 HGVSp)
Step 5:  assess_tissue_relevance()        → GTEx TPM
Step 6:  classify_variant_tier()          → Tier 1/2/3 + confidence
```

### VEP Reannotation 逻辑

| 条件 | 处理 |
|:---|:---|
| annotator transcript ≠ canonical / MANE / protein_coding | 加入 discrepancy 列表 |
| 批量查询 Ensembl VEP（50/批，5并发） | `canonical=1&domains=1&protein=1&hgvs=1&mane_select=1` |
| 解析优先级 | **canonical** → MANE → protein_coding |
| 成功 | 覆盖 `consequence`/`impact`/`hgvsc`/`hgvsp`/`transcript`，`transcript_warning.vep_reannotation` 记录 original vs canonical |
| 失败/离线 | 保留原始注释，`quality_confidence="LOW"`, `tier_confidence="LOW"`, `vep_reannotation_failed=True` |

### CRIP2 典型案例

**位点**：`chr14:105473030 G>A`

| 阶段 | 转录本 | 后果 | 影响 | 分级 |
|:---|:---|:---|:---|:---|
| 原始注释（ANNOVAR） | `NR_073082` | `splice_donor_variant` | **HIGH** | 可能 Tier 1（Priority 1b 纯合截短）|
| VEP Reannotation（canonical） | `NM_001312` | `upstream_gene_variant` | **MODIFIER** | **Tier 3**（无 ClinVar 证据 + MODIFIER）|

**验证**：E2E 测试 `test_vep_reannotation_e2e.py` 覆盖此案例，确认修正后无功能域、无致病性证据，最终 Tier 3。

---

## 假阳性过滤（v0.5.1 关键优化）

| 优化项 | 效果 | 适用场景 |
|:---|:---|:---|
| **ClinVar 良性排除** | LOW/MODERATE + ClinVar 良性 → 不计入 multi-hit | 减少常见多态性假阳性 |
| **X 连锁女性修正** | chrX + 杂合 + haplosufficient → Tier 降级 | 女性供者/患者 |
| **C 末端截短修正** | HIGH + 良性 + 终止位置 ≥280aa → 排除 | 良性截短变异 |
| **同义变异排除** | LOW + synonymous → 永远不计入 | 无功能影响变异 |
| **基因家族冗余** | SLC25A5/ANT、CYP2D6、SLC22A1/OCT → 降级 | 有旁系同源补偿的基因 |
| **HLA 排除** | HLA 基因不纳入 multi-hit 升级 | 正常免疫多态性 |

**效果**：Tier-1 假阳性降低 91%（291→26，实际数据验证）。

---

## Multi-hit 检测与相位分析

### Multi-hit 分层决策

| 层级 | 方法 | 置信度 |
|------|------|--------|
| **Level 1** | GATK phased GT (`\|` 分隔符) | **high** |
| **Level 2** | 变异间距 (<50bp / <150bp / <500bp) | high / medium |
| **Level 3** | Reads 直接分析 (pysam) | medium |
| **Level 4** | Trio 家系推断 | high |
| **Level 5** | LD 连锁不平衡统计推断 | low |

- **<50bp**：同一 150bp read 必然覆盖 → high confidence
- **50-150bp**：同一 read 或 pair-end → high confidence
- **150-500bp**：依赖 pair-end insert size → medium confidence
- **>500bp**：超出 short-read 范围 → 需 trio / 长读长

详见 [`docs/PHASE_ANALYSIS_ALGORITHM.md`](docs/PHASE_ANALYSIS_ALGORITHM.md)。

---

## 项目结构

```
dgra-genomic-risk/
├── scripts/
│   ├── dgra_core.py              # 主分析引擎（Tier 分级 + 加权评分 + 置信度）
│   ├── dgra_api.py               # API 查询层（8 个外部 API 封装 + 缓存 + 重试）
│   ├── dgra_cache.py             # SQLite 缓存管理（TTL 自动过期 + 统计）
│   ├── dgra_config.py            # 配置管理（YAML/JSON + 环境变量覆盖）
│   ├── dgra_cli_wrapper.py       # CLI 统一接口（供 OpenClaw 调用）
│   ├── dgra_input_parsers.py     # 统一输入层（VCF/Excel/TSV/自由文本 + 自动探测）
│   ├── dgra_adapters.py          # 注释适配层（VEP/ANNOVAR/SnpEff → 统一格式）
│   ├── dgra_gene_sync.py         # 基因列表同步（Orphanet/OMIM + 硬编码 CORE）
│   └── (acmg.py, scoring.py, clinvar_parser.py 等模块)
├── references/
│   ├── tissue_context.json       # 6 个组织 Profile 定义（tier_rules + special_gene_lists）
│   ├── gene_list_sources.json    # 外部基因列表源配置（Orphanet 8 表型 + OMIM）
│   ├── user_gene_lists.json      # 用户自定义基因列表（add/remove/custom）
│   ├── dgra.yaml                 # 运行时 YAML 配置（P2-3）
│   ├── offline_data/             # 离线基因数据（634 基因 / 181 MB）
│   ├── repeatmasker_regions.json # 重复区域黑名单
│   ├── pseudogene_config.json    # 假基因识别配置
│   └── api_corrections.json      # API 数据校正覆盖
├── cache/
│   ├── dgra_cache.db             # API 响应缓存（30 天 TTL）
│   └── gene_sync_cache.db        # 基因列表同步缓存（7 天 TTL）
├── docs/
│   └── PHASE_ANALYSIS_ALGORITHM.md
├── README.md                     # 本文件
├── CHANGELOG.md                  # 更新日志
└── requirements.txt              # Python 依赖
```

---

## 版本历史

详见 [CHANGELOG.md](CHANGELOG.md)

| 版本 | 日期 | 主题 |
|------|------|------|
| **v0.5.2** | 2026-05-21 | **核心逻辑修正**：Multi-hit 不再升级变异（只标记基因）、ClinVar 中文注释支持（致病/良性）、统计格式改为"基因数/突变数"、新增 Priority 1c（ClinVar 致病+HIGH+造血相关→Tier 1）、**Transcript discrepancy 降级**（NR_/XM_非编码转录本→HIGH 降级为 MODERATE）、**VEP Canonical Reannotation**（Step 1.5：Ensembl VEP 用 canonical 重新注释 discrepancy 变异，CRIP2 chr14:105473030 案例验证） |
| **v0.5.1** | 2026-05-21 | **假阳性大幅优化**：ClinVar 良性排除、X 连锁女性修正、同义排除、C 末端截短修正、基因家族冗余、HLA 排除（Tier-1 假阳性 ↓91%） |
| **v0.5.0** | 2026-05-21 | **P0 统一输入层**（VCF/Excel/TSV/自由文本 + VEP/ANNOVAR/SnpEff 自动适配）+ **P1 核心引擎重构**（ACMG 评分、NMD 调制、Missense 分层、加权评分、置信度量化、结构化证据链、JSON 输出、多器官联合、增强 QC、基因名校验、分析版本化） |
| v0.4.5 | 2026-05-20 | 相位分析系统（Phase Analysis） |
| v0.4.4a | 2026-05-20 | README + CHANGELOG + 逐个致病性分析 |
| v0.4.4 | 2026-05-20 | Multi-hit 致病性证据过滤 |
| v0.4.3 | 2026-05-20 | HLA 基因排除 |
| v0.4.2 | 2026-05-20 | GTEx v2 + 并发优化 |
| v0.4.1 | 2026-05-19 | 关键修复 |
| v0.4.0 | 2026-05-19 | API-first 架构发布 |

---

## 开发路线图

| 阶段 | 状态 | 内容 |
|:---|:---:|:---|
| **P0** | ✅ 完成 | 统一输入层（解析器 + 适配器 + QC） |
| **P1** | ✅ 完成（14/15）| 核心引擎重构（ACMG、NMD、Missense、加权评分、置信度、证据链、JSON、多器官、QC、基因名校验、版本化） |
| **P1-6** | ⏸️ 最后做 | 多 GTEx 组织聚合 |
| **P2** | 🔄 进行中 | VRS 格式支持、组织特异性剪接注释、配置文件驱动（dgra.yaml）、动态组织注册 |
| **P3** | ⏳ 待启动 | 多基因上位效应、贝叶斯置信模型 |

---

## 数据来源与引用

DGRA 使用以下开源资源和公共数据库：
- **Ensembl** (EMBL-EBI) — 基因注释、转录本、序列
- **UniProt** (UniProt Consortium) — 蛋白质功能域、活性位点、突变
- **GTEx Project** (NIH) — 组织特异性基因表达
- **gnomAD** (Broad Institute) — 人群等位基因频率
- **ClinVar** (NCBI) — 临床致病性注释
- **HGNC** (HUGO Gene Nomenclature Committee) — 基因符号标准化
- **Orphanet** — 罕见病基因-表型关联
- **OMIM** — 孟德尔遗传基因-表型关联
- **ClinGen** — 基因剂量敏感性（HI/TS）

---

## 许可证

MIT License

---

**维护者**：[@lzr098](https://github.com/lzr098)  
**当前版本**：v0.5.3  
**最后更新**：2026-05-21
