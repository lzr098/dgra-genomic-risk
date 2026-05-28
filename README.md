<h1 align="center">
  <code>GPA</code> · Genomic Phenotype Association
</h1>

<p align="center">
  <strong>个体基因组变异 → 表型关联 → 三级风险分级</strong><br>
  API-first · 离线容灾 · 组织感知 · 证据链可溯
</p>


---

## 一句话

GPA 接收任何格式的基因组变异数据（VCF / Excel / TSV / 自由文本），通过 8 个公共 API 实时注释 + 5 层离线容灾，按组织特异性动态加权评分，输出带完整证据链的 Tier 1/2/3 分级报告。

**和直接看 VEP/ClinVar 的区别**：GPA 不只是"翻译"注释——它把功能影响 × 人群频率 × 致病性证据 × 组织表达 × 基因约束 五个维度加权融合，同一变异在不同临床场景下可以分级完全不同。

---

## 核心数据

| 指标 | 数值 |
|------|------|
| 核心代码 | 16,669 行（20 个模块） |
| 测试代码 | 5,605 行（3 层：单元 + 集成 + E2E） |
| 外部 API | 8 个（Ensembl / UniProt / GTEx / gnomAD / ClinVar / HGNC / Orphanet / OMIM） |
| 假基因对 | 51 个临床相关对 |
| Tier-1 假阳性下降 | **↓ 91%**（v0.4→v0.5 实测：291→26） |
| 输入格式 | VCF / VCF.gz / BCF / Excel / TSV / CSV / 自由文本 |
| 注释适配 | VEP / ANNOVAR / SnpEff 自动探测 |

---

## 快速开始

```bash
git clone https://github.com/lzr098/dgra-genomic-risk.git
cd dgra-genomic-risk
pip install -r requirements.txt   # aiohttp, cyvcf2, openpyxl
```

### 30 秒跑通

```bash
# 在线模式：110K 变异的 VCF → 自动分批 → 约 3 分钟出报告
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue hematopoietic \
  --filter-preset clinical

# 离线模式：无网也能跑（本地归档 + 缓存）
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue general \
  --offline

# Raw VCF（无注释）：自动调 VEP REST API 实时注释
python scripts/dgra_cli_wrapper.py \
  --input-file raw_variants.vcf \
  --tissue neurological \
  --disease-description "acute myeloid leukemia"

# 开启 SpliceAI 剪接预测（默认关闭）
python scripts/dgra_cli_wrapper.py \
  --input-file variants.tsv \
  --tissue cardiovascular \
  --spliceai

# 多器官联合评估
python scripts/dgra_cli_wrapper.py \
  --input-file variants.tsv \
  --multi-organ hematopoietic,cardiovascular,hepatic
```

---

## 架构

```
┌──────────────────────────────────────────────────────┐
│                    Output Layer                       │
│  Markdown (证据链) · JSON (结构化) · 多器官联合报告     │
├──────────────────────────────────────────────────────┤
│                  Scoring Layer                        │
│  Tier 1/2/3 加权评分 · ACMG 证据 · SpliceAI 调制      │
│  NMD 逃逸 · Missense 5 层 · ClinVar 星级置信度         │
├──────────────────────────────────────────────────────┤
│                Annotation Layer                       │
│  8 API 实时查询 · SQLite 缓存(30d) · 5 层离线容灾      │
├──────────────────────────────────────────────────────┤
│                 Adapter Layer                         │
│  VEP / ANNOVAR / SnpEff → 统一内部格式                 │
├──────────────────────────────────────────────────────┤
│               Input + QC + Filter                     │
│  VCF·Excel·TSV·自由文本 · strict/clinical/broad       │
│  Raw VCF → VEP REST 实时注释 · 疾病感知转录本选择       │
└──────────────────────────────────────────────────────┘
```

详细设计文档：[DESIGN_v0.4.md](DESIGN_v0.4.md)

---

## 关键能力

### 🧬 组织感知分级

同一变异，不同临床场景，分级可能完全不同：

| 基因 | 变异 | hematopoietic | cardiovascular | neurological |
|------|------|:---:|:---:|:---:|
| RUNX1 | frameshift | **Tier 1** | Tier 2 | Tier 3 |
| MYH11 | splice_region | Tier 2 | **Tier 1** | Tier 3 |
| HTT | CAG repeat | Tier 3 | Tier 3 | **Tier 1** |

6 个组织 Profile：`general` · `hematopoietic` · `cardiovascular` · `hepatic` · `renal` · `neurological`

### 🛡️ 5 层离线容灾

```
在线查询 → 内存缓存 → SQLite(30d TTL) → 离线归档 → 硬编码安全列表 → 保守规则
```

断网/限流/被封 IP → GPA 仍能产出报告（置信度标记为 LOW），不会卡死。离线归档覆盖高频查询基因，确保关键变异在无网环境下仍可评估。

### 🔬 SpliceAI 剪接验证（v0.8.0+）

VEP 标注 `splice_donor_variant` / `HIGH` 不一定真影响剪接。SpliceAI delta score 作为独立证据：

- delta = 0 → **降级**（VEP 过调用，假阳性剪接）
- delta ≥ 0.5 → **升级**（强剪接变化，MODERATE → Tier 2+）
- 双通道：Broad Institute API（主） + Ensembl VEP REST SpliceAI plugin（fallback）

### 🧪 Raw VCF 端到端（v0.9.0+）

丢进来一个完全没注释的 VCF → GPA 自动检测 → VEP REST API 实时注释 → 疾病感知转录本选择 → 分析出报告。**零预处理。**

转录本选择四层评分：组织表达 → 后果严重性 → canonical/MANE → 位置偏好。歧义时支持 LLM 辅助。

### ⚡ 假阳性治理

| 版本 | 机制 | Tier-1 FP |
|------|------|:---------:|
| v0.4.3 | HLA 排除 | 334→238 |
| v0.4.4 | Multi-hit 致病证据过滤 | 238→28 |
| v0.5.1 | ClinVar 良性排除 + X 连锁修正 + 同义排除 | 291→**26** |
| v0.5.2 | VEP Canonical Reannotation（NR_→NM_ 修正） | — |
| v0.7.1 | ClinVar 冲突注释检测 + 预过滤 | — |
| v0.7.2 | ClinVar 星级置信度（single_submitter=0.40） | — |
| v0.9.1 | gnomAD 频率守卫（AF>1% → Tier 3） | — |

### 🏗️ 工程稳定性

- **指数退避重试**：HTTP 429/502/503/504 + Timeout → 1s→2s→4s
- **流式下载 + 断点续传**：GTF 大文件 chunk 8KB + HTTP Range 206
- **构建状态持久化**：`.dgra_build_state.json` 原子化步骤记录，崩溃可恢复
- **自动分批**：>500 变异自动分批，每批独立超时
- **Proxy 自适应**：`global_config.proxy` 统一推导 `trust_env`，不再硬编码
- **大数据集自动路由**：>5000 变异强制走 VEP API 而非本地 subprocess

---

## 输入格式

| 格式 | 扩展名 | 自动探测 | 说明 |
|:---|:---|:---:|:---|
| VCF (VEP 注释) | `.vcf` `.vcf.gz` `.bcf` | ✅ | 解析 INFO/CSQ，提取 GT/DP/GQ/VAF |
| Raw VCF (无注释) | `.vcf` | ✅ | 自动触发 VEP REST API 实时注释 |
| Excel | `.xlsx` `.xlsm` | ✅ | pandas 读取，自动探测 sheet |
| TSV / CSV | `.tsv` `.csv` | ✅ | 自动探测分隔符 |
| 自由文本 | `.txt` `.md` 任意 | ✅ | 识别 "GENE p.Pro123Leu" 格式 |

注释工具自动适配：VEP · ANNOVAR · SnpEff

---

## CLI 参数速查

```
必选：
  --input-file PATH          输入文件路径
  --tissue PROFILE           组织类型：general|hematopoietic|cardiovascular|hepatic|renal|neurological

可选：
  --offline                  离线模式（缓存 + 归档，不查询 API）
  --filter-preset PRESET     预过滤：strict|clinical|broad（默认 clinical）
  --spliceai                 启用 SpliceAI 剪接预测
  --somatic                  肿瘤体细胞模式
  --multi-organ P1,P2,...    多器官联合评估
  --disease-description TXT  疾病描述（触发疾病感知转录本选择）
  --output-json PATH         输出 JSON 文件
  --config PATH              YAML 配置文件
  --batch-size N             分批大小（默认 500）
  --timeout N                每批超时秒数（默认 300）
```

---

## Python API

```python
import asyncio
from dgra_core import run_gpa_pipeline, GPAConfig
from dgra_input_parsers import parse_input

# 自动探测格式 + 注释适配
variants = parse_input("patient_variants.vcf.gz")

# 配置
config = GPAConfig(
    tissue_profile="hematopoietic",
    filter_preset="clinical",
    spliceai_enabled=True,
    disease_description="acute myeloid leukemia",
)

# 运行
results = asyncio.run(run_gpa_pipeline(variants, config=config))
report_md = results["report"]
```

---

## 项目结构

```
dgra-genomic-risk/
├── scripts/
│   ├── dgra_core.py                # 核心引擎（Variant/GPAConfig/dataclass）
│   ├── gpa_pipeline.py             # Pipeline 主流程（async 编排）
│   ├── gpa_tier_classifier.py      # 三级分级 + 证据链
│   ├── gpa_report.py               # Markdown/JSON 报告生成
│   ├── gpa_vcf_annotator.py        # Raw VCF → VEP REST 注释
│   ├── gpa_transcript_selector.py  # 疾病感知转录本选择器
│   ├── gpa_phenotype_match.py      # LLM 语义表型匹配
│   ├── gpa_phaser.py               # 相位分析（GATK/distance/trio）
│   ├── gpa_multi_hit.py            # 多基因命中检测
│   ├── gpa_qc.py                   # QC 检查
│   ├── gpa_i18n.py                 # 中英文术语映射
│   ├── gpa_input.py                # 输入格式探测
│   ├── dgra_api.py                 # 8 API 封装 + 重试 + 缓存
│   ├── dgra_adapters.py            # VEP/ANNOVAR/SnpEff 适配
│   ├── dgra_input_parsers.py       # 统一输入解析
│   ├── dgra_cli_wrapper.py         # CLI 入口
│   ├── dgra_cache.py               # SQLite 缓存管理
│   ├── dgra_config.py              # YAML/JSON 配置
│   ├── dgra_splice_predictor.py    # SpliceAI 查询（Broad + VEP REST）
│   ├── dgra_variant_filter.py      # 预过滤模块
│   ├── dgra_gene_sync.py           # 基因列表同步
│   ├── dgra_myvariant.py           # MyVariant.info 批量查询
│   ├── dgra_batch_runner.py        # 批量运行器
│   └── dgra_pseudogene_sync.py     # GENCODE 假基因同步
├── references/
│   ├── tissue_context.json          # 6 组织 Profile 定义
│   ├── dgra.yaml                    # 运行时配置
│   ├── gene_phenotype_map.json      # 基因-表型映射
│   ├── pseudogene_lookup.json       # 51 假基因对
│   └── offline_data/                # 离线查询归档
├── tests/                           # 5,605 行测试
│   ├── test_a_layer.py              # A-Layer 稳定性回归
│   ├── test_v09.py                  # v0.9 功能测试
│   ├── test_v091.py                 # v0.9.1 热修复测试
│   ├── test_v095_hotfix.py          # v0.9.5 热修复测试
│   ├── test_v0_9_3_e2e.py           # v0.9.3 E2E
│   ├── test_spliceai.py             # SpliceAI 模块测试
│   └── ...                          # L1~L5 分层测试
├── docs/
│   └── PHASE_ANALYSIS_ALGORITHM.md
├── README.md
├── CHANGELOG.md
├── DESIGN_v0.4.md
└── requirements.txt
```

---

## 迭代时间线

**8 天，10 个版本，0 个周末停更。**

```
v0.4.0  ──→  v0.5.0  ──→  v0.6.0  ──→  v0.7.0  ──→  v0.8.0  ──→  v0.9.0  ──→  v0.10.0
5/19       5/21        5/22        5/23        5/23        5/23        5/25
API-first   P0/P1       假基因       表型关联      SpliceAI     Raw VCF     God Module
架构发布     引擎重构     架构升级     LLM语义      剪接验证      端到端       拆分重构
```

| 版本 | 日期 | 主题 |
|------|------|------|
| **v0.10.0** | 2026-05-25 | **God Module 拆分**：dgra_core.py 2098行 → 6 个独立模块（pipeline/tier/report/phaser/multi_hit/qc） |
| **v0.9.5** | 2026-05-24 | **P0/P1 热修**：MAD2L2/OR2B11 假阳性、SpliceAI API 404 迁移、pairwise 相位 Bug |
| **v0.9.3** | 2026-05-23 | **P0 热修**：gnomAD GraphQL schema 变更、VEP 批量串行→并发、proxy 自适应 |
| **v0.9.2** | 2026-05-23 | **P0 热修**：gnomAD chr 前缀缺失、variant_id 误判 |
| **v0.9.1** | 2026-05-23 | **P0 热修**：DDX3X 常见多态误判 Tier 1（三 Bug 致命链） |
| **v0.9.0** | 2026-05-23 | **Raw VCF 端到端**：VEP REST 实时注释 + 疾病感知转录本选择 |
| **v0.8.0** | 2026-05-23 | **SpliceAI 剪接预测**：Broad API + VEP REST fallback + delta 0/0.5 降升级 |
| **v0.7.2** | 2026-05-23 | **ClinVar 星级置信度**：practice_guideline 0.95 / single_submitter 0.40 |
| **v0.7.1** | 2026-05-23 | **预过滤 + 中文兼容**：strict/clinical/broad + gpa_i18n |
| **v0.7.0** | 2026-05-23 | **表型关联**：LLM 语义匹配 + Top 100 罕见病基因 |
| **v0.6.0** | 2026-05-22 | **假基因架构**：51 对 + VAF 模式检测 + confidence 降级原则 |
| v0.5.x | 2026-05-21 | **P0/P1 交付**：统一输入层 + ACMG + NMD + 加权评分 + 假阳性↓91% |
| v0.4.0 | 2026-05-19 | **API-first 架构发布** |

---

## 组织 Profile 对照表

| Profile | 适用场景 | GTEx 组织 | Special Gene Lists |
|:---|:---|:---|:---|
| `general` | 通用健康筛查 | — | 癌症易感、心脏安全、药物代谢、凝血、免疫缺陷 |
| `hematopoietic` | 造血/血液肿瘤 | Bone Marrow, Whole Blood, Spleen, Thymus | 药物代谢、凝血、FA DNA 修复、KIR 簇 |
| `cardiovascular` | 心血管/心肌病 | Heart, Aorta | 心肌病、离子通道、主动脉病、心律失常 |
| `hepatic` | 肝脏 | Liver, Small Intestine | 胆红素代谢、CYP450、胆汁淤积、血色病 |
| `renal` | 肾脏 | Kidney Cortex/Medulla, Bladder | 肾小球、肾小管、囊肿、补体 |
| `neurological` | 神经系统 | Brain Cortex, Cerebellum, Hippocampus | 三核苷酸重复、运动神经元、帕金森、周围神经病 |

---

## 数据来源

Ensembl · UniProt · GTEx · gnomAD · ClinVar · HGNC · Orphanet · OMIM · ClinGen

---

## 许可证

MIT-0

---

**维护者**：[@lzr098](https://github.com/lzr098)
**当前版本**：v0.10.0
**最后更新**：2026-05-26
