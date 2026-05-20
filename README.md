# DGRA - 个体基因组风险评估系统

**DGRA (Dynamic Genomic Risk Assessment)** 是一个用于个体基因组变异致病性评估的自动化分析系统。基于 VEP 注释后的 VCF/CSV 数据进行多维度注释和分级，帮助识别可能影响特定组织/器官功能的遗传变异。

---

## 适用场景

DGRA 根据"目标组织"动态调整分析权重，同一套 VCF 数据可针对不同临床场景生成针对性的风险评估：

- **血液系统** — 造血干细胞、红细胞、免疫细胞相关基因
- **神经系统** — 脑、脊髓、神经元相关基因  
- **心血管系统** — 心肌、血管内皮相关基因
- **肝脏** — 肝细胞代谢、解毒相关基因
- **肾脏** — 肾小管、肾小球相关基因
- **药物基因组学** — 药物代谢酶、转运蛋白相关基因

**核心优势**：不局限于单一疾病模型，而是基于"组织特异性表达 × 蛋白质功能域影响 × 人群频率 × 致病性证据"的多维度动态评估。

---

## 核心功能

| 功能 | 说明 |
|------|------|
| **多源数据库整合** | Ensembl（基因注释）、UniProt（蛋白质功能域）、GTEx（组织表达）、gnomAD（人群频率）、ClinVar（致病性） |
| **动态组织背景** | 支持造血、神经、心血管、肝脏、肾脏等多种临床场景 |
| **三级分类体系** | Tier 1（需关注）→ Tier 2（需知情）→ Tier 3（无风险） |
| **Multi-hit 检测** | 识别同一基因的多个致病性变异，评估复合杂合风险 |
| **相位分析** | 基于 GATK GT 格式和变异间距判断 cis/trans 关系 |
| **离线模式** | 支持无网络环境下的分析，使用本地缓存数据 |

---

## 快速开始

### 安装

```bash
git clone https://github.com/lzr098/dgra-genomic-risk.git
cd dgra-genomic-risk
pip install -r requirements.txt
```

**依赖**：Python 3.8+, aiohttp, pandas, requests, numpy

### 基本用法

#### 命令行

```bash
# 在线模式（推荐）— 自动查询 Ensembl/UniProt/GTEx
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic

# 离线模式（无网络）— 使用本地缓存
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic --offline

# 自定义配置
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic \
    --config custom_config.json --output report.md
```

#### 支持的 tissue 参数

| tissue 值 | 适用场景 |
|-----------|----------|
| `hematopoietic` | 造血系统 / 血液疾病 |
| `neural` | 神经系统 / 神经退行性疾病 |
| `cardiovascular` | 心血管系统 / 心肌病 |
| `hepatic` | 肝脏 / 代谢疾病 |
| `renal` | 肾脏 / 肾病 |
| `pulmonary` | 肺部 / 呼吸系统 |

#### 输入文件格式

DGRA 接受 VEP 注释后的 CSV 文件，必须包含以下列：

```csv
CHROM,POS,REF,ALT,GENE,Feature,EXON,IMPACT,Consequence,HGVSp,HGVSc,CLIN_SIG,DP,GQ,GT,VAF,gnomAD_AF
```

对于中文标注的 VEP 输出，DGRA 会自动映射列名。

#### Python API

```python
import asyncio
from dgra_core import run_dgra_pipeline, DGRAConfig

# 准备变异数据（从 VEP 输出解析）
variants_data = [
    {
        "CHROM": "chr1",
        "POS": 123456,
        "REF": "A",
        "ALT": "G",
        "GENE": "BRCA1",
        "Feature": "NM_007294",
        "IMPACT": "MODERATE",
        "Consequence": "missense_variant",
        "HGVSp": "NP_009225.1:p.Ser746Asn",
        "HGVSc": "c.2237G>A",
        "CLIN_SIG": "Pathogenic",
        "GT": "0/1",
        "VAF": 0.45,
        "gnomAD_AF": 0.0001
    }
]

# 配置分析
config = DGRAConfig(
    tissue_profile="hematopoietic",  # 组织背景
    offline_mode=False               # 在线/离线模式
)

# 运行分析
results = asyncio.run(run_dgra_pipeline(variants_data, config=config))

# 获取报告
report_md = results["report"]
json_data = results["json"]

# 保存报告
with open("report.md", "w") as f:
    f.write(report_md)
```

### 配置文件

```json
{
  "tissue_profile": "hematopoietic",
  "offline_mode": false,
  "cache_ttl_days": 30
}
```

---

## 三级分类说明

### Tier 1 - 需关注（Action Required）

需要临床关注的变异：
- ClinVar 致病性/可能致病性（Pathogenic / Likely Pathogenic）
- 影响关键蛋白质功能域 + 目标组织表达
- 同一基因多个致病性变异（multi-hit，需确认 cis/trans 相位）
- 剪切位点改变

**处理建议**：临床确认、家族史调查、可能需针对性监测

### Tier 2 - 需知情（Inform & Monitor）

应知情的变异：
- 药物代谢相关（CYP 家族、ABCB1 等）
- 可能影响治疗期用药剂量和效果
- 不构成本身疾病风险，但需医生知情

**处理建议**：告知医生调整用药方案

### Tier 3 - 无风险（No Concern）

无需关注的变异：
- ClinVar Benign / Likely Benign
- 目标组织中不表达（GTEx TPM < 1.0）
- 不影响功能域的常见多态性

---

## Multi-hit 相位分析

DGRA v0.4.5 新增相位分析系统，用于判断同一基因内多个致病性变异是否位于同一等位基因（cis）或不同等位基因（trans）。

### 分层决策逻辑

| 层级 | 方法 | 置信度 | 适用条件 |
|------|------|--------|----------|
| **Level 1** | GATK phased GT (`\|` 分隔符) | **high** | GT 明确区分单倍型 |
| **Level 2** | 变异间距 (<50bp / <150bp / <500bp) | high / medium | 短 reads 物理覆盖可行性 |
| **Level 3** | Reads 直接分析 (pysam) | medium | 需要 BAM/CRAM |
| **Level 4** | Trio 家系推断 | high | 需要父母数据 |
| **Level 5** | LD 连锁不平衡统计推断 | low | 人群数据支持 |

### 间距与可行性映射

- **<50bp**：同一 150bp read 必然覆盖 → high confidence
- **50-150bp**：同一 read（靠近 3' 端）或 pair-end → high confidence
- **150-500bp**：依赖 pair-end insert size → medium confidence
- **>500bp**：超出 short-read 范围 → 需 trio / 长读长

详见 [`docs/PHASE_ANALYSIS_ALGORITHM.md`](docs/PHASE_ANALYSIS_ALGORITHM.md)。

---

## 数据来源

| 数据库 | 用途 | 版本 |
|--------|------|------|
| **Ensembl REST API** | 基因注释、转录本校正 | v1.0 |
| **UniProt REST API** | 蛋白质功能域映射 | 当前 |
| **GTEx API v2** | 组织特异性表达 | v8 |
| **gnomAD** | 人群等位基因频率 | v2.1/v4 |
| **ClinVar** | 致病性注释 | 当前 |

---

## 项目结构

```
dgra-genomic-risk/
├── scripts/
│   ├── dgra_core.py       # 主分析引擎
│   ├── dgra_api.py        # API 查询层（Ensembl/UniProt/GTEx）
│   ├── dgra_cache.py      # SQLite 缓存管理
│   └── dgra_config.py     # 配置管理
├── references/
│   └── offline_data/      # 离线数据存档（JSON）
├── cache/
│   └── dgra_cache.db      # SQLite 缓存数据库
├── config.json            # 用户配置文件
├── README.md              # 本文件
├── CHANGELOG.md           # 更新日志
├── docs/
│   └── PHASE_ANALYSIS_ALGORITHM.md  # 相位分析算法文档
└── requirements.txt       # Python 依赖
```

---

## 分析流程

1. **数据解析**：读取 VEP 注释的 CSV，提取关键字段
2. **转录本校正**：Ensembl REST API → canonical transcript
3. **假基因检测**：VAF 偏差分析识别已知假基因对
4. **gnomAD 频率分类**：常见（AF>1%）/ 罕见（AF<0.1%）/ 未捕获
5. **蛋白质功能域映射**：UniProt REST API → DOMAIN/REGION 特征
6. **组织相关性评估**：GTEx API → 目标组织表达量（TPM）
7. **三级分类**：综合 ClinVar、Impact、Domain、Tissue 信息
8. **Multi-hit 检测**：同一基因的多个致病性变异
9. **相位分析**：GATK GT + 间距判断 → cis/trans 相位状态
10. **报告生成**：Markdown 格式，包含详细变异列表和致病性分析

---

## 版本历史

详见 [CHANGELOG.md](CHANGELOG.md)

| 版本 | 日期 | 主题 |
|------|------|------|
| v0.4.5 | 2026-05-20 | 相位分析系统（Phase Analysis） |
| v0.4.4a | 2026-05-20 | README + CHANGELOG + 逐个致病性分析 |
| v0.4.4 | 2026-05-20 | Multi-hit 致病性证据过滤 |
| v0.4.3 | 2026-05-20 | HLA 基因排除 |
| v0.4.2 | 2026-05-20 | GTEx v2 + 并发优化 |
| v0.4.1 | 2026-05-19 | 关键修复 |
| v0.4.0 | 2026-05-19 | API-first 架构发布 |

---

## 引用

DGRA 使用以下开源资源和公共数据库：
- Ensembl (EMBL-EBI)
- UniProt (UniProt Consortium)
- GTEx Project (NIH)
- gnomAD (Broad Institute)
- ClinVar (NCBI)

---

## 许可证

MIT License

---

**维护者**：@lzr098  
**当前版本**：v0.4.5  
**最后更新**：2026-05-20
