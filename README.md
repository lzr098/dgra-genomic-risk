# DGRA - 供者基因组风险评估系统

**DGRA (Donor Genomic Risk Assessment)** 是一个用于供者基因组风险评估的自动化分析系统。通过对全外显子组测序（WES）数据进行多维度注释和分级，帮助临床医生识别供者基因组中可能影响移植安全性或治疗效果的遗传变异。

---

## 适用场景

DGRA 不仅限于造血干细胞移植，可应用于任何需要供者基因组评估的场景：

- **器官移植**（肾、肝、心、肺等）— 评估供者遗传病风险
- **造血干细胞移植** — 评估供者造血相关基因变异
- **组织移植** — 评估供者对特定组织的遗传风险
- **生殖细胞/胚胎捐赠** — 评估遗传病携带状态
- **生物样本库** — 样本提供者的遗传风险评估

**核心优势**：根据"目标组织"动态调整分析权重，同一套 WES 数据，换用不同组织背景可得到针对性的风险报告。

---

## 核心功能

| 功能 | 说明 |
|------|------|
| **多源数据库整合** | Ensembl（基因注释）、UniProt（蛋白质功能域）、GTEx（组织表达）、gnomAD（人群频率）、ClinVar（致病性） |
| **动态组织背景** | 支持造血系统、神经系统、心血管系统、肝脏、肾脏等多种临床场景 |
| **三级分类体系** | Tier 1（需关注）→ Tier 2（需知情）→ Tier 3（无风险） |
| **Multi-hit 检测** | 识别同一基因的多个致病性变异，评估复合杂合风险 |
| **患者-供者交叉核对** | 检查患者体细胞驱动突变是否存在于供者基因组 |
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
| `hematopoietic` | 造血系统 / 造血干细胞移植 |
| `neural` | 神经系统 / 脑组织 |
| `cardiovascular` | 心血管系统 / 心脏 |
| `hepatic` | 肝脏 / 肝移植 |
| `renal` | 肾脏 / 肾移植 |
| `pulmonary` | 肺部 / 肺移植 |

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

需要临床干预的变异：
- ClinVar 致病性/可能致病性（Pathogenic / Likely Pathogenic）
- 影响关键蛋白质功能域 + 目标组织表达
- 同一基因多个致病性变异（multi-hit，需确认 cis/trans 相位）
- 剪切位点改变

**处理建议**：移植前确认、家族史调查、可能需更换供者

### Tier 2 - 需知情（Inform & Monitor）

供者应知情的变异：
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
8. **Multi-hit 检测**：同一基因的多个致病性变异 → 相位确认提醒
9. **患者-供者交叉核对**：体细胞驱动突变的遗传检测
10. **报告生成**：Markdown 格式，包含详细变异列表和致病性分析

---

## 版本历史

详见 [CHANGELOG.md](CHANGELOG.md)

| 版本 | 日期 | 主题 |
|------|------|------|
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
**当前版本**：v0.4.4a  
**最后更新**：2026-05-20
