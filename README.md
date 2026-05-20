# DGRA - 供者基因组风险评估系统

**DGRA (Donor Genomic Risk Assessment)** 是一个用于造血干细胞移植供者基因组风险评估的分析系统。通过整合多源生物信息学数据库，对外显子组测序（WES）数据进行三级分类，帮助临床医生识别供者基因组中可能影响移植安全性的变异。

---

## 核心功能

| 功能 | 说明 |
|------|------|
| **多源数据库整合** | Ensembl（基因注释）、UniProt（蛋白质功能域）、GTEx（组织表达）、gnomAD（人群频率） |
| **动态组织背景** | 支持造血系统、神经系统、心血管系统等多种临床场景 |
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

### 基本用法

```bash
# 在线模式（推荐）
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic

# 离线模式（无网络）
python scripts/dgra_core.py --input variants.csv --tissue hematopoietic --offline
```

### Python API

```python
import asyncio
from dgra_core import run_dgra_pipeline, DGRAConfig

config = DGRAConfig(tissue_profile="hematopoietic")
results = asyncio.run(run_dgra_pipeline(variants_data, config=config))

# 获取报告
report_md = results["report"]
json_data = results["json"]
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
- 可能影响移植期用药剂量和效果
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
│   ├── dgra_api.py        # API 查询层
│   ├── dgra_cache.py      # 缓存管理
│   └── dgra_config.py     # 配置管理
├── references/
│   └── offline_data/      # 离线数据存档
├── cache/
│   └── dgra_cache.db      # SQLite 缓存
├── config.json            # 用户配置
└── CHANGELOG.md           # 更新日志
```

---

## 版本历史

详见 [CHANGELOG.md](CHANGELOG.md)

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
**当前版本**：v0.4.4
