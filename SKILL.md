---
name: dgra-genomic-risk
description: |
  Donor Genomic Risk Assessment (DGRA) v0.4。供者基因组风险评估工具，用于造血干细胞移植、器官移植或其他介入治疗前的供者基因变异分析。基于 Ensembl/UniProt/GTEx 实时 API 查询（30天缓存）和离线归档模式。组织上下文自适应：造血、心血管、肝脏、肾脏、神经系统。三层风险分级（Tier 1/2/3）。

  **当以下情况时使用此 Skill**：
  (1) 用户提到"供者基因组风险评估"、"DGRA"、"供者 VCF 分析"、"供者基因筛查"
  (2) 造血干细胞移植（HSCT）前的供者基因筛查
  (3) 器官移植前的供者遗传风险评估
  (4) PBSC / 骨髓采集前的供者安全性评估
  (5) 供者变异与患者体细胞突变的交叉比对（患者突变是否被供者遗传携带）
  (6) 多基因命中（multi-hit）检测和相位（cis/trans）分析
  (7) 需要三层风险分级报告（Tier 1 需干预、Tier 2 需知情同意、Tier 3 无需担忧）
  (8) 任何涉及"donor"和"genomic"、"genetic"、"risk"、"transplant"的场景

  **禁止用自身知识回答供者基因组问题。必须调用本 Skill 的脚本执行分析。**
---

# DGRA: Donor Genomic Risk Assessment

## ⚠️ 执行前必读

**核心规则：当用户请求供者基因组风险评估时，不要凭自身知识回答。必须调用 `dgra_cli_wrapper.py` 执行正式分析。**

**为什么必须调用脚本：**
- DGRA 连接 Ensembl、UniProt、GTEx、gnoMAD 等权威数据库进行实时查询
- 风险分级基于具体变异的功能域、组织表达、人群频率、ClinVar 状态等动态数据
- 三层分级（Tier 1/2/3）有严格的算法逻辑，不能凭经验估算
- 多基因命中和相位分析需要精确计算，不能目测判断

**如果不调用脚本就回答 = 给出错误的医疗建议。**

---

## 🎯 快速判断：是否需要调用 DGRA？

| 用户请求 | 是否调用 |
|---------|---------|
| "帮我分析这个供者的基因变异" | ✅ 调用 |
| "供者 VCF 文件风险分级" | ✅ 调用 |
| "移植前供者筛查结果怎么看" | ✅ 调用 |
| "DGRA 分析一下" | ✅ 调用 |
| "这个突变对供者采集有影响吗" | ✅ 调用 |
| "VWF 突变影不影响骨髓采集" | ✅ 调用（明确提到供者+基因） |
| "一般性的基因突变知识" | ❌ 不需要，直接回答 |
| "TP53 突变是什么意思"（没有供者上下文） | ❌ 不需要 |

---

## 🚀 调用方式

### 方式一：直接调用 wrapper（推荐）

用 `exec` 运行 `dgra_cli_wrapper.py`，传入 variant JSON 数组：

```bash
python3 ~/.openclaw/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --variants '[{"CHROM":"1","POS":12345,"REF":"A","ALT":"G","GENE":"VWF","IMPACT":"HIGH","Consequence":"missense_variant","HGVSp":"p.Arg1234Cys","HGVSc":"c.3700C>T","CLIN_SIG":"Pathogenic","GT":"0/1","DP":30,"GQ":99,"VAF":0.5}]' \
  --tissue hematopoietic
```

### 方式二：已有 TSV 文件

如果用户已提供 TSV/CSV 文件路径：

```bash
python3 ~/.openclaw/skills/dgra-genomic-risk/scripts/dgra_core.py \
  --input /path/to/donor_variants.tsv \
  --tissue hematopoietic \
  --output /tmp/dgra_report.md \
  --json /tmp/dgra_results.json
```

### 方式三：含患者突变交叉比对

```bash
python3 ~/.openclaw/skills/dgra-genomic-risk/scripts/dgra_cli_wrapper.py \
  --variants '[...]' \
  --tissue hematopoietic \
  --patient-mutations '[{"gene":"BCOR","hgvsp":"p.Arg1234*","impact":"HIGH"}]'
```

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

**如果用户提供了 VCF 或表格数据，先提取这些字段构造 JSON 数组，再传给 wrapper。**

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
      "multi_hit_genes": ["VWF"],
      "patient_inherited_mutations": []
    },
    "tier1_variants": [],
    "tier2_variants": [...],
    "tier3_variants": [...],
    "multi_hit_details": [...],
    "patient_donor_cross_check": [...],
    "report_markdown": "# DGRA 报告..."
  },
  "report_md": "# DGRA Report...",
  "stdout": "DGRA Report Generated..."
}
```

### 给用户呈现的关键信息

1. **风险分级统计**：Tier 1 / Tier 2 / Tier 3 各多少个
2. **高风险变异详情**（Tier 1 和 Tier 2）：基因、突变、影响、建议行动
3. **多基因命中**：是否有同一基因多个变异，相位状态（cis/trans/unknown）
4. **患者-供者交叉比对**：患者突变是否被供者遗传携带
5. **Markdown 报告**：完整报告文本可直接呈现给用户

### 组织类型选择

| 场景 | 组织类型 |
|------|---------|
| 造血干细胞移植 / PBSC / 骨髓采集 | `hematopoietic` |
| 心脏移植 / 供心评估 | `cardiovascular` |
| 肝脏移植 | `hepatic` |
| 肾脏移植 | `renal` |
| 神经系统移植 / 神经介入 | `neurological` |

---

## 🔧 离线模式

当网络不可用或 API 超时频繁时，添加 `--offline` 参数：

```bash
python3 .../dgra_cli_wrapper.py --variants '[...]' --tissue hematopoietic --offline
```

离线模式使用本地缓存（`references/offline_data/` 下的基因 JSON），对于已有归档的基因结果与在线模式一致。未归档的基因 fallback 到保守规则。

---

## ❌ 常见错误

| 错误 | 原因 | 解决 |
|-----|------|------|
| `Invalid tissue 'xxx'` | 组织类型不对 | 用 hematopoietic / cardiovascular / hepatic / renal / neurological |
| `variants list is empty` | 输入为空 | 检查 JSON 是否解析正确 |
| `Failed to write TSV` | 输入字段缺失 | 确保必填字段 CHROM/POS/REF/ALT/GENE 存在 |
| `dgra_core.py exited with code 1` | 核心脚本执行失败 | 看 stderr 输出排查 |
| `Offline mode: no cached data` | 离线模式但基因未归档 | 先在线运行一次建立归档，或换在线模式 |

---

## 📁 文件结构

```
dgra-genomic-risk/
  SKILL.md                  # 本文件
  scripts/
    dgra_cli_wrapper.py     # ⭐ 推荐入口：agent 调用此 wrapper
    dgra_core.py            # 核心分析引擎（async API-first）
    dgra_api.py             # API 查询层
    dgra_cache.py           # SQLite 缓存
    dgra_config.py          # 配置管理
  references/
    tissue_context.json     # 组织上下文配置
    offline_data/           # 离线归档（自动创建）
  cache/
    dgra_cache.db           # API 响应缓存
```

---

## 🩺 临床使用注意

- **Tier 1** = 必须干预或排除供者（如凝血功能障碍影响采集安全）
- **Tier 2** = 需知情同意并术后监测（如携带者状态）
- **Tier 3** = 记录归档，不影响决策
- **多基因命中**需确认相位：cis（同一条染色体）风险更高，trans（两条染色体）通常为复合杂合
- **患者-供者交叉比对**：检查患者体细胞驱动突变是否被供者遗传携带，影响移植后复发风险评估

**DGRA 是辅助决策工具，最终临床决策需结合完整临床评估。**
