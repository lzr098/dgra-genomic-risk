<h1 align="center">
  <code>GPA</code> · Genomic Phenotype Association
</h1>

<p align="center">
  <a href="https://github.com/lzr098/dgra-genomic-risk"><img src="https://img.shields.io/badge/version-0.10.16-blue" alt="version"></a>
  <a href="#"><img src="https://img.shields.io/badge/python-3.10%2B-green" alt="python"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT--0-lightgrey" alt="license"></a>
</p>

<p align="center">
  <strong>Individual genomic variants → Phenotype association → Tier 1/2/3 risk classification</strong><br>
  API-first · Offline resilient · Tissue-aware · Evidence-traceable
</p>

<p align="center">
  <a href="#中文文档">中文</a> · <a href="#english-documentation">English</a>
</p>

---

<h2 id="中文文档">📖 中文文档</h2>

### 一句话

GPA 接收任何格式的基因组变异数据（VCF / Excel / TSV / 自由文本），通过 8 个公共 API 实时注释 + 5 层离线容灾，按组织特异性动态加权评分，输出带完整证据链的 Tier 1/2/3 分级报告。

**和直接看 VEP/ClinVar 的区别**：GPA 不只是"翻译"注释——它把功能影响 × 人群频率 × 致病性证据 × 组织表达 × 基因约束 五个维度加权融合，同一变异在不同临床场景下可以分级完全不同。

---

### 📋 目录

- [功能特性](#功能特性)
- [系统要求](#系统要求)
- [部署方式](#部署方式)
  - [WorkBuddy 部署](#workbuddy-部署)
  - [CodeX / OpenClaw 部署](#codex--openclaw-部署)
  - [通用部署（直接 clone）](#通用部署直接-clone)
- [依赖安装](#依赖安装)
- [快速开始](#快速开始)
- [CLI 使用指南](#cli-使用指南)
- [Python API 使用](#python-api-使用)
- [输入格式](#输入格式)
- [组织类型](#组织类型)
- [两阶段管线优化](#两阶段管线优化)
- [离线模式](#离线模式)
- [常见问题](#常见问题)
- [测试](#测试)

---

### 功能特性

| 特性 | 说明 |
|------|------|
| 🧬 **组织感知分级** | 同一变异在不同组织场景下分级不同（9 个 Profile） |
| 🔌 **8 API 实时注释** | Ensembl / UniProt / GTEx / gnomAD / ClinVar / HGNC / Orphanet / OMIM |
| 🛡️ **5 层离线容灾** | 在线 → 内存缓存 → SQLite(30d) → 离线归档 → 硬编码安全列表 → 保守规则 |
| 🔬 **SpliceAI 剪接验证** | Broad Institute API + Ensembl VEP REST fallback |
| 🧪 **Raw VCF 端到端** | 无注释 VCF 自动调用 VEP REST API 实时注释 + 疾病感知转录本选择 |
| ⚡ **两阶段管线优化** | 大型 VCF API 调用量减少 50-200x |
| 🩺 **表型 Rescue 搜索** | 自动分级未发现 Tier 1 时，根据患者表型动态构建基因集，救援被遗漏的候选变异 |
| 🛡️ **Preflight 健康检查** | 分析前自动检查依赖就绪状态 |
| 🌍 **中英文兼容** | 输入支持中文/英文 consequence 术语自动映射 |

---

### 系统要求

| 项目 | 要求 |
|------|------|
| Python | ≥ 3.10 |
| 操作系统 | macOS / Linux / Windows (WSL) |
| 网络 | 在线模式需访问 8 个公共 API（可代理） |
| 磁盘 | ~150MB（代码 + 缓存 + 离线归档 + gnomAD local） |

---

### 部署方式

#### WorkBuddy 部署

WorkBuddy 是 macOS 上的 AI 助手桌面应用，支持 Skill 扩展。

**步骤 1：打开 Skill 目录**

```bash
# 在终端中打开 WorkBuddy Skill 目录
open ~/.workbuddy/skills/
```

**步骤 2：克隆仓库到 Skill 目录**

```bash
cd ~/.workbuddy/skills/
git clone https://github.com/lzr098/dgra-genomic-risk.git
```

**步骤 3：重启 WorkBuddy 或刷新 Skill 列表**

- WorkBuddy 会自动扫描 `~/.workbuddy/skills/` 目录下的 Skill
- 重启应用或等待自动刷新后，GPA Skill 即可使用

**步骤 4：验证部署**

在 WorkBuddy 对话中输入：

```
分析我的基因组变异数据
```

如果 GPA Skill 被正确加载，系统会提示确认分析目的和样本身份。

> 💡 **提示**：WorkBuddy 的 Skill 安装路径为 `~/.workbuddy/skills/`。你也可以通过 WorkBuddy 的 Skill 市场搜索 `gpa-genomic-phenotype` 直接安装。

#### CodeX / OpenClaw 部署

CodeX（或 OpenClaw）是命令行/IDE 集成的 AI 编程助手。

**步骤 1：找到 Skill 目录**

CodeX 的 Skill 目录通常位于：

```bash
# 默认路径
~/.codex/skills/
# 或
~/.openclaw/skills/
```

如果目录不存在，手动创建：

```bash
mkdir -p ~/.codex/skills
```

**步骤 2：克隆仓库**

```bash
cd ~/.codex/skills/
git clone https://github.com/lzr098/dgra-genomic-risk.git
```

**步骤 3：配置 CodeX 识别 Skill**

CodeX 通常通过 `config.json` 自动识别 Skill。确保 `dgra-genomic-risk/config.json` 存在：

```json
{
  "id": "gpa-genomic-phenotype",
  "name": "GPA",
  "version": "0.10.15",
  "description": "GPA (Genomic Phenotype Association) with dynamic tissue-context analysis",
  "license": "MIT-0",
  "requires": {
    "python": ">=3.10",
    "packages": ["aiohttp"]
  },
  "entry": "scripts/dgra_cli_wrapper.py",
  "profiles": ["general", "hematopoietic", "cardiovascular", "hepatic", "renal", "neurological", "endocrine", "metabolic", "ophthalmic"]
}
```

**步骤 4：安装依赖**

```bash
cd dgra-genomic-risk
pip install -r requirements.txt
```

**步骤 5：重启 CodeX**

重启 IDE 或 CodeX 扩展，GPA Skill 即可在对话中使用。

#### 通用部署（直接 clone）

如果你不使用任何 AI 助手平台，也可以独立部署 GPA：

```bash
# 克隆仓库
git clone https://github.com/lzr098/dgra-genomic-risk.git
cd dgra-genomic-risk

# 安装依赖
pip install -r requirements.txt

# 验证安装
python scripts/dgra_cli_wrapper.py --help
```

---

### 依赖安装

**方式一：venv（推荐）**

```bash
python3 -m venv venv
source venv/bin/activate   # Linux/macOS
# 或 venv\Scripts\activate  # Windows

pip install -r requirements.txt
```

**方式二：直接安装（macOS/Linux）**

```bash
pip install -r requirements.txt --break-system-packages
```

**依赖列表**（`requirements.txt`）：

```
aiohttp>=3.9.0      # 异步 HTTP 客户端，API 查询核心依赖
vcfpy               # VCF 文件解析
openpyxl>=3.1.0     # Excel 文件读取
chardet>=5.0.0      # 编码自动探测
```

**可选依赖**：

```bash
# 如果需要 LLM 辅助转录本选择
pip install openai    # 或任意 OpenAI-compatible API 客户端

# 开发测试依赖
pip install -r tests/requirements-dev.txt
```

---

### 快速开始

#### 30 秒跑通

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

#### 预检健康检查（推荐）

首次使用前，建议先运行预检：

```bash
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue hematopoietic \
  --preflight
```

预检会验证：Python 依赖、8 个 API 连通性、本地文件、磁盘空间、网络代理。

---

### 🩺 表型 Rescue 搜索

当 GPA 自动分级**未发现 Tier 1** 或 Tier 1/2 与患者表型不匹配时，启动 Rescue 模块：

```bash
# Step 1: 根据表型动态构建基因集（OMIM + HPO）
python scripts/gpa_gene_set_builder.py \
  --phenotypes "joubert,polydactyly,epilepsy" \
  --omim-db ~/.workbuddy/data/omim/omim.db \
  --output genes.txt \
  --max-genes 80

# Step 2: 在 VCF 中搜索候选变异
python scripts/gpa_phenotype_rescue.py \
  --vcf patient.vep.vcf.gz \
  --gene-list genes.txt \
  --output rescue.tsv \
  --patient-sex male \
  --min-impact MODERATE \
  --max-af 0.01
```

**Rescue 能解决什么问题？**
- MODERATE  impact 变异被自动分级低估
- ClinVar 标签保守（"likely benign" / VUS）的新兴基因
- GTEx 组织数据缺失导致组织评分失败
- X 连锁男性半合子被忽略
- 数据库滞后、尚未收录的新致病基因

**典型场景**：患者有明显临床综合征（如多指 + 小脑蚓部发育不全），GPA 因数据库原因未检出 Tier 1，Rescue 通过 ciliopathy 通路动态构建基因集，在 VCF 中定位到 OFD1 半合子致病变异。

---

### CLI 使用指南

#### 必选参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `--input-file PATH` | 输入文件路径 | `patient.vcf.gz` |
| `--tissue PROFILE` | 组织类型 | `general` / `hematopoietic` / `cardiovascular` / `hepatic` / `renal` / `neurological` / `endocrine` / `metabolic` / `ophthalmic` |

#### 可选参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--offline` | 离线模式（不查询 API） | 否 |
| `--filter-preset PRESET` | 预过滤：`strict` / `clinical` / `broad` | `clinical` |
| `--spliceai` | 启用 SpliceAI 剪接预测 | 否 |
| `--somatic` | 肿瘤体细胞模式 | 否 |
| `--multi-organ P1,P2,...` | 多器官联合评估 | 否 |
| `--disease-description TXT` | 疾病描述（触发疾病感知转录本选择） | 否 |
| `--two-phase` | 启用两阶段管线（大型 VCF 优化） | 否 |
| `--preflight` | 分析前执行健康检查 | 否 |
| `--output-json PATH` | 输出 JSON 文件 | 否 |
| `--config PATH` | YAML 配置文件 | 否 |
| `--report-detail-level LEVEL` | 报告详细程度：`minimal` / `standard` / `full` | `minimal` |
| `--timeout N` | 每批超时秒数 | 300 |

#### 完整示例

```bash
# 生成详细程度可控的报告（minimal / standard / full）
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue general \
  --report-detail-level minimal

# 大型 WES VCF，两阶段优化
python scripts/dgra_cli_wrapper.py \
  --input-file large_wes.vcf.gz \
  --tissue general \
  --two-phase \
  --filter-preset clinical \
  --spliceai \
  --output-json report.json
```

---

### Python API 使用

```python
import asyncio
from scripts.dgra_core import run_gpa_pipeline, GPAConfig
from scripts.dgra_input_parsers import parse_input

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
print(report_md)
```

---

### 输入格式

| 格式 | 扩展名 | 自动探测 | 说明 |
|:---|:---|:---:|:---|
| VCF (VEP 注释) | `.vcf` `.vcf.gz` `.bcf` | ✅ | 解析 INFO/CSQ，提取 GT/DP/GQ/VAF |
| Raw VCF (无注释) | `.vcf` | ✅ | 自动触发 VEP REST API 实时注释 |
| Excel | `.xlsx` `.xlsm` | ✅ | pandas 读取，自动探测 sheet |
| TSV / CSV | `.tsv` `.csv` | ✅ | 自动探测分隔符 |
| 自由文本 | `.txt` `.md` 任意 | ✅ | 识别 "GENE p.Pro123Leu" 格式 |

注释工具自动适配：VEP · ANNOVAR · SnpEff

#### Variant JSON 格式（直接传入）

```python
variants = [
    {
        "CHROM": "12",
        "POS": 6126538,
        "REF": "G",
        "ALT": "A",
        "GENE": "VWF",
        "IMPACT": "HIGH",
        "Consequence": "stop_gained",
        "HGVSp": "p.Gln1311Ter",
        "HGVSc": "c.3931C>T",
        "CLIN_SIG": "Pathogenic",
        "GT": "0/1",
        "DP": 30,
        "GQ": 99,
        "VAF": 0.5
    }
]
```

必填字段：`CHROM`, `POS`, `REF`, `ALT`, `GENE`

---

### 组织类型

| Profile | 适用场景 |
|:---|:---|
| `general` | 通用健康筛查（默认） |
| `hematopoietic` | 造血/血液肿瘤 |
| `cardiovascular` | 心血管/心肌病 |
| `hepatic` | 肝脏 |
| `renal` | 肾脏 |
| `neurological` | 神经系统 |
| `endocrine` | 内分泌系统/代谢 |
| `metabolic` | 代谢疾病 |
| `ophthalmic` | 眼科/视网膜 |

同一变异在不同组织下分级可能完全不同：

| 基因 | 变异 | hematopoietic | cardiovascular | neurological |
|------|------|:---:|:---:|:---:|
| RUNX1 | frameshift | **Tier 1** | Tier 2 | Tier 3 |
| MYH11 | splice_region | Tier 2 | **Tier 1** | Tier 3 |

---

### 两阶段管线优化

针对大型 VCF（>5,000 变异）的 API 调用优化：

```bash
python scripts/dgra_cli_wrapper.py \
  --input-file large_wes.vcf.gz \
  --tissue general \
  --two-phase
```

**Phase 1**（本地，<30 秒）：VEP 注释 + 本地基因列表 → 过滤 >95% 常见 SNP/低影响变异

**Phase 2**（API 仅查候选）：gnomAD、SpliceAI、表型 LLM 仅对 Tier 1/2 候选变异执行

典型 germline VCF 的 API 调用量减少 **50-200x**。

---

### 离线模式

当网络不可用或 API 超时频繁时：

```bash
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue general \
  --offline
```

离线模式使用本地缓存（`references/offline_data/` 下的基因 JSON）。未归档的基因 fallback 到保守规则。

---

### 常见问题

| 问题 | 原因 | 解决 |
|------|------|------|
| `Invalid tissue 'xxx'` | 组织类型不对 | 使用 general / hematopoietic / cardiovascular / hepatic / renal / neurological / endocrine / metabolic / ophthalmic |
| `variants list is empty` | 输入为空 | 检查文件路径和格式 |
| `Failed to write TSV` | 输入字段缺失 | 确保必填字段 CHROM/POS/REF/ALT/GENE 存在 |
| `dgra_core.py exited with code 1` | 核心脚本执行失败 | 看 stderr 输出排查 |
| `Offline mode: no cached data` | 离线模式但基因未归档 | 先在线运行一次建立归档 |
| API 超时频繁 | 网络不稳定 | 使用 `--offline` 或检查代理设置 |

---

### 测试

```bash
# 安装开发依赖
cd tests && pip install -r requirements-dev.txt

# 运行全部测试
pytest

# 运行指定分层
pytest -m l2          # 单元测试
pytest -m l3          # 集成测试

# 覆盖率报告
pytest --cov=scripts --cov-report=html
```

---

<hr>

<h2 id="english-documentation">📖 English Documentation</h2>

### One-Liner

GPA accepts genomic variant data in any format (VCF / Excel / TSV / free text), performs real-time annotation through 8 public APIs with 5-layer offline resilience, dynamically weights scores by tissue specificity, and outputs Tier 1/2/3 classification reports with full evidence chains.

**The difference from reading VEP/ClinVar directly**: GPA doesn't just "translate" annotations—it fuses functional impact × population frequency × pathogenic evidence × tissue expression × gene constraint into a weighted score. The same variant can have completely different tiers in different clinical contexts.

---

### Table of Contents

- [Features](#features)
- [System Requirements](#system-requirements)
- [Deployment](#deployment)
  - [WorkBuddy](#workbuddy-deployment)
  - [CodeX / OpenClaw](#codex--openclaw-deployment)
  - [Generic Deployment](#generic-deployment)
- [Dependency Installation](#dependency-installation)
- [Quick Start](#quick-start)
- [CLI Usage Guide](#cli-usage-guide)
- [Python API Usage](#python-api-usage)
- [Input Formats](#input-formats)
- [Tissue Profiles](#tissue-profiles)
- [Two-Phase Pipeline](#two-phase-pipeline)
- [Offline Mode](#offline-mode)
- [FAQ](#faq)
- [Testing](#testing)

---

### Features

| Feature | Description |
|---------|-------------|
| 🧬 **Tissue-Aware Classification** | Same variant, different tiers across 9 tissue profiles |
| 🔌 **8-API Real-Time Annotation** | Ensembl / UniProt / GTEx / gnomAD / ClinVar / HGNC / Orphanet / OMIM |
| 🛡️ **5-Layer Offline Resilience** | Online → Memory cache → SQLite(30d) → Offline archive → Hardcoded safety list → Conservative rules |
| 🔬 **SpliceAI Splice Verification** | Broad Institute API + Ensembl VEP REST fallback |
| 🧪 **Raw VCF End-to-End** | Unannotated VCF auto-detected → VEP REST API annotation → Disease-aware transcript selection |
| ⚡ **Two-Phase Pipeline** | API calls reduced 50-200x for large VCFs |
| 🩺 **Phenotype Rescue Search** | When automated tiering finds no Tier 1, dynamically build gene set from phenotypes and rescue missed candidates |
| 🛡️ **Preflight Health Check** | Auto-check dependency readiness before analysis |
| 🌍 **Bilingual Support** | Chinese/English consequence term auto-mapping |

---

### System Requirements

| Item | Requirement |
|------|-------------|
| Python | ≥ 3.10 |
| OS | macOS / Linux / Windows (WSL) |
| Network | Online mode requires access to 8 public APIs (proxy supported) |
| Disk | ~150MB (code + cache + offline archive + gnomAD local) |

---

### Deployment

#### WorkBuddy Deployment

WorkBuddy is an AI assistant desktop app for macOS that supports Skill extensions.

**Step 1: Open the Skill directory**

```bash
# Open WorkBuddy Skill directory in terminal
open ~/.workbuddy/skills/
```

**Step 2: Clone the repository into the Skill directory**

```bash
cd ~/.workbuddy/skills/
git clone https://github.com/lzr098/dgra-genomic-risk.git
```

**Step 3: Restart WorkBuddy or refresh the Skill list**

- WorkBuddy auto-scans `~/.workbuddy/skills/` for Skills
- Restart the app or wait for auto-refresh, then GPA Skill is ready

**Step 4: Verify deployment**

In WorkBuddy chat, type:

```
Analyze my genomic variant data
```

If GPA Skill is loaded correctly, the system will prompt to confirm analysis purpose and sample identity.

> 💡 **Tip**: WorkBuddy's Skill installation path is `~/.workbuddy/skills/`. You can also search for `gpa-genomic-phenotype` in WorkBuddy's Skill marketplace to install directly.

#### CodeX / OpenClaw Deployment

CodeX (or OpenClaw) is a command-line/IDE-integrated AI programming assistant.

**Step 1: Locate the Skill directory**

CodeX's Skill directory is typically at:

```bash
# Default path
~/.codex/skills/
# or
~/.openclaw/skills/
```

Create if it doesn't exist:

```bash
mkdir -p ~/.codex/skills
```

**Step 2: Clone the repository**

```bash
cd ~/.codex/skills/
git clone https://github.com/lzr098/dgra-genomic-risk.git
```

**Step 3: Configure CodeX to recognize the Skill**

CodeX typically auto-detects Skills via `config.json`. Ensure `dgra-genomic-risk/config.json` exists:

```json
{
  "id": "gpa-genomic-phenotype",
  "name": "GPA",
  "version": "0.10.15",
  "description": "GPA (Genomic Phenotype Association) with dynamic tissue-context analysis",
  "license": "MIT-0",
  "requires": {
    "python": ">=3.10",
    "packages": ["aiohttp"]
  },
  "entry": "scripts/dgra_cli_wrapper.py",
  "profiles": ["general", "hematopoietic", "cardiovascular", "hepatic", "renal", "neurological", "endocrine", "metabolic", "ophthalmic"]
}
```

**Step 4: Install dependencies**

```bash
cd dgra-genomic-risk
pip install -r requirements.txt
```

**Step 5: Restart CodeX**

Restart the IDE or CodeX extension, and GPA Skill is available in chat.

#### Generic Deployment

If you don't use any AI assistant platform, you can deploy GPA standalone:

```bash
# Clone repository
git clone https://github.com/lzr098/dgra-genomic-risk.git
cd dgra-genomic-risk

# Install dependencies
pip install -r requirements.txt

# Verify installation
python scripts/dgra_cli_wrapper.py --help
```

---

### Dependency Installation

**Option 1: venv (Recommended)**

```bash
python3 -m venv venv
source venv/bin/activate   # Linux/macOS
# or venv\Scripts\activate  # Windows

pip install -r requirements.txt
```

**Option 2: Direct install (macOS/Linux)**

```bash
pip install -r requirements.txt --break-system-packages
```

**Dependencies** (`requirements.txt`):

```
aiohttp>=3.9.0      # Async HTTP client, core API dependency
vcfpy               # VCF file parser
openpyxl>=3.1.0     # Excel file reader
chardet>=5.0.0      # Encoding auto-detection
```

**Optional dependencies**:

```bash
# For LLM-assisted transcript selection
pip install openai    # or any OpenAI-compatible API client

# Development test dependencies
pip install -r tests/requirements-dev.txt
```

---

### Quick Start

#### Get Started in 30 Seconds

```bash
# Online mode: 110K variant VCF → auto-batching → ~3 min report
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue hematopoietic \
  --filter-preset clinical

# Offline mode: works without internet (local archive + cache)
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue general \
  --offline

# Raw VCF (no annotation): auto VEP REST API annotation
python scripts/dgra_cli_wrapper.py \
  --input-file raw_variants.vcf \
  --tissue neurological \
  --disease-description "acute myeloid leukemia"

# Enable SpliceAI splice prediction (default: off)
python scripts/dgra_cli_wrapper.py \
  --input-file variants.tsv \
  --tissue cardiovascular \
  --spliceai

# Multi-organ joint assessment
python scripts/dgra_cli_wrapper.py \
  --input-file variants.tsv \
  --multi-organ hematopoietic,cardiovascular,hepatic
```

#### Preflight Health Check (Recommended)

Before first use, run preflight:

```bash
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue hematopoietic \
  --preflight
```

Preflight verifies: Python dependencies, 8 API connectivity, local files, disk space, network proxy.

---

### 🩺 Phenotype Rescue Search

When GPA automated tiering **finds no Tier 1** or Tier 1/2 variants do not match the patient's phenotype, trigger the Rescue module:

```bash
# Step 1: Dynamically build gene set from phenotypes (OMIM + HPO)
python scripts/gpa_gene_set_builder.py \
  --phenotypes "joubert,polydactyly,epilepsy" \
  --omim-db ~/.workbuddy/data/omim/omim.db \
  --output genes.txt \
  --max-genes 80

# Step 2: Search VCF for candidate variants
python scripts/gpa_phenotype_rescue.py \
  --vcf patient.vep.vcf.gz \
  --gene-list genes.txt \
  --output rescue.tsv \
  --patient-sex male \
  --min-impact MODERATE \
  --max-af 0.01
```

**What problems does Rescue solve?**
- MODERATE impact variants underestimated by automated tiering
- Emerging genes with conservative ClinVar labels ("likely benign" / VUS)
- Tissue scoring failures due to missing GTEx data
- X-linked hemizygosity in males overlooked
- Database lag — newly discovered disease genes not yet annotated

**Typical scenario**: A patient presents with a recognizable syndrome (e.g., polydactyly + cerebellar vermis hypoplasia). GPA finds no Tier 1 due to database limitations. Rescue dynamically builds a ciliopathy gene set, scans the VCF, and identifies a hemizygous OFD1 pathogenic variant.

---

### CLI Usage Guide

#### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--input-file PATH` | Input file path | `patient.vcf.gz` |
| `--tissue PROFILE` | Tissue type | `general` / `hematopoietic` / `cardiovascular` / `hepatic` / `renal` / `neurological` / `endocrine` / `metabolic` / `ophthalmic` |

#### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--offline` | Offline mode (no API queries) | No |
| `--filter-preset PRESET` | Pre-filter: `strict` / `clinical` / `broad` | `clinical` |
| `--spliceai` | Enable SpliceAI splice prediction | No |
| `--somatic` | Tumor somatic mode | No |
| `--multi-organ P1,P2,...` | Multi-organ joint assessment | No |
| `--disease-description TXT` | Disease description (triggers disease-aware transcript selection) | No |
| `--two-phase` | Enable two-phase pipeline (large VCF optimization) | No |
| `--preflight` | Run health check before analysis | No |
| `--output-json PATH` | Output JSON file | No |
| `--config PATH` | YAML config file | No |
| `--report-detail-level LEVEL` | Report detail: `minimal` / `standard` / `full` | `minimal` |
| `--timeout N` | Timeout per batch (seconds) | 300 |

#### Complete Example

```bash
# Generate report with controllable detail level (minimal / standard / full)
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue general \
  --report-detail-level minimal

# Large WES VCF with two-phase optimization
python scripts/dgra_cli_wrapper.py \
  --input-file large_wes.vcf.gz \
  --tissue general \
  --two-phase \
  --filter-preset clinical \
  --spliceai \
  --output-json report.json
```

---

### Python API Usage

```python
import asyncio
from scripts.dgra_core import run_gpa_pipeline, GPAConfig
from scripts.dgra_input_parsers import parse_input

# Auto-detect format + annotation adapter
variants = parse_input("patient_variants.vcf.gz")

# Configure
config = GPAConfig(
    tissue_profile="hematopoietic",
    filter_preset="clinical",
    spliceai_enabled=True,
    disease_description="acute myeloid leukemia",
    report_detail_level="minimal",  # minimal / standard / full
)

# Run
results = asyncio.run(run_gpa_pipeline(variants, config=config))
report_md = results["report"]
print(report_md)
```

---

### Input Formats

| Format | Extensions | Auto-detect | Notes |
|:---|:---|:---:|:---|
| VCF (VEP annotated) | `.vcf` `.vcf.gz` `.bcf` | ✅ | Parse INFO/CSQ, extract GT/DP/GQ/VAF |
| Raw VCF (unannotated) | `.vcf` | ✅ | Auto-trigger VEP REST API annotation |
| Excel | `.xlsx` `.xlsm` | ✅ | pandas read, auto-detect sheet |
| TSV / CSV | `.tsv` `.csv` | ✅ | Auto-detect delimiter |
| Free text | `.txt` `.md` any | ✅ | Recognize "GENE p.Pro123Leu" format |

Annotation tools auto-adapt: VEP · ANNOVAR · SnpEff

#### Variant JSON Format (Pass Directly)

```python
variants = [
    {
        "CHROM": "12",
        "POS": 6126538,
        "REF": "G",
        "ALT": "A",
        "GENE": "VWF",
        "IMPACT": "HIGH",
        "Consequence": "stop_gained",
        "HGVSp": "p.Gln1311Ter",
        "HGVSc": "c.3931C>T",
        "CLIN_SIG": "Pathogenic",
        "GT": "0/1",
        "DP": 30,
        "GQ": 99,
        "VAF": 0.5
    }
]
```

Required fields: `CHROM`, `POS`, `REF`, `ALT`, `GENE`

---

### Tissue Profiles

| Profile | Use Case |
|:---|:---|
| `general` | General health screening (default) |
| `hematopoietic` | Hematopoiesis / hematologic malignancy |
| `cardiovascular` | Cardiovascular / cardiomyopathy |
| `hepatic` | Liver |
| `renal` | Kidney |
| `neurological` | Nervous system |
| `endocrine` | Endocrine system / metabolism |
| `metabolic` | Metabolic diseases |
| `ophthalmic` | Ophthalmology / retina |

Same variant, different tiers across tissues:

| Gene | Variant | hematopoietic | cardiovascular | neurological |
|------|---------|:---:|:---:|:---:|
| RUNX1 | frameshift | **Tier 1** | Tier 2 | Tier 3 |
| MYH11 | splice_region | Tier 2 | **Tier 1** | Tier 3 |

---

### Two-Phase Pipeline

API call optimization for large VCFs (>5,000 variants):

```bash
python scripts/dgra_cli_wrapper.py \
  --input-file large_wes.vcf.gz \
  --tissue general \
  --two-phase
```

**Phase 1** (local, <30s): VEP annotation + local gene list → filter >95% common SNPs / low-impact variants

**Phase 2** (API for candidates only): gnomAD, SpliceAI, phenotype LLM only for Tier 1/2 candidates

Typical germline VCF API call reduction: **50-200x**.

---

### Offline Mode

When network is unavailable or APIs time out frequently:

```bash
python scripts/dgra_cli_wrapper.py \
  --input-file patient.vcf.gz \
  --tissue general \
  --offline
```

Offline mode uses local cache (`references/offline_data/` gene JSONs). Genes not archived fallback to conservative rules.

---

### FAQ

| Issue | Cause | Solution |
|-------|-------|----------|
| `Invalid tissue 'xxx'` | Wrong tissue type | Use general / hematopoietic / cardiovascular / hepatic / renal / neurological / endocrine / metabolic / ophthalmic |
| `variants list is empty` | Empty input | Check file path and format |
| `Failed to write TSV` | Missing input fields | Ensure required fields CHROM/POS/REF/ALT/GENE exist |
| `dgra_core.py exited with code 1` | Core script failure | Check stderr output |
| `Offline mode: no cached data` | Offline mode but gene not archived | Run online once to build archive |
| Frequent API timeouts | Unstable network | Use `--offline` or check proxy settings |

---

### Testing

```bash
# Install dev dependencies
cd tests && pip install -r requirements-dev.txt

# Run all tests
pytest

# Run specific tier
pytest -m l2          # Unit tests
pytest -m l3          # Integration tests

# Coverage report
pytest --cov=scripts --cov-report=html
```

---

## Architecture

```
┌──────────────────────────────────────────────────────┐
│                    Output Layer                       │
│  Markdown (evidence chain) · JSON (structured)        │
│  Multi-organ joint report                             │
├──────────────────────────────────────────────────────┤
│                  Scoring Layer                        │
│  Tier 1/2/3 weighted scoring · ACMG evidence         │
│  SpliceAI modulation · NMD escape · Missense 5-layer │
│  ClinVar star confidence                              │
├──────────────────────────────────────────────────────┤
│                Annotation Layer                       │
│  8 API real-time queries · SQLite cache(30d)         │
│  5-layer offline resilience                           │
├──────────────────────────────────────────────────────┤
│                 Adapter Layer                         │
│  VEP / ANNOVAR / SnpEff → unified internal format    │
├──────────────────────────────────────────────────────┤
│               Input + QC + Filter                     │
│  VCF·Excel·TSV·free text · strict/clinical/broad    │
│  Raw VCF → VEP REST annotation · Disease-aware tx    │
│  Preflight health check · Two-phase pipeline          │
└──────────────────────────────────────────────────────┘
```

---

## Project Structure

```
dgra-genomic-risk/
├── scripts/                        # Core scripts (~30 modules)
│   ├── dgra_cli_wrapper.py         # ⭐ CLI entry point
│   ├── dgra_core.py                # Core engine (~200 lines, backward compat)
│   ├── dgra_api.py                 # 8 API wrappers + retry + cache
│   ├── dgra_cache.py               # SQLite cache management
│   ├── dgra_gnomad_local.py        # gnomAD local archive query
│   ├── build_gnomad_local.py       # gnomAD local archive builder
│   ├── gpa_pipeline.py             # Pipeline orchestration
│   ├── gpa_tier_classifier.py      # Tier 1/2/3 classification
│   ├── gpa_report.py               # Markdown/JSON report generation
│   ├── gpa_vcf_annotator.py        # Raw VCF → VEP REST annotation
│   ├── gpa_transcript_selector.py  # Disease-aware transcript selection
│   ├── gpa_phenotype_rescue.py     # Phenotype-driven VCF rescue search
│   ├── gpa_gene_set_builder.py     # Dynamic gene set builder (OMIM + HPO)
│   ├── gpa_phenotype_match.py      # LLM semantic phenotype matching
│   ├── gpa_phaser.py               # Phasing analysis
│   ├── gpa_multi_hit.py            # Multi-gene hit detection
│   ├── gpa_qc.py                   # QC checks
│   ├── gpa_i18n.py                 # Bilingual term mapping
│   ├── gpa_two_phase.py            # Two-phase pipeline optimization
│   ├── gpa_preflight.py            # Preflight health check
│   ├── gpa_workflow.py             # Workflow-as-Code definition
│   └── ...                         # Additional modules
├── references/                     # Reference data
│   ├── tissue_context.json         # 9 tissue profiles
│   ├── dgra.yaml                   # Runtime configuration
│   └── offline_data/               # Offline query archive
├── tests/                          # ~6,500 lines of tests
│   ├── e2e/                        # End-to-end tests
│   ├── test_l1_unit.py             # L1 unit tests
│   ├── test_l2_integration.py      # L2 integration tests
│   ├── test_l3_functional.py       # L3 functional tests
│   ├── test_l4_performance.py      # L4 performance tests
│   └── test_l5_edge_boundary.py    # L5 edge/boundary tests
├── docs/                           # Documentation
├── README.md                       # This file
├── CHANGELOG.md                    # Version history
├── config.json                     # Skill metadata
├── pyproject.toml                  # Python project config
└── requirements.txt                # Python dependencies
```

---

## Changelog Highlights

| Version | Date | Highlights |
|---------|------|------------|
| **v0.10.16** | 2026-06-10 | Phenotype Rescue workflow: dynamic gene set building + VCF rescue search for cases with no Tier 1 |
| **v0.10.15** | 2026-06-10 | VCF direct API · 9 tissue profiles · Report detail levels · gnomAD local archive · SQLite integrity fallback · Exact SO term matching |
| **v0.10.0** | 2026-05-25 | God Module split: dgra_core.py 2098 lines → 6 independent modules |
| **v0.9.0** | 2026-05-23 | Raw VCF end-to-end: VEP REST annotation + disease-aware transcript selection |
| **v0.8.0** | 2026-05-23 | SpliceAI splice prediction: Broad API + VEP REST fallback |
| **v0.7.2** | 2026-05-23 | ClinVar star confidence: practice_guideline 0.95 / single_submitter 0.40 |
| **v0.7.1** | 2026-05-23 | Pre-filtering + Chinese/English compatibility |
| **v0.6.0** | 2026-05-22 | Pseudogene architecture: 51 pairs + VAF pattern detection |
| **v0.5.x** | 2026-05-21 | Unified input layer + ACMG + NMD + weighted scoring |
| **v0.4.0** | 2026-05-19 | API-first architecture release |

See [CHANGELOG.md](CHANGELOG.md) for full details.

---

## Data Sources

Ensembl · UniProt · GTEx · gnomAD · ClinVar · HGNC · Orphanet · OMIM · ClinGen

---

## License

MIT-0

---

---

## Related Skills · 相关技能

| Skill · 技能 | Repo · 仓库 | Purpose · 用途 |
|---|---|---|
| **GPA Filter** | [lzr098/GPA-Filter](https://github.com/lzr098/GPA-Filter) | Genomic region pre-filter |
| **variant-impact** | [lzr098/variant-impact](https://github.com/lzr098/variant-impact) | Single variant ACMG classification |
| **disease-risk-query** | [lzr098/Disease-Risk-Query](https://github.com/lzr098/Disease-Risk-Query) | Disease-specific genetic risk |
| **sensory-genomics** | [lzr098/sensory-genomics-skill](https://github.com/lzr098/sensory-genomics-skill) | Five-sense genetic analysis |
| **PGS/PRS** | [lzr098](https://github.com/lzr098) | Polygenic risk scores from VCF |

---

**Maintainer**: [@lzr098](https://github.com/lzr098)  
**Current Version**: v0.10.15  
**Last Updated**: 2026-06-10
