# GPA 全面测试方案 (Comprehensive Test Plan)

> 版本: v0.10.0+
> 框架: pytest + pytest-asyncio + pytest-cov
> 更新日期: 2026-05-31

---

## 一、概述

本文档定义 GPA (Genomic Phenotype Association) 系统的全面测试策略，覆盖从单元测试到端到端测试的完整金字塔。测试方案采用 **L0~L6 分层 + E2E** 架构，配合 **录制-回放 (Record-Replay)** 机制处理外部 API 依赖。

### 1.1 设计原则

| 原则 | 说明 |
|------|------|
| **分层隔离** | L0~L5 按测试粒度和目的分层，避免重复和遗漏 |
| **录制-回放** | 外部 API 首次真实调用并录制，后续回放，兼顾真实性和稳定性 |
| **优先级驱动** | P0 测试必须全绿才能发布，P1 允许少量失败，P2 仅记录 |
| **模块全覆盖** | 所有 ~30 个模块均有对应测试文件，核心模块覆盖率 ≥ 80% |
| **向后兼容** | L6 回归测试确保版本间行为一致，防止修复引入新问题 |

### 1.2 测试范围

```
                    ┌─────────────────┐
                    │   E2E 端到端     │  ← 完整临床场景
                    │  4 个场景文件    │
                    └────────┬────────┘
                             │
              ┌──────────────┼──────────────┐
              │              │              │
       ┌──────▼──────┐ ┌────▼─────┐ ┌─────▼──────┐
       │ L6 回归测试  │ │ L5 边界  │ │ L4 性能    │
       │  向后兼容    │ │ 异常安全 │ │ 基准压力   │
       └─────────────┘ └──────────┘ └────────────┘
              │              │              │
       ┌──────▼──────┐ ┌────▼─────┐ ┌─────▼──────┐
       │ L3 集成测试  │ │ L2 单元  │ │ L1 静态    │
       │ 模块交互     │ │ 单模块   │ │ 导入结构   │
       └─────────────┘ └──────────┘ └────────────┘
              │
       ┌──────▼──────┐
       │ L0 契约测试  │
       │ 输入输出格式 │
       └─────────────┘
```

---

## 二、测试分层定义

### L0 — 契约测试 (Contract Tests)

验证公共 API 的输入/输出契约不因实现变更而破坏。

| 测试目标 | 验证内容 |
|---------|---------|
| 输入 Schema | Variant dict 必填字段、类型、取值范围 |
| 输出 Schema | JSON 报告字段完整性、Markdown 报告章节结构 |
| CLI 参数契约 | `dgra_cli_wrapper.py` 参数名称、类型、默认值 |
| 序列化契约 | `Variant` dataclass → JSON → `Variant` 往返无损 |

### L1 — 静态测试 (Static Tests)

不执行业务逻辑，验证代码结构和环境。

| 测试目标 | 验证内容 |
|---------|---------|
| 模块导入 | 所有 ~30 个模块可成功导入，无 `ImportError` |
| 数据结构 | `Variant`, `Evidence`, `GPAConfig` 默认值正确 |
| 循环依赖 | 模块间无循环导入导致的 `RecursionError` |
| 配置完整性 | `tissue_context.json` 格式正确，`FILTER_PRESETS` 存在 |
| i18n 映射 | `CHINESE_COLUMN_MAP` 非空，关键列存在 |

### L2 — 单元测试 (Unit Tests)

单个模块/函数的独立测试，**纯 mock**，不调用外部 API。

| 测试目标 | 验证内容 |
|---------|---------|
| 分级规则 | Tier 1/2/3 各触发条件正确 |
| 频率守卫 | EAS AF > 50% → Tier 3，Global AF > 80% → Tier 3 |
| ClinVar 处理 | Pathogenic 升级、Benign 降级、Conflicting 标记 |
| 星级加权 | `practice_guideline` weight=0.95，`single_submitter` weight=0.40 |
| NMD 预测 | escape/sensitive/possible_escape 判定 |
| SpliceAI | delta=0 降级、delta≥0.5 升级、阈值分类 |
| 转录本选择 | canonical/MANE 优先、歧义检测、LLM 辅助（mock） |
| 过滤预设 | strict/clinical/broad 各保留/排除规则 |
| 缓存 | set/get/TTL/并发 |
| 配置 | 默认值、环境变量覆盖、文件加载 |
| 输入解析 | VCF/TSV/CSV/Excel 解析、中文表头映射 |
| 适配器 | VEP/ANNOVAR/CSV/ClinVar 适配逻辑 |
| i18n | 中文 consequence → SO term、IMPACT 推断 |
| 批处理 | 分批逻辑、合并去重、最高 tier wins |
| 预检 | Python 依赖、本地工具、API 连通性、磁盘空间 |
| 相位 | GATK/distance/trio 相位计算 |
| 多基因命中 | cis/trans 检测、基因家族冗余 |
| QC | flag 设置、质量降级 |

### L3 — 集成测试 (Integration Tests)

验证模块间交互和完整 pipeline。

| 测试目标 | 验证内容 |
|---------|---------|
| Pipeline E2E | `run_dgra_pipeline()` 完整流程（offline） |
| CLI Wrapper | `run_gpa()` / `run_gpa_from_file()` 参数传递、输出格式 |
| VCF Pipeline | Raw VCF → VEP 注释 → 转录本选择 → 分级 → 报告 |
| 全量分析 | 多变异混合场景：TP53(T1) + DDX3X(T3) + BRCA1(T3) |
| 多器官联合 | `run_multi_organ_assessment()` 两个组织联合评估 |
| 两阶段 | Phase 1 triage → Phase 2 enrichment → 合并结果 |
| 预检集成 | `run_preflight_check()` → pipeline 入口接入 |

### L4 — 性能测试 (Performance Tests)

验证系统在大数据量下的表现。

| 测试目标 | 验证内容 | 阈值 |
|---------|---------|------|
| 100 变异 E2E | 100 变异完整 pipeline | < 60s |
| 1000 变异内存 | 1000 变异内存占用 | < 500MB |
| 缓存吞吐量 | SQLite 缓存 100 次查询 | > 2 req/s |
| 批量加速 | 批量查询 vs 单条查询 | 批量不慢于单条 2x |
| 批处理合并 | 3x1000 变异合并 | < 0.5s |

### L5 — 边界/异常测试 (Edge/Boundary Tests)

验证系统在极端和畸形输入下的鲁棒性。

| 测试目标 | 验证内容 |
|---------|---------|
| 空输入 | `[]` 空变异列表 → 报告生成，不崩溃 |
| 畸形 VCF | 缺失 CHROM/POS、POS="abc"、负坐标 |
| 极端 AF | AF=999.0、AF=-0.1、AF=1.0 |
| 极端 VAF | VAF=-0.1、VAF=1.5、VAF="NaN" |
| 编码问题 | Unicode 表型、中文基因型、UTF-8 BOM CSV |
| 缺失字段 | 缺少 IMPACT/Consequence/CLIN_SIG → 优雅降级 |
| 多等位基因 | ALT="A,C" → 正确处理 |
| 结构变异 | ALT="<DEL>" → 不崩溃 |
| chrMT 标准化 | "chrMT"/"MT" → "M" |
| 循环引用 | JSON 序列化不 `RecursionError` |
| DP=0 | 零深度 → 过滤或低质量标记 |
| 并发安全 | 多线程/异步并发操作缓存 |

### L6 — 回归测试 (Regression Tests)

验证新版本不破坏旧版本行为。

| 测试目标 | 验证内容 |
|---------|---------|
| v0.9 兼容 | v0.9.x 核心接口行为一致 |
| v0.9.5 热修 | aselect()、batch dedup、proxy config、i18n map |
| v0.10 重构 | God Module 拆分后公共 API 行为一致 |
| v0.10.1 预检 | Preflight 集成不改变分析结果 |

### E2E — 端到端测试

模拟真实用户场景的完整分析流程。

| 场景 | 输入 | 预期 |
|------|------|------|
| 产前遗传筛查 | 患者 VCF + general tissue | 致病/携带者变异报告 |
| 肿瘤体细胞驱动 | 肿瘤 VCF + somatic mode | 驱动突变 + 可干预性分级 |
| 药物基因组学 | 药物代谢基因子集 | CYP450/TPMT/DPYD 多态性报告 |
| 罕见病诊断 | 患者 VCF + 表型描述 + disease-aware | 候选致病变异排序 |

---

## 三、录制-回放机制 (Record-Replay)

### 3.1 设计目标

外部 API（Ensembl/UniProt/GTEx/gnomAD/NCBI 等）的测试需要真实响应，但：
- 真实 API 调用慢（>100ms/次）、不稳定（429/503）、有配额限制
- 纯 mock 无法验证 API 变更（如 gnomAD GraphQL schema 变更）

**录制-回放** 结合两者优势：
- **首次运行**：真实调用 API，响应保存到 `tests/recording/`
- **后续运行**：从录制文件加载响应，秒级完成
- **API 变更检测**：定期用 `--record-mode=refresh` 重新录制，diff 发现 schema 变更

### 3.2 实现机制

```
API Request
    │
    ▼
Generate cache key: "<api_name>_<variant_signature>_<endpoint_hash>"
    │
    ├── 录制文件存在 ──→ 加载 JSON 响应 ──→ 返回
    │
    └── 录制文件不存在
            │
            ├── record_mode=record ──→ 真实 API 调用 ──→ 保存 JSON ──→ 返回
            │
            └── record_mode=strict ──→ 抛出 MissingRecordingError
```

### 3.3 录制文件格式

```json
{
  "meta": {
    "api_name": "gnomad",
    "variant_id": "1-100000-A-G",
    "recorded_at": "2026-05-31T12:00:00Z",
    "api_version": "4.1",
    "recorded_by": "test_gnomad_api_rare_variant"
  },
  "request": {
    "url": "https://gnomad.broadinstitute.org/api",
    "query": "...",
    "headers": {}
  },
  "response": {
    "status": 200,
    "headers": {},
    "body": {...}
  }
}
```

### 3.4 使用方式

```bash
# 安装开发依赖
cd tests && pip install -r requirements-dev.txt

# 运行纯 mock 测试（默认，最快）
pytest -m "mock and not recording"

# 运行录制回放测试
pytest -m "recording"

# 重新录制所有 API 响应（需要网络）
pytest -m "recording" --record-mode=refresh

# 录制新响应（仅录制缺失的）
pytest -m "recording" --record-mode=record

# 严格模式：缺少录制则失败
pytest -m "recording" --record-mode=strict
```

### 3.5 pytest fixture (conftest.py)

```python
# conftest.py 中提供:
# - record_mode: 命令行 --record-mode 参数
# - api_recorder: 录制/回放上下文管理器
# - mock_gnomad / mock_ensembl / mock_gtex / mock_uniprot: 预录制 fixture

@pytest.fixture
def api_recorder(record_mode, tmp_path):
    """Provide a recording-replay context for API calls."""
    recorder = APIRecorder(
        recording_dir=Path(__file__).parent / "recording",
        mode=record_mode,
    )
    yield recorder
    recorder.save_index()
```

---

## 四、模块测试矩阵

### 4.1 核心模块（覆盖率目标 ≥ 80%）

#### gpa_tier_classifier.py (~650 行)

| 测试ID | 场景 | 输入 | 预期 | 优先级 |
|--------|------|------|------|--------|
| TIER-01 | ClinVar Pathogenic + HIGH + primary → T1 | RUNX1, Pathogenic, HIGH, gt=1/1 | Tier 1 | P0 |
| TIER-02 | ClinVar Pathogenic + HIGH + het → T2 | RUNX1, Pathogenic, HIGH, gt=0/1 | Tier 2 | P0 |
| TIER-03 | ClinVar Likely_pathogenic + HIGH → T2 | BRCA1, Likely_pathogenic, HIGH | Tier 2 | P0 |
| TIER-04 | ClinVar Benign → T3 | CFTR, Benign | Tier 3 | P0 |
| TIER-05 | ClinVar Conflicting → 不升级 + flag | BRCA1, Conflicting | qc_flag=CONFLICTING | P0 |
| TIER-06 | EAS AF > 50% → T3 | OR2B11, AF=0.55 EAS | Tier 3 | P0 |
| TIER-07 | Global AF > 80% → T3 | MAD2L2, AF=0.85 | Tier 3 | P0 |
| TIER-08 | 纯合 LoF + primary + rare → T1 | CFTR, 1/1, HIGH, NOT_CAPTURED | Tier 1 | P0 |
| TIER-09 | 杂合 LoF + primary → T2 | RUNX1, 0/1, HIGH, NOT_CAPTURED | Tier 2 | P0 |
| TIER-10 | gnomAD API_FAILED → 降级 T1→T2 | DDX3X, API_FAILED, 1/1, HIGH | Tier 2 | P0 |
| TIER-11 | NOT_CAPTURED + 纯合 → T1 | RUNX1, NOT_CAPTURED, 1/1, HIGH | Tier 1 | P0 |
| TIER-12 | phenotype_match < 0.6 → 降级 | BRCA1, Pathogenic, score=0.3 | Tier 2 | P1 |
| TIER-13 | UNKNOWN impact → 保守 HIGH | TEST, impact=UNKNOWN | 按 HIGH | P1 |
| TIER-14 | practice_guideline → weight=0.95 | BRCA1, practice_guideline | 高置信 | P1 |
| TIER-15 | single_submitter → weight=0.40 | BRCA1, single_submitter | 低置信 | P1 |
| TIER-16 | SpliceAI delta=0 → 降级 | splice_donor, delta=0 | Tier 降级 | P1 |
| TIER-17 | SpliceAI delta=0.8 → 升级 | splice_region, delta=0.8 | Tier 升级 | P1 |
| TIER-18 | NMD escape → PVS1 不适用 | frameshift, last exon | 不为 T1 | P1 |
| TIER-19 | NMD sensitive → PVS1 适用 | frameshift, internal exon | Tier 1 | P1 |
| TIER-20 | NMD possible_escape → 降级 | frameshift, penultimate | Tier 2 | P1 |
| TIER-21 | 假基因干扰 → confidence=LOW | VWF, VAF=0.15 | tier 不变 | P0 |
| TIER-22 | 转录本歧义 (NR_/XM_) → 降级 | NR_001, HIGH | 不为 T1 | P1 |
| TIER-23 | X-linked 女性杂合 → 下调 | X, gt=0/1, haplosufficient | Tier 下调 | P1 |
| TIER-24 | 基因家族完全代偿 → 降级 | HLA-A, complete | T1→T2 | P1 |
| TIER-25 | Somatic VAF>0.5 → T3 | TP53, VAF=0.98, somatic | Tier 3 | P1 |
| TIER-26 | Somatic TSG LOF → T1 | TP53, HIGH, is_tsg, somatic | Tier 1 | P1 |
| TIER-27 | 药物代谢基因 → T2 | CYP2D6, MODERATE | Tier 2 | P1 |
| TIER-28 | 无组织相关性 fast-track → T3 | BRCA1, no relevance | Tier 3 | P0 |
| TIER-29 | 证据链完整性 | 任意 variant | evidence_chain 非空 | P0 |
| TIER-30 | qc_flags 设置 | Conflicting/Benign/common | 对应 flag | P0 |

#### gpa_pipeline.py (~500 行)

| 测试ID | 场景 | 预期 | 优先级 |
|--------|------|------|--------|
| PIPE-01 | 空变异列表 | report_markdown + json_report 存在 | P0 |
| PIPE-02 | 单变异完整报告 | Markdown + JSON 完整 | P0 |
| PIPE-03 | 3 变异混合 | TP53=T1, DDX3X=T3, BRCA1=T3 | P0 |
| PIPE-04 | offline 模式 | 无 API 调用，仅用缓存+本地 | P0 |
| PIPE-05 | somatic 模式 | TSG 逻辑生效，VAF>0.5 过滤 | P1 |
| PIPE-06 | phenotype 参数 | phenotype_match_score 计算 | P1 |
| PIPE-07 | 预过滤 strict | 仅 HIGH/MODERATE 进入分级 | P1 |
| PIPE-08 | 预过滤 clinical | splice_region LOW 保留 | P1 |
| PIPE-09 | 多器官联合 | joint_risk_matrix + 联合报告 | P1 |
| PIPE-10 | 100 变异批量 | 自动分批，合并结果 | P1 |

#### gpa_report.py (~450 行)

| 测试ID | 场景 | 预期 | 优先级 |
|--------|------|------|--------|
| RPT-01 | Markdown Tier 章节 | 含 T1/T2/T3 变异详情 | P0 |
| RPT-02 | Markdown 方法学附录 | 含分级规则说明 | P0 |
| RPT-03 | Markdown 转录本评估 | 歧义时含 alternatives | P1 |
| RPT-04 | Markdown QC 标记 | 含所有 qc_flags | P0 |
| RPT-05 | Markdown 多基因命中 | 含 cis/trans 判断 | P1 |
| RPT-06 | JSON 字段完整性 | meta/summary/variants/qc | P0 |
| RPT-07 | JSON summary 计数 | tier1_gene_count 等正确 | P0 |
| RPT-08 | 空报告 | 0 变异时结构完整 | P0 |
| RPT-09 | 中文输出 | 中文表头/术语正确 | P1 |

#### dgra_api.py (~836 行)

| 测试ID | 场景 | API | 预期 | 优先级 |
|--------|------|-----|------|--------|
| API-01 | Ensembl VEP REST 单条 | Ensembl | 200 + 注释结果 | P0 |
| API-02 | Ensembl VEP REST 批量 | Ensembl | 200 + 批量结果 | P0 |
| API-03 | UniProt 蛋白域 | UniProt | 200 + domain info | P1 |
| API-04 | GTEx 表达数据 | GTEx | 200 + TPM | P1 |
| API-05 | gnomAD 罕见变异 | gnomAD | 200 + AF rare | P0 |
| API-06 | gnomAD 常见变异 | gnomAD | 200 + AF common | P0 |
| API-07 | gnomAD API_FAILED | gnomAD | status=API_FAILED | P0 |
| API-08 | NCBI ClinVar | NCBI | 200 + ClinVar 数据 | P1 |
| API-09 | HGNC 基因符号 | HGNC | 200 + 标准化符号 | P1 |
| API-10 | MyVariant.info | MyVariant | 200 + 补充注释 | P1 |
| API-11 | HTTP 429 重试 | Ensembl | Retry-After 等待 | P1 |
| API-12 | HTTP 503 退避 | Ensembl | 指数退避 1→2→4s | P1 |
| API-13 | 超时重试 | Ensembl | 3 次重试 | P1 |
| API-14 | 缓存命中 | 任意 | 不触发真实请求 | P0 |
| API-15 | 速率限制 | 任意 | 不超过限制 | P1 |

#### gpa_vcf_annotator.py (~398 行)

| 测试ID | 场景 | 预期 | 优先级 |
|--------|------|------|--------|
| VCF-01 | Raw VCF 检测（无 CSQ） | 识别为 RAW_VCF | P0 |
| VCF-02 | Annotated VCF 检测（有 CSQ） | 识别为 ANNOTATED_VCF | P0 |
| VCF-03 | VEP REST 注释 | 生成注释字段 | P0 |
| VCF-04 | 批量注释（>500） | 自动分批并发 | P1 |
| VCF-05 | disease-aware 选择 | 触发转录本选择器 | P1 |
| VCF-06 | 歧义检测 | is_ambiguous=True | P1 |
| VCF-07 | LLM 辅助选择（mock） | method=llm | P2 |
| VCF-08 | 大数据集 >5000 | 强制走 VEP API | P1 |
| VCF-09 | 代理配置 | trust_env 正确 | P1 |

#### dgra_cli_wrapper.py (~311 行)

| 测试ID | 场景 | 预期 | 优先级 |
|--------|------|------|--------|
| CLI-01 | 空 variants 列表 | success=False, "empty" | P0 |
| CLI-02 | 无效 tissue | success=False, "Invalid tissue" | P0 |
| CLI-03 | 单变异成功 | success=True, report_md 存在 | P0 |
| CLI-04 | 从文件加载 VCF | 正确解析并分析 | P0 |
| CLI-05 | 从文件加载 TSV | 正确解析并分析 | P0 |
| CLI-06 | 自动分批（>500） | 分多批处理 | P1 |
| CLI-07 | somatic 标志 | somatic_mode 传递 | P1 |
| CLI-08 | spliceai 标志 | spliceai_enabled 传递 | P1 |
| CLI-09 | filter-preset | 预过滤生效 | P1 |
| CLI-10 | 输出保存 | --output 文件写入 | P0 |
| CLI-11 | multi_organ 参数 | 多器官评估 | P1 |
| CLI-12 | database_version | 传递正确 | P1 |

---

### 4.2 其他模块（覆盖率目标 ≥ 60%）

#### gpa_preflight.py (~16641 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| PFL-01 | Python 依赖齐全 | all required = True |
| PFL-02 | Python 依赖缺失 | blocker 包含缺失包 |
| PFL-03 | 本地文件齐全 | tissue_context.json 存在 |
| PFL-04 | API 全部连通 | all api = reachable |
| PFL-05 | API 部分不通 | warning 包含不通 API |
| PFL-06 | 磁盘空间充足 | disk_space = True |
| PFL-07 | 磁盘空间不足 | disk_space = False |
| PFL-08 | 建议 continue | 无 blocker |
| PFL-09 | 建议 offline | API 不通但其他 OK |
| PFL-10 | 建议 abort | Python 依赖缺失 |

#### gpa_two_phase.py (~4784 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| TPH-01 | Phase 1 triage | 本地快速筛选候选变异 |
| TPH-02 | Phase 2 enrichment | API 补充注释候选变异 |
| TPH-03 | 阈值过滤 | AF/IMPACT 阈值正确应用 |
| TPH-04 | 结果合并 | Phase1 + Phase2 结果合并 |
| TPH-05 | 无候选变异 | Phase2 不执行 |

#### gpa_workflow*.py (3 个文件)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| WFL-01 | 工作流编排 | 步骤按序执行 |
| WFL-02 | PM 决策 | 分析目的识别正确 |
| WFL-03 | runner 执行 | 异步任务调度 |

#### gpa_transcript_selector.py (~16641 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| TXS-01 | 空列表 | method=none |
| TXS-02 | 单转录本 | method=single |
| TXS-03 | canonical 优先 | canonical+HIGH 选中 |
| TXS-04 | MANE Select 优先 | MANE 评分更高 |
| TXS-05 | 歧义检测 | 差距≤10分 → is_ambiguous |
| TXS-06 | LLM 辅助（mock） | method=llm |
| TXS-07 | 异步上下文 | 不 RuntimeError |

#### dgra_input_parsers.py (~26113 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| PAR-01 | VCF 解析 | CHROM/POS/REF/ALT 提取 |
| PAR-02 | VCF 无 CSQ | 标记为 RAW_VCF |
| PAR-03 | VCF 有 CSQ | 解析 CSQ 字段 |
| PAR-04 | TSV 解析 | 制表符分隔 |
| PAR-05 | CSV 解析 | 逗号分隔 |
| PAR-06 | Excel 解析 | .xlsx 读取 |
| PAR-07 | 中文表头 | translate_chinese_headers |
| PAR-08 | 字段映射 | 标准字段名映射 |

#### dgra_adapters.py (~18830 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| ADP-01 | VEPAdapter | VEP 字段 → 标准字段 |
| ADP-02 | ANNOVARAdapter | ANNOVAR 字段 → 标准字段 |
| ADP-03 | CSVAdapter | 通用 CSV 字段映射 |
| ADP-04 | ClinVarAdapter | ClinVar VCF 解析 |
| ADP-05 | 中文 exonic func | 移码变异 → frameshift_variant |

#### gpa_i18n.py (~13977 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| I18-01 | 中文表头 13 列 | 全部正确映射 |
| I18-02 | 错义变异 | → missense_variant |
| I18-03 | 无义变异 | → stop_gained |
| I18-04 | IMPACT 推断 | 高 → HIGH |
| I18-05 | ClinVar 中文 | 致病 → Pathogenic |

#### dgra_cache.py (~13904 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| CCH-01 | set/get | 数据存取正确 |
| CCH-02 | TTL 过期 | 过期后返回 None |
| CCH-03 | 并发安全 | 多线程不冲突 |
| CCH-04 | JSON 序列化 | dict/list 正确序列化 |
| CCH-05 | 大量数据 | 性能不下降 |

#### dgra_config.py (~12713 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| CFG-01 | GPAConfig 默认值 | min_dp=20, min_gq=90.0 |
| CFG-02 | 环境变量覆盖 | HTTPS_PROXY 生效 |
| CFG-03 | 文件加载 | dgra.yaml 解析 |
| CFG-04 | DGRAGlobalConfig | proxy=None 默认 |
| CFG-05 | __DIRECT__ 处理 | 跳过代理 |

#### dgra_variant_filter.py (~9981 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| FLT-01 | strict 预设 | 仅 HIGH/MODERATE |
| FLT-02 | clinical 预设 | splice_region LOW 保留 |
| FLT-03 | broad 预设 | HIGH/MODERATE/LOW 全保留 |
| FLT-04 | 统计信息 | input/output/excluded 计数 |

#### dgra_batch_runner.py (~17498 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| BAT-01 | 分批逻辑 | >500 自动分多批 |
| BAT-02 | 去重 | 同变异只出现一次 |
| BAT-03 | 最高 tier wins | T2+T1 → 保留 T1 |
| BAT-04 | 失败批次跳过 | 失败批次不影响整体 |
| BAT-05 | 性能 | 3x1000 合并 <0.5s |

#### dgra_splice_predictor.py (~21600 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| SPL-01 | should_query | splice_acceptor → True |
| SPL-02 | should_query | missense → False |
| SPL-03 | delta=0 → none | canonical threshold |
| SPL-04 | delta=0.55 → strong | canonical threshold |
| SPL-05 | delta=0.25 → strong | splice_region threshold |
| SPL-06 | API 查询（mock） | 返回 SpliceAIResult |
| SPL-07 | 缓存 | 重复查询命中缓存 |

#### gpa_phenotype_match.py (~9247 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| PHM-01 | 顿号分隔 | 3 个表型 |
| PHM-02 | 逗号分隔 | 3 个表型 |
| PHM-03 | 句号分隔 | 3 个表型 |
| PHM-04 | 混合分隔 | 4 个表型 |
| PHM-05 | LLM 匹配（mock） | 返回 score + explanation |

#### dgra_myvariant.py (~19852 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| MYV-01 | 查询 | 返回注释数据 |
| MYV-02 | 解析 | AF/CADD/ClinVar 提取 |
| MYV-03 | 错误 | API 失败 graceful |

#### gpa_phaser.py (~8282 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| PHA-01 | GATK 相位 | phased/unphased 判断 |
| PHA-02 | 距离相位 | <50bp → likely_phased |
| PHA-03 | trio 相位 | 父母子代传递 |
| PHA-04 | 边界修正 | <50bp 不误判 unphased |

#### gpa_multi_hit.py (~6252 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| MHT-01 | 多基因命中 | 检测同一基因多个变异 |
| MHT-02 | cis 判断 | 同染色体 → cis |
| MHT-03 | trans 判断 | 不同染色体 → trans |
| MHT-04 | 基因家族冗余 | complete compensation |

#### gpa_qc.py (~4784 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| QC-01 | 质量检查 | DP<GQ 阈值 → flag |
| QC-02 | flag 添加 | qc_flags 追加 |
| QC-03 | 质量降级 | low confidence 标记 |

#### dgra_gene_sync.py (~22341 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| GSY-01 | 基因列表获取 | 返回 tissue_genes |
| GSY-02 | 合并 | 多源合并去重 |

#### dgra_pseudogene_sync.py (~21600 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| PSY-01 | 假基因列表 | 返回 pseudogene set |
| PSY-02 | 同步 | GTF 下载+解析 |

#### dgra_build_state.py (~4331 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| BLD-01 | 保存/加载 | JSON 持久化 |
| BLD-02 | 步骤状态 | pending/in_progress/complete/failed |
| BLD-03 | 上下文管理器 | BuildStep 自动状态 |
| BLD-04 | 恢复 | 删除后重新执行 |

#### convert_csv_to_gpa.py (~6404 行)

| 测试ID | 场景 | 预期 |
|--------|------|------|
| CSV-01 | CSV 转 GPA | 字段映射正确 |
| CSV-02 | 中文表头 | 自动翻译 |

---

## 五、覆盖率目标

### 5.1 核心模块（≥ 80%）

| 模块 | 行数 | 目标 | 当前状态 |
|------|------|------|----------|
| gpa_tier_classifier.py | ~650 | ≥ 80% | 待测 |
| gpa_pipeline.py | ~500 | ≥ 80% | 待测 |
| gpa_report.py | ~450 | ≥ 80% | 待测 |
| dgra_api.py | ~836 | ≥ 80% | 待测 |
| gpa_vcf_annotator.py | ~398 | ≥ 80% | 待测 |
| dgra_cli_wrapper.py | ~311 | ≥ 80% | 待测 |

### 5.2 其他模块（≥ 60%）

| 模块 | 行数 | 目标 | 当前状态 |
|------|------|------|----------|
| gpa_preflight.py | ~400 | ≥ 60% | 待测 |
| gpa_two_phase.py | ~350 | ≥ 60% | 待测 |
| gpa_workflow*.py | ~600 | ≥ 60% | 待测 |
| gpa_transcript_selector.py | ~350 | ≥ 60% | 待测 |
| dgra_input_parsers.py | ~500 | ≥ 60% | 待测 |
| dgra_adapters.py | ~400 | ≥ 60% | 待测 |
| gpa_i18n.py | ~300 | ≥ 60% | 待测 |
| dgra_cache.py | ~300 | ≥ 60% | 待测 |
| dgra_config.py | ~250 | ≥ 60% | 待测 |
| dgra_variant_filter.py | ~200 | ≥ 60% | 待测 |
| dgra_batch_runner.py | ~300 | ≥ 60% | 待测 |
| dgra_splice_predictor.py | ~350 | ≥ 60% | 待测 |
| gpa_phenotype_match.py | ~200 | ≥ 60% | 待测 |
| dgra_myvariant.py | ~300 | ≥ 60% | 待测 |
| gpa_phaser.py | ~200 | ≥ 60% | 待测 |
| gpa_multi_hit.py | ~150 | ≥ 60% | 待测 |
| gpa_qc.py | ~100 | ≥ 60% | 待测 |
| dgra_gene_sync.py | ~350 | ≥ 60% | 待测 |
| dgra_pseudogene_sync.py | ~300 | ≥ 60% | 待测 |
| dgra_build_state.py | ~150 | ≥ 60% | 待测 |
| convert_csv_to_gpa.py | ~150 | ≥ 60% | 待测 |

---

## 六、运行方式

### 6.1 安装

```bash
cd /Users/zhaorongli/.workbuddy/skills/dgra-genomic-risk
pip install -r requirements.txt
pip install -r tests/requirements-dev.txt
```

### 6.2 快速运行

```bash
# 运行所有测试
cd tests && pytest

# 运行指定分层
pytest -m l2          # 单元测试
pytest -m l3          # 集成测试
pytest -m "l2 or l3"  # 单元+集成

# 运行指定优先级
pytest -m p0          # 仅 P0
pytest -m "p0 or p1"  # P0+P1

# 运行纯 mock（无网络）
pytest -m "mock and not recording"

# 运行录制回放
pytest -m recording

# 运行性能测试
pytest -m l4

# 运行回归测试
pytest -m l6

# 运行端到端
pytest -m e2e

# 运行指定模块
pytest l2_unit/test_tier_classifier.py

# 运行指定组织上下文
pytest -m hematopoietic
```

### 6.3 覆盖率报告

```bash
# 生成覆盖率报告
pytest --cov=scripts --cov-report=html

# 查看 HTML 报告
open tests/htmlcov/index.html

# 仅查看未覆盖行
pytest --cov=scripts --cov-report=term-missing
```

### 6.4 录制-回放模式

```bash
# 录制缺失的 API 响应（需要网络）
pytest -m recording --record-mode=record

# 重新录制所有（刷新已有录制）
pytest -m recording --record-mode=refresh

# 严格模式：缺少录制则失败
pytest -m recording --record-mode=strict

# 回放模式（默认，不需要网络）
pytest -m recording --record-mode=playback
```

---

## 七、文件清单

### 7.1 新增文件

| 路径 | 说明 |
|------|------|
| `tests/pytest.ini` | pytest 配置 |
| `tests/requirements-dev.txt` | 开发依赖 |
| `tests/TEST_PLAN.md` | 本文档 |
| `tests/recording/` | API 录制响应目录 |
| `tests/recording/gnomad/` | gnomAD 录制 |
| `tests/recording/ensembl/` | Ensembl 录制 |
| `tests/recording/gtex/` | GTEx 录制 |
| `tests/recording/uniprot/` | UniProt 录制 |
| `tests/recording/clinvar/` | ClinVar 录制 |
| `tests/recording/myvariant/` | MyVariant 录制 |
| `tests/recording/hgnc/` | HGNC 录制 |
| `tests/recording/ncbi/` | NCBI 录制 |
| `tests/l0_contract/test_input_schema.py` | 输入契约测试 |
| `tests/l0_contract/test_output_schema.py` | 输出契约测试 |
| `tests/l1_static/test_imports.py` | 导入测试 |
| `tests/l1_static/test_dataclasses.py` | 数据结构测试 |
| `tests/l1_static/test_circular_deps.py` | 循环依赖测试 |
| `tests/l2_unit/test_tier_classifier.py` | 分级器单元测试 |
| `tests/l2_unit/test_pipeline.py` | Pipeline 单元测试 |
| `tests/l2_unit/test_report.py` | 报告单元测试 |
| `tests/l2_unit/test_api.py` | API 单元测试 |
| `tests/l2_unit/test_vcf_annotator.py` | VCF 注释单元测试 |
| `tests/l2_unit/test_transcript_selector.py` | 转录本选择测试 |
| `tests/l2_unit/test_phaser.py` | 相位分析测试 |
| `tests/l2_unit/test_multi_hit.py` | 多基因命中测试 |
| `tests/l2_unit/test_qc.py` | QC 测试 |
| `tests/l2_unit/test_input_parsers.py` | 输入解析测试 |
| `tests/l2_unit/test_adapters.py` | 适配器测试 |
| `tests/l2_unit/test_i18n.py` | 国际化测试 |
| `tests/l2_unit/test_cache.py` | 缓存测试 |
| `tests/l2_unit/test_config.py` | 配置测试 |
| `tests/l2_unit/test_variant_filter.py` | 过滤测试 |
| `tests/l2_unit/test_batch_runner.py` | 批处理测试 |
| `tests/l2_unit/test_splice_predictor.py` | SpliceAI 测试 |
| `tests/l2_unit/test_phenotype_match.py` | 表型匹配测试 |
| `tests/l2_unit/test_preflight.py` | 预检测试 |
| `tests/l2_unit/test_two_phase.py` | 两阶段测试 |
| `tests/l2_unit/test_workflow.py` | 工作流测试 |
| `tests/l2_unit/test_pseudogene_sync.py` | 假基因同步测试 |
| `tests/l2_unit/test_build_state.py` | 构建状态测试 |
| `tests/l3_integration/test_pipeline_e2e.py` | Pipeline 集成 |
| `tests/l3_integration/test_cli_wrapper.py` | CLI 集成 |
| `tests/l3_integration/test_full_analysis.py` | 全量分析 |
| `tests/l3_integration/test_vcf_pipeline.py` | VCF Pipeline |
| `tests/l4_performance/test_benchmark.py` | 基准测试 |
| `tests/l4_performance/test_memory.py` | 内存测试 |
| `tests/l4_performance/test_cache_throughput.py` | 缓存吞吐 |
| `tests/l5_edge/test_empty_input.py` | 空输入边界 |
| `tests/l5_edge/test_malformed_data.py` | 畸形数据 |
| `tests/l5_edge/test_encoding.py` | 编码问题 |
| `tests/l5_edge/test_concurrency.py` | 并发安全 |
| `tests/l5_edge/test_extreme_values.py` | 极端值 |
| `tests/l6_regression/test_v09_compat.py` | v0.9 兼容 |
| `tests/l6_regression/test_v095_hotfix.py` | v0.9.5 回归 |
| `tests/l6_regression/test_v010_regression.py` | v0.10 回归 |
| `tests/e2e/test_clinical_scenarios.py` | 临床场景 |
| `tests/e2e/test_germline_risk.py` | 胚系风险 |
| `tests/e2e/test_somatic_driver.py` | 体细胞驱动 |
| `tests/e2e/test_pharmacogenomics.py` | 药物基因组 |

### 7.2 修改文件

| 路径 | 修改内容 |
|------|----------|
| `tests/conftest.py` | 重构为 pytest fixtures，添加录制-回放基础设施 |
| `SKILL.md` | 添加"测试方案"章节引用 |
| `requirements.txt` | 可选：添加测试依赖（或独立 requirements-dev.txt） |

---

## 八、里程碑

| 阶段 | 内容 | 预计用例数 |
|------|------|-----------|
| Phase 1 | L0+L1 静态测试 + 基础设施搭建 | ~20 |
| Phase 2 | L2 核心模块单元测试（tier/pipeline/report/api/vcf/cli） | ~80 |
| Phase 3 | L2 其他模块单元测试 | ~60 |
| Phase 4 | L3 集成测试 + 录制-回放 | ~20 |
| Phase 5 | L4 性能 + L5 边界 + L6 回归 | ~30 |
| Phase 6 | E2E 端到端 | ~16 |
| **合计** | | **~226** |

---

## 九、维护指南

### 新增模块时

1. 在 `l2_unit/` 下创建 `test_<module>.py`
2. 在本文档"模块测试矩阵"中添加对应测试用例
3. 如模块涉及外部 API，在 `conftest.py` 中添加 mock fixture

### 新增 API 时

1. 在 `tests/recording/` 下创建对应子目录
2. 在 `conftest.py` 的 `APIRecorder` 中注册 API 名称和 base URL
3. 首次运行 `--record-mode=record` 生成录制文件

### 版本发布前

```bash
# 1. 运行全部 P0 测试
pytest -m p0

# 2. 运行回归测试
pytest -m l6

# 3. 检查核心模块覆盖率
pytest --cov=scripts --cov-report=term-missing

# 4. 确认核心模块 ≥ 80%
# 5. 全部通过后方可发布
```
