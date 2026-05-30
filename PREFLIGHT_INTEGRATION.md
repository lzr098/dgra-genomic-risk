# GPA Preflight 集成方案

## 1. 设计目标

每次接到全新分析任务时，先执行一次**前置可用性检查**，确认所有依赖就绪后再启动耗时较长的分析流程。避免分析到一半才发现 API 不通、工具缺失或磁盘空间不足。

## 2. 模块文件

- **新建**: `scripts/gpa_preflight.py` — 预检核心模块
- **修改**: `scripts/gpa_pipeline.py` — 在 pipeline 入口接入
- **修改**: `scripts/dgra_core.py` — 在 analyze_variants 入口接入（可选）
- **修改**: `scripts/gpa_two_phase.py` — 在两阶段分析入口接入（可选）

## 3. 检查范围（六大类）

| 分类 | 检查项 | 必须/可选 | 说明 |
|------|--------|-----------|------|
| **python_deps** | aiohttp | 必须 | 所有 API 调用的基础 |
| | cyvcf2 | 必须 | VCF/BCF 解析 |
| | yaml (PyYAML) | 必须 | 配置解析 |
| **local_tools** | vep | 可选 | 仅 `--annotator vep_local` 时需要 |
| | git | 可选 | 仅版本信息收集 |
| **api_connectivity** | Ensembl VEP REST | 必须 | 核心注释 API |
| | UniProt REST | 必须 | 蛋白域/功能注释 |
| | GTEx Portal | 必须 | 组织表达数据 |
| | gnomAD API | 必须 | 人群频率 |
| | NCBI EUtils | 必须 | ClinVar / 文献 |
| | HGNC | 必须 | 基因符号标准化 |
| | MyVariant.info | 必须 | 补充注释 (AF/CADD/ClinVar) |
| | SpliceAI Broad | 可选 | 有 VEP REST fallback |
| **local_files** | references/ 目录 | 必须 | 配置文件存放地 |
| | tissue_context.json | 必须 | Tier 分类必需 |
| | cache/ 目录 | 必须 | 自动创建 |
| | dgra_cache.db | 可选 | 首次运行自动创建 |
| | dgra.yaml | 可选 | 自定义配置 |
| **disk_space** | 磁盘可用空间 | 必须 | 默认阈值 500 MB |
| **network_env** | 直连 NCBI | 可选 | 网络环境探测 |
| | 代理可用性 | 可选 | 自动探测常见代理端口 |

## 4. 决策流程

```
接到新任务
    │
    ▼
调用 run_preflight_check(config)
    │
    ├── Python 包缺失 ──→ abort（必须 pip install）
    ├── 本地文件缺失 ──→ abort（必须修复目录结构）
    ├── 磁盘空间不足 ──→ abort（必须清理/扩容）
    │
    ├── API 不通 ──→ 建议 offline 模式 或 配置代理
    ├── 可选工具缺失 ──→ 继续（运行时自动跳过）
    │
    ▼
全部通过 ──→ 直接进入分析
有 blocker ──→ 输出交互提示，等待用户选择:
              [1] 离线模式 — 跳过 API，仅用缓存+本地
              [2] 忽略警告继续 — 缺失功能自动跳过
              [3] 中止 — 修复环境后重试
```

## 5. 核心 API

```python
from gpa_preflight import run_preflight_check, interactive_prompt, suggest_action

# 执行检查
report = await run_preflight_check(config)

# 查询结果
if report.is_ready():
    # 全部就绪
    pass
else:
    blockers = report.blockers()   # 必须修复的项
    warnings = report.warnings()   # 可选缺失的项

# 获取建议动作（自动决策，无用户交互时）
action = suggest_action(report)   # "continue" | "offline" | "abort"

# 生成交互提示文本（有用户交互时）
prompt_text = interactive_prompt(report)
# 返回文本中包含 [1]/[2]/[3] 选项
```

## 6. 集成到 gpa_pipeline.py

在 `run_gpa_pipeline()` 函数开头添加：

```python
from gpa_preflight import run_preflight_check, suggest_action

async def run_gpa_pipeline(variants, config, client, tissue_profile="general"):
    # v0.10.1: Preflight health check
    from gpa_preflight import run_preflight_check, suggest_action
    preflight = await run_preflight_check(config)
    if not preflight.is_ready():
        action = suggest_action(preflight)
        if action == "abort":
            print("[GPA Preflight] 环境检查未通过，中止分析。")
            print(preflight.to_markdown())
            return {"error": "Preflight failed", "report": preflight.to_dict()}
        elif action == "offline":
            print("[GPA Preflight] 切换到离线模式（跳过所有 API 调用）")
            config.offline_mode = True
        else:
            print("[GPA Preflight] 忽略警告，继续分析（部分功能可能不可用）")
    # ... 原有 pipeline 逻辑
```

## 7. 集成到 dgra_core.py

在 `analyze_variants()` 函数开头添加：

```python
async def analyze_variants(...) -> Dict[str, Any]:
    # v0.10.1: Preflight check
    from gpa_preflight import run_preflight_check, suggest_action
    preflight = await run_preflight_check(config)
    if not preflight.is_ready():
        action = suggest_action(preflight)
        if action == "abort":
            return {
                "tier_report": "# 分析中止\n\n前置环境检查未通过。",
                "summary": {"error": "Preflight failed", "report": preflight.to_dict()},
            }
        elif action == "offline":
            config.offline_mode = True
    # ... 原有分析逻辑
```

## 8. 集成到 gpa_two_phase.py

在 `run_two_phase_analysis()` 函数开头添加：

```python
async def run_two_phase_analysis(...):
    from gpa_preflight import run_preflight_check, suggest_action
    preflight = await run_preflight_check(config)
    if not preflight.is_ready():
        action = suggest_action(preflight)
        if action == "abort":
            return {"error": "Preflight failed", "report": preflight.to_dict()}
        elif action == "offline":
            config.offline_mode = True
    # ... 原有 Phase 1 / Phase 2 逻辑
```

## 9. 离线模式行为

当 `config.offline_mode = True` 时：
- `run_preflight_check()` 自动跳过所有 `api_connectivity` 检查
- Pipeline 中的 API 调用被缓存命中或返回 None
- Tier 分类仍可进行（依赖本地文件 `tissue_context.json`）
- 报告会标注 "离线模式"，提示部分证据缺失

## 10. 扩展检查项

如需新增检查（例如 `bcftools`、`tabix`、`pandoc`），只需修改 `gpa_preflight.py`：

```python
def check_local_tools() -> List[CheckItem]:
    return [
        _check_cli_tool("vep", required=False),
        _check_cli_tool("git", required=False),
        _check_cli_tool("bcftools", required=False),   # ← 新增
        _check_cli_tool("tabix", required=False),      # ← 新增
    ]
```

如需新增 API 检查，添加到 `_API_CHECKS` 字典即可：

```python
_API_CHECKS = {
    # ... 现有检查
    "orphanet": (
        "https://api.orphacode.org/nomenclature/orphanumber/166001/genes",
        5,
        None,
    ),
}
```

## 11. 命令行自测

```bash
cd scripts
python gpa_preflight.py
```

输出 Markdown 格式的检查报告，包括延迟、建议和总体结论。
