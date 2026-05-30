# GPA 审查标准自动化遵守机制

**版本**: v1.0  
**生效日期**: 2026-05-30  

---

## 机制概览

本机制确保每次修改 GPA skill 代码时，**无法绕过**审查标准。采用"三道防线"设计：

```
┌─────────────────────────────────────────────────────────────────────┐
│                        三道防线模型                                  │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  第一道: 开发者自检                                                   │
│  ┌─────────────────┐    运行: python3 gpa_review_gate.py            │
│  │ 修改代码        │ ──▶ 发现 Blocker → 强制修复后才能继续           │
│  └─────────────────┘    发现 Critical → 建议修复                      │
│                                                                      │
│  第二道: 专家审查 (AI)                                                │
│  ┌─────────────────┐    @gpa-code-reviewer 专家                       │
│  │ 提交前询问      │ ──▶ "审查这些修改"                                │
│  └─────────────────┘    专家运行完整 AST 分析 + 逻辑审查               │
│                                                                      │
│  第三道: 自动化门禁 (CI)                                              │
│  ┌─────────────────┐    GitHub Actions / 本地 git hook               │
│  │ git commit/push │ ──▶ 自动运行 gpa_review_gate.py                 │
│  └─────────────────┘    Blocker 存在则阻止提交                        │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 防线一: 开发者自检脚本

### 脚本位置
```
~/.workbuddy/skills/dgra-genomic-risk/scripts/gpa_review_gate.py
```

### 用法

```bash
# 检查 git 工作区中修改的文件
python3 gpa_review_gate.py

# 检查所有文件
python3 gpa_review_gate.py --all

# 检查指定文件
python3 gpa_review_gate.py dgra_core.py gpa_pipeline.py

# 查看帮助
python3 gpa_review_gate.py --help
```

### 退出码

| 退出码 | 含义 | 行动 |
|--------|------|------|
| 0 | 无问题 | 可提交 |
| 1 | 有 Blocker | **禁止提交，必须修复** |
| 2 | 有 Critical | 建议修复后可提交 |
| 3 | 仅 Nit | 可选修复 |

### 自检流程（每次修改代码时执行）

```bash
# 1. 修改代码...

# 2. 运行门禁
python3 gpa_review_gate.py

# 3. 如果返回 1，按报告逐项修复
#    如果返回 0/2/3，可以继续

# 4. 提交
```

---

## 防线二: AI 专家审查

### 激活专家

在 WorkBuddy 中直接使用：

```
@gpa-code-reviewer 审查 dgra_core.py 的修改
```

或首次对话：
```
审查这个 GPA skill 的代码修改，找出所有 blocker 级别的 bug
```

### 专家能力

- 运行自动化脚本 (`gpa_review_gate.py`)
- AST 深层分析（循环导入、未定义变量、类型不一致）
- 逻辑审查（业务错误、边界条件）
- 生成修复补丁

### 何时调用专家

| 场景 | 建议 |
|------|------|
| 修改超过 3 个文件 | 调用专家做全面审查 |
| 修改核心模块 (dgra_core.py, dgra_api.py) | **必须**调用专家 |
| 修复 Blocker 后 | 调用专家验证修复 |
| 不确定某写法是否合规 | 调用专家确认 |

---

## 防线三: 自动化门禁 (Git Hook)

### 安装 Git Pre-commit Hook

```bash
# 进入 skill 仓库
cd ~/.workbuddy/skills/dgra-genomic-risk

# 创建 pre-commit hook
cat > .git/hooks/pre-commit << 'EOF'
#!/bin/bash
# GPA Code Review Gate — Pre-commit Hook

echo "🔍 Running GPA Code Review Gate..."

python3 scripts/gpa_review_gate.py
EXIT_CODE=$?

if [ $EXIT_CODE -eq 1 ]; then
    echo ""
    echo "❌ COMMIT BLOCKED: Blocker issues found."
    echo "   Fix the issues above before committing."
    echo "   Run: python3 scripts/gpa_review_gate.py --all"
    exit 1
fi

if [ $EXIT_CODE -eq 2 ]; then
    echo ""
    echo "⚠️  Warning: Critical issues found."
    echo "   Consider fixing before committing."
    echo "   To bypass (not recommended): git commit --no-verify"
fi

exit 0
EOF

chmod +x .git/hooks/pre-commit
```

### 效果

每次 `git commit` 时：
- 自动运行 `gpa_review_gate.py`
- **Blocker 存在 → commit 被阻止**
- Critical 存在 → 警告但允许提交（需确认）
- 无问题 → 直接通过

### 临时绕过（紧急修复时）

```bash
git commit --no-verify  # 跳过 pre-commit hook，不推荐
```

---

## 防线三 (备选): GitHub Actions CI

### 配置 `.github/workflows/gpa-review.yml`

```yaml
name: GPA Code Review

on: [pull_request]

jobs:
  review-gate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Run GPA Review Gate
        run: |
          python3 scripts/gpa_review_gate.py --all
        working-directory: .

      - name: Check for Blockers
        if: failure()
        run: |
          echo "::error::Blocker issues found. Fix before merging."
          exit 1
```

---

## 审查触发清单

### 修改代码前必读

```markdown
## 修改前检查清单

- [ ] 我理解了要修改的功能和上下文
- [ ] 我已阅读相关函数的 docstring
- [ ] 我知道这个修改可能影响哪些其他模块

## 修改中遵守

- [ ] 未引入新的 `except Exception:`
- [ ] `open()` 都有 `encoding='utf-8'`
- [ ] 无可变默认参数
- [ ] 日志使用 lazy % 而非 f-string
- [ ] 新函数有类型注解

## 修改后验证

- [ ] 运行 `python3 scripts/gpa_review_gate.py`
- [ ] 所有 Blocker 已修复
- [ ] 相关测试通过
- [ ] 专家审查通过（如修改核心模块）
```

---

## 问题分级速查

| 问题 | 级别 | 自动检测 | 自动修复 |
|------|------|---------|---------|
| 语法错误 | 🔴 Blocker | ✅ py_compile | ❌ |
| bare except | 🔴 Blocker | ✅ AST | ❌ |
| except Exception | 🔴 Blocker | ✅ AST | ❌ |
| 可变默认参数 | 🔴 Blocker | ✅ AST | ❌ |
| open() 无 encoding | 🔴 Blocker | ✅ AST | ✅ 可脚本化 |
| eval/exec | 🔴 Blocker | ✅ AST | ❌ |
| 集合/字典重复项 | 🔴 Blocker | ✅ AST | ✅ 可脚本化 |
| `is` 用于非 None | 🔴 Blocker | ✅ AST | ❌ |
| 未使用导入 | 🟡 Critical | ✅ AST | ✅ 可脚本化 |
| 日志 f-string | 🟡 Critical | ✅ 正则 | ✅ 可脚本化 |
| 无意义 f-string | 💭 Nit | ✅ AST | ✅ 可脚本化 |
| 循环导入 | 🔴 Blocker | ✅ 自定义 | ❌ |
| 全局可变状态 | 🟡 Critical | ✅ AST | ❌ |
| 缺少类型注解 | 💭 Nit | ❌ | ❌ |

---

## 责任矩阵

| 角色 | 自检 | 专家审查 | 自动化门禁 |
|------|------|---------|-----------|
| 开发者 (你) | ✅ 必须 | 按需 | 被动触发 |
| AI 专家 | ❌ | ✅ 执行 | 可集成 |
| CI/Git Hook | ❌ | ❌ | ✅ 强制执行 |

---

## 运行记录

| 日期 | 版本 | 检查文件数 | Blocker | Critical | Nit |
|------|------|-----------|---------|----------|-----|
| 2026-05-30 | v0.10.10 | 35 | 121 | — | — |

---

## 附录: 快速修复命令

```bash
# 修复 open() 编码（批量）
find scripts -name "*.py" -exec sed -i '' "s/open(\([^,]*\))/open(\\1, encoding='utf-8')/g" {} \;

# 查找所有 except Exception
 grep -rn "except Exception:" scripts/

# 查找所有 bare except
 grep -rn "except:" scripts/ | grep -v "except Exception"

# 运行完整审查
 python3 scripts/gpa_review_gate.py --all
```
