# LDT-NSGA-II

轻量级分布迁移动态多目标优化器（**L**ightweight **D**istribution-**T**ransfer NSGA-II），
统一基于 pymoo 的 NSGA-II 内核实现。

本仓库是学位论文《面向动态多目标优化的轻量级分布迁移 NSGA-II：算法设计、实现与实验研究》
的**配套研究代码 + 全部实验数据**。基于**同一批数据**产出**两份学位论文**：

| 交付物 | 规格 | 与本文数据的关系 |
|---|---|---|
| **本科毕业论文** | 28 页 · 7 表 · 2 图 | 精简重写，保留全部核心结论与已知问题清单 |
| **硕士学位论文（实验论文）** | 177 页 · 51 表 · 11 图 · 1255 公式对象 | 含独立文献综述章与附录 A–E（附录 D 逐字嵌入本仓库源码） |
| **答辩演示 V2.5** | 15 页 · 9 张动画图表 | 数字全部由 `scripts/export_defense_data.py` 从 `results/` 导出 |

文档都在 [`../文档/`](../文档/) 目录下，本目录只放代码与数据。

- 当前代码版本：**v3.0**
- 实验数据：**6030 次运行记录，0 失败**（v2.4 四组 3780 次 + v2.5 两组 2250 次）
- 单元测试：**77 个，全部通过**（`python -m pytest -q`）
- CI：`.github/workflows/ci.yml`（Linux + Windows × Python 3.9/3.11）

> ⚠️ **v3.0 只改了工程可信度，没有改搜索行为。** 全部 v3.0 改动都是元数据、分组键、
> 路径、依赖与状态语义，**不动任何随机数消耗顺序**。已用 FDA4/τ=25/static 复现核对：
> 新代码 seeds 0–5 稳态 MIGD 均值 0.1037，与论文记录的 n=30 均值 0.1020（sd 0.0035）一致。
> 因此既有 6030 次运行数据**继续有效**，无需重跑。

---

## 目录

- [1. 这次更新了什么（v2.5 → v3.0）](#1-这次更新了什么v25--v30)
- [2. 项目文件夹地图：每个目录是干什么的](#2-项目文件夹地图每个目录是干什么的)
- [3. 数据可信度层：run_id / config_hash / 完整性](#3-数据可信度层run_id--config_hash--完整性)
- [4. 快速上手：五分钟跑通](#4-快速上手五分钟跑通)
- [5. 算法与实验设计](#5-算法与实验设计)
- [6. 实验结果（v3.0 口径不变）](#6-实验结果v30-口径不变)
- [7. 研究路线图：方案已规划、但尚未实现的部分](#7-研究路线图方案已规划但尚未实现的部分)
- [8. 文档在哪儿](#8-文档在哪儿)
- [9. 已知问题清单](#9-已知问题清单)
- [10. 复现实验的完整命令](#10-复现实验的完整命令)
- [11. 环境与版本说明](#11-环境与版本说明)

---

## 1. 这次更新了什么（v2.5 → v3.0）

v3.0 的改动全部来自《LDT-NSGA-II 项目优化与精进方案 V2.0》
（该方案基于**公开仓库 main 分支的静态审阅**写成；本工作区在若干项上已经领先，
下面逐项标注"本工作区原本已修 / v3.0 新修"）。**v3.0 不改变任何实验结果**，
只把"实验数据能不能被第三方追溯和复现"这件事做扎实。

### 1.1 P0 版本戳造假：逐代 CSV 的 `code_version` 是硬编码的

**这是本轮最严重的一个真实缺陷。** `algorithms/nsga2_core.py` 的 `DynamicDriver`
里写着 `self.code_version = "v2.4"`，而写进逐代 CSV 的正是这个字段——
**无论实际跑的是哪个版本的代码，CSV 都标 v2.4**。同时各驱动脚本又各自硬编码了
自己的 `code_version`（`run_comparison.py` 写的还是 `"v2.2"`）。

后果：论文分析脚本无法仅凭 `code_version` 区分 v2.2/v2.4/v2.5 数据。磁盘上
`results/summary_v25/integrity.json` 只出现 `v2.2` 与 `v2.4`、从不出现 `v2.5`，
原因就在这里。

**v3.0 新修**：

- `DynamicDriver.__init__` 增加 `code_version` 参数，**由运行配置注入**；
- 默认值为 `"unknown"`——缺失时明确表示"不知道"，**不再谎称 v2.4**；
- 四个驱动脚本统一改为 `code_version="v3.0"`；
- 逐代 CSV 新增 `run_id` 列，把每一行绑定到唯一一次运行。

> 历史数据的 `code_version` 戳保持原样（不重跑）。**因此判断一份旧 CSV 的真实版本，
> 应当依据其所在目录 + `run_meta.jsonl` + git 历史，而不是该字段。**

### 1.2 P0 运行唯一键：不同实验来源会被合并计数

`scripts/thesis_analysis.py` 的 `RUN_KEYS` 只有
`algorithm, problem, tau_t, population_size, n_t, seed`。敏感性实验里同一个
算法/问题/tau/N/n_t/seed 可以对应不同的 `alpha` 或 `rho`，
**仅凭这组字段做唯一键会把不同配置合并**。（本工作区此前已加 `_kind` 字段缓解，
但最小唯一键本身仍然不足。）

**v3.0 新修**：

- `RunConfig` 新增 `experiment_group`（主比较 / 补充 / 消融 / 敏感性 / alpha 专项）
  与显式的 `alpha` / `rho` 扫描槽位，不适用时记为 null 而非默认值；
- 新增 `RunConfig.config_dict()` 与 `config_hash()`（SHA-256 前 16 位）；
- 新增 `make_run_id()`：由 `(config_hash, algorithm, problem, seed)` 决定，
  **与时钟和随机源无关**，因此重跑同一命令得到同一个 id、重复可检测；
- `run_meta.jsonl` 写入完整 `config` 字典与 `run_id`。

### 1.3 P0 HTML 报告依赖本机路径

`build_report.py` 里写着 `ROOT = Path(r"D:\23-兴湘计科\LDT-NSGA-II")`，
换台机器或换个目录就构建不了。

**v3.0 新修**：改为 `ROOT = Path(__file__).resolve().parent`；
启动时校验 5 个必需输入（`doc.css`、`echarts.min.js`、`_chart_data.json`、
`_recovery.json`、`_boxplot.json`），缺失就**一次性列出全部缺失文件并以退出码 2 结束**，
而不是抛一个看不懂的 traceback 或悄悄生成一份没有图表的报告。
已验证：在仓库内直接 `python build_report.py` 正常产出 1.15 MB 自包含 HTML。

### 1.4 P1 补充实验默认只跑两个问题

`run_supplementary.py` 的 `--problems` 默认值是 `DMOP1,DMOP2`，而 README 的完整命令
传三个问题；`main.py` 调用时不传参数，**因此默认只覆盖两个问题**，
表面却像跑完了。

**v3.0 新修**：默认改为 `DMOP1,DMOP2,DMOP3` 并导出为模块级 `DEFAULT_PROBLEMS`；
运行开始前打印"预期组合数 = 问题数 × 算法数 × 种子数"；
结束时若 `complete` 为假则以**非零退出**。
项目变体 `DMOP1/2/3` 与权威实现 `DMOP1A/2A/3A` 继续严格区分。

### 1.5 P1 依赖清单与实际导入不一致

`requirements.txt` 没有 pandas，但 `evaluation/statistics.py` 与
`scripts/thesis_analysis.py` 都 `import pandas`——照 requirements 装完，分析路径就跑不起来。

**v3.0 新修**：补齐 pandas，并把全部依赖**固定到实测安装的版本**
（`pymoo==0.6.1.5`、`numpy==2.0.2`、`scipy==1.13.1`、`scikit-learn==1.6.1`、
`pandas==2.3.3`、`matplotlib==3.9.4`、`PyYAML==6.0.3`、`tqdm==4.70.1`、`pytest==8.4.2`）。

> 说明：这 9 个版本号是**逐个 import 实测**得到的，不是照抄上限。

### 1.6 P1 首次环境变化的 D_k 初始化语义

`ldmt_nsga2._detect_change()` 在 `_prev_probe_F is None` 时把当前探针值存下来
**并立即返回 0.0**。因为该方法首次调用发生在**环境已经切换之后**，
所以这个 0.0 把"无法比较"伪装成了"变化幅度为零"，污染了首次变化时的
自适应迁移宽度诊断（`adaptive_transfer` 为可选机制，默认关闭）。

**v3.0 新修**：

- 新增 `LDTNSGA2.establish_baseline(problem)`：在**每个环境开始时**（含环境 0）
  记录共享固定探针的目标值。探针是问题自有的固定点，不消耗决策变量、
  不影响搜索行为。
- `_detect_change()` 在无基线时返回 **`None`**，语义明确为"不可测量"，
  与"测量到、无变化"的 `0.0` 严格区分；`info["Dk"]` 相应可为 null。
- `adaptive_transfer` 遇到 `Dk is None` 时跳过该次统计，而不是按 0 计入均值。

### 1.7 P0/P1 完整性统计：截断的运行不再算"完成"

**v3.0 新修**（对应方案 5.2）：

- 一次完整运行应写 `n_gen` 行；`_existing_seeds()` 现在**只把行数达标的种子算完成**；
- `_partial_seeds()` 报告被截断的种子；
- `_drop_partial_seeds()` 在重跑前**清掉该种子的残缺行**，避免追加写出重复行；
- 汇总新增 `missing` / `partial` / `complete` 三栏，并把 `experiment_group`、
  `config_hash`、`code_version` 写进汇总；
- `run_comparison.py` / `run_supplementary.py` / `run_v25_alpha.py` 在
  `complete` 为假时**非零退出**；`run_ablation.py` 不再用裸 `continue` 吞掉失败，
  改为计数并写 `failed_runs.jsonl`。

### 1.8 P1 敏感性/消融实验也走统一元数据

`run_sensitivity.py` 以前**完全不写元数据**，并且直接用 `DynamicDriver`
（于是逐代行标 `unknown`）；`run_ablation.py` 也是自己一套续跑逻辑。
这正是"敏感性 tau25 与主比较 tau25 混在一起"的根源。

**v3.0 新修**：两个脚本都改用 `ExperimentRunner` 的完整性规则与元数据通道，
并把**被扫描的变量**（`sweep` / `N` / `n_t` / `alpha` / `rho`）显式写进
`experiment_group` 与 `config`，使每次运行可唯一识别。

### 1.9 P2 新增自动化 CI

新增 `.github/workflows/ci.yml`：矩阵 Linux + Windows × Python 3.9/3.11，
依次执行依赖安装 → 单元测试 → 权威 dMOP 自检 → 隔离目录冒烟实验 →
元数据链断言 → `analyze.py` + `build_report.py` → 报告自包含断言。
**全量长实验不进 CI**，冒烟结果写入 `results/_ci_smoke`（与正式数据隔离）。

### 1.10 新增 v3.0 回归测试

`tests/test_v30_fixes.py`（13 个测试）逐条锁定上面各项，包括：
版本戳必须是被注入的值、缺省必须是 `unknown`、`run_id` 必须逐行存在、
`config_hash` 必须随行为参数变化、不同实验组/扫描值必须得到不同 `run_id`、
截断种子不算完成且重跑前被清理、`build_report.py` 不得含盘符硬编码、
干净目录下必须给出可操作的缺文件报错、补充实验默认必须覆盖三个问题。

### 1.11 测试数量变化

63 → **77 个单元测试**（新增 13 个 v3.0 回归测试，另 1 个原有测试按新语义拆分为 2 个）。

---

## 2. 项目文件夹地图：每个目录是干什么的

### 2.1 顶层：两个目录，职责分明

```text
轻量级分布迁移动态多目标优化器/
├── LDT-NSGA-II/     ← 代码 + 实验数据（本目录）
└── 文档/             ← 论文、开题报告、答辩、审计报告
```

### 2.2 `LDT-NSGA-II/`（代码与数据）

| 目录 / 文件 | 作用 |
|---|---|
| **`algorithms/`** | 6 个算法 + 驱动器 + 共享数学设施 |
| ├ `nsga2_core.py` | `make_nsga2` 算子工厂（SBX 0.9/η15、PM η20、去重）、`DynamicDriver` 环境分段驱动、pymoo 评估计数兼容；**v3.0：`code_version` 参数注入 + 逐行 `run_id` + 每环境 `establish_baseline`** |
| ├ `ldmt_nsga2.py` | **LDT-NSGA-II 主算法**（归一化矩特征 → 趋势中心 → 高斯采样 → 50/20/30 重组）；含 `LDTNSGA2Alpha1`；**v3.0：`establish_baseline()` 与 `Dk=None` 语义** |
| ├ `distribution_features.py` | `normalize` / `denormalize` / `extract_features`（含 `SIGMA_MIN=1e-3`）/ 非支配集判定 |
| ├ `environment_manager.py` | 时间规则 + `ChangeDetector`（变化幅度 D_k） |
| ├ `transfer_strategy.py` | 可插拔策略层（默认实验未启用） |
| ├ `ar_nsga2.py` / `pps_nsga2.py` / `gmm_transfer_nsga2.py` | 三个预测/迁移基线（K1 修复后：先预留移民名额再拼接） |
| └ `dnsga2_response.py` | 静态基线 + DNSGA2 型响应基线 |
| **`problems/`** | 9 个测试问题 + 解析参考前沿 |
| ├ `base.py` | `DynamicProblem`：时间规则 + 共享固定探针（10 点、seed 2026、真实边界） |
| ├ `fda1.py` / `fda3.py` / `fda4.py` | FDA 系列（与 PlatEMO 逐式核对一致） |
| ├ `dmop1.py` / `dmop2.py` / `dmop3.py` | 三个**项目变体**（DMOP3 为静态对照，不得称标准 dMOP） |
| ├ `dmop1a.py` / `dmop2a.py` / `dmop3a.py` | **权威** dMOP1/2/3（覆盖 Type III/II/I） |
| └ `true_fronts.py` | 9 个问题的解析 Pareto 前沿；按 (问题, M, 点数, t) 缓存；FDA4 去重后 507 点 |
| **`evaluation/`** | 指标与统计 |
| ├ `igd.py` / `hv.py` | IGD、HV（HV 含容斥回退） |
| ├ `statistics.py` | Holm 逐步下降校正 + 配对/独立 Wilcoxon（**依赖 pandas**） |
| └ `friedman_test.py` | Friedman 卡方（平均秩）+ Nemenyi 双档 q 表与临界差 |
| **`experiments/`** | 实验驱动（**v3.0 全部统一元数据通道**） |
| ├ `common.py` | `RunConfig`（含 `experiment_group`/`alpha`/`rho`/`config_hash`）、`RunConfig.config_dict()`、`make_run_id()`、`git_commit()`、`git_dirty()`、`environment_meta()`、`ExperimentRunner`（完整性统计 + 断点续跑 + 截断清理） |
| ├ `run_comparison.py` | 主比较实验 |
| ├ `run_supplementary.py` | DMOP 补充实验（**默认三个项目变体**） |
| ├ `run_ablation.py` | 消融实验（4 个变体） |
| ├ `run_sensitivity.py` | 敏感性实验（N / n_t / ρ / α 四组扫描） |
| └ `run_v25_alpha.py` | α=1.0 对照 + 权威 dMOP 补充（两组原语料，现标 v3.0/组名） |
| **`scripts/`** | 分析与自检工具 |
| ├ `check_dmop_authoritative.py` | 权威 dMOP 解析自检 |
| ├ `thesis_analysis.py` | 论文级汇总与统计（从原始 CSV 从零重算） |
| ├ `generate_dmop_section.py` | 生成论文 5.9.2 节结果表 |
| ├ `export_defense_data.py` | 导出答辩演示所需的全部真实数据 |
| └ `make_thesis_figures.py` | 生成论文插图（`results/figures_thesis/`） |
| **`tests/`** | 10 个测试文件、**77 个单元测试**（含 `test_v30_fixes.py`） |
| **`.github/workflows/ci.yml`** | **v3.0 新增**：Linux+Windows × py3.9/3.11 质量门 |
| **`results/`** | **实验数据**（详见下表） |
| `configs/` | YAML 配置镜像（⚠️ **已过时且未被任何脚本读取**，见 K16） |
| `ldmt-experiment-report/` | `build_report.py` 生成的 HTML 报告 |
| `analyze.py` | v2.4 冻结汇总与图表（写 `results/summary/`、`results/figures/`） |
| `build_report.py` / `serve_report.py` | HTML 报告生成与本地预览 |
| `main.py` | 一键流水线（pytest → analyze → build_report → 打开报告） |
| `conftest.py` | pytest 配置（注册 `smoke` marker） |
| `requirements.txt` | 依赖清单（**v3.0 固定版本 + 补 pandas**） |

### 2.3 `results/` 各子目录（哪些数据参与论文统计）

| 目录 | 内容 | 参与论文统计？ |
|---|---|---|
| `raw/` | v2.4 四组实验逐代 CSV（127 文件）+ `run_meta.jsonl` | ✅ 主比较、补充、消融、敏感性 |
| `raw_v25/` | v2.5 两阶段实验逐代 CSV + `_logs_v25/` | ✅ α 专项、权威 dMOP |
| `raw_v25_badauthoritative/` | 首版**错误** dMOP 定义下的 12 个 CSV | ❌ 仅施工留痕 |
| `raw_legacy_v21/` | v2.1 时期旧数据（178 文件） | ❌ |
| `raw_legacy_v22_k1/` | K1 修复前（三基线移民被截断）的数据 | ❌ 归档备查 |
| `summary/` | v2.4 冻结汇总（由 `analyze.py` 生成，**不要手改**） | 交叉验证用 |
| `summary_v25/` | 论文级汇总（由 `scripts/thesis_analysis.py` 生成） | ✅ 论文表格直接来源 |
| `figures/` | 汇总图表（`analyze.py`） | 论文未直接使用 |
| `figures_thesis/` | 论文插图（11 张 PNG） | ✅ 正文插图 |

### 2.4 `../文档/`（论文与材料）

| 目录 / 文件 | 内容 |
|---|---|
| `论文材料/` | 两份学位论文 + 开题报告 + 答辩演示的 DOCX / Markdown / HTML |
| ├ `LDT-NSGA-II本科毕业论文.docx` / `.md` | **本科毕业论文（28 页 / 7 表 / 2 图）** |
| ├ `LDT-NSGA-II实验论文v3.0.docx` / `.md` | **硕士（实验）论文 v3.0** |
| ├ `LDT-NSGA-II开题报告.docx` / `.md` | 开题报告（29 页 / 12 表） |
| ├ `LDT-NSGA-II V2.5答辩.html` | 答辩演示（15 页 / 9 张动画图表） |
| ├ `_源文件/` | 文档构建系统：分章节 Markdown 源 + 构建脚本 |
| ├ `_源文件/_本科版/` | 本科版的独立入口与分章节源 |
| ├ `_audit/` | 独立审计与复核报告（代码审计、数值复核、文献核查、dMOP 定义考证） |
| └ `参考归档/` | 中期检查材料：V2.5 论文与旧版、v2.4 草稿、开发日志、精进方案 |

**改论文的正确姿势**：编辑 `文档/论文材料/__源文件/chapterN_*.md`，然后运行
`python 文档/论文材料/__源文件/build_docx.py` 重新生成 DOCX 与交付用 Markdown。
本科版改 `_源文件/_本科版/*.md`，跑 `_源文件/_本科版/build_ug_docx.py`。

---

## 3. 数据可信度层：run_id / config_hash / 完整性

v3.0 的核心增量。**每条独立运行在启动时获得一个不可变 `run_id`**：

```text
run_id = sha256(config_hash + algorithm + problem + seed)[:20]
config_hash = sha256(规范化 config_dict)[:16]
```

`config_dict` 覆盖全部会影响行为的参数：问题、`tau_t`、`n_t`、`N`、`n_gen`、
算法列表、种子区间、`code_version`、`experiment_group`、`alpha`、`rho` 及 `extra`。
因为不含时钟与随机源，**重跑同一命令得到同一个 `run_id`**——重复可以直接检测。

`run_meta.jsonl` 每次运行追加一条记录：

| 字段 | 含义 |
|---|---|
| `run_id` / `config_hash` / `config` | 身份与完整配置 |
| `experiment_group` | 主比较 / 补充 / 消融 / 敏感性 / alpha 专项 |
| `algorithm` / `problem` / `tau_t` / `n_t` / `N` / `n_gen` / `seed` | 运行坐标 |
| `n_gen_rows` | 实际写出的逐代行数（完整性） |
| `n_eval_total` / `n_eval_per_env` | 函数评估次数（预算公平性） |
| `git_commit` / `git_dirty` | 代码版本；**dirty 树无法追溯到单一提交，如实记录** |
| `code_version` / `status` / `environment` | 版本、状态、Python 与依赖版本 |

逐代 CSV 通过 `run_id` 与元数据关联，并新增 `run_id` 列。

**完整性判定**：一次完整运行应为 `n_gen` 行。行数不足的种子**不算完成**、
会在重跑前被清理掉，汇总中计入 `partial`；完全没有数据的计入 `missing`。
`complete` 为假时驱动脚本非零退出。

**一句话**：以后拿任意一行 CSV，都能沿 `run_id → run_meta.jsonl → config + git_commit`
查回它是哪次运行、跑的什么配置、哪份代码、什么环境。

---

## 4. 快速上手：五分钟跑通

```powershell
cd LDT-NSGA-II

# 1) 安装依赖（版本已固定；实测 Python 3.9.13）
python -m pip install -r requirements.txt

# 2) 77 个单元测试（约 4 秒）
python -m pytest -q

# 3) 权威 dMOP 解析自检（应输出 FINAL: ALL CHECKS PASSED）
python -m scripts.check_dmop_authoritative

# 4) 单次冒烟（分钟级；务必写入独立目录！）
python -m experiments.run_comparison --problems FDA3 --taus 25 --seeds 0 `
    --out results/_smoke --n-gen 50 --limit-seeds 1

# 5) HTML 报告（v3.0：路径全部相对脚本自身，可干净克隆构建）
python analyze.py
python build_report.py
```

> ⚠️ **两条必须知道的坑**
> 1. 冒烟数据用的是缩小预算（50 代而非 300 代），一旦写进 `results/raw` 就会**污染正式统计**。
>    v3.0 的完整性检查会把这类文件里的种子判为 `partial` 并要求重跑——这是刻意的。
> 2. v3.0 起 `run_supplementary.py` 默认跑三个项目变体；若要复现论文的历史 540 次，
>    请显式传 `--problems DMOP1,DMOP2,DMOP3`。

---

## 5. 算法与实验设计

### 5.1 LDT-NSGA-II 核心

把 DMOP 求解显式拆成**内核**与**策略**两个正交部分：内核是未改动的 pymoo NSGA-II；
策略实现 `respond(prev_X, prev_F, problem, env_index, rng)` 接口，在环境切换时生成新环境初始种群。

设上一环境末种群归一化后为 `Z ∈ [0,1]^{N×D}`，取其非支配子集（若规模 `< max(10, N/4)` 则回退全种群）：

```
μe = mean(ND)                    逐维均值
σe = std(ND, ddof=0)  且 σe ≥ 1e-3   逐维标准差（归一化单位下限）
μ̂  = clip(μe + α·(μe − μprev), 0, 1)   趋势中心外推，默认 α = 0.5
```

新环境初始种群按 **50% 迁移 + 20% 旧解保留 + 30% 随机移民** 拼接：

1. 迁移个体：`x_norm = clip(μ̂ + z ⊙ σe, 0, 1)`，`z ~ N(0, I)`，**反归一化恰好一次**；
2. 旧解保留：从上一种群随机抽取（`rng.permutation`）；
3. 随机移民：真实边界上的均匀采样，作为预测失准时的多样性兜底。

**"轻量"的量化**：单次响应复杂度 `O(MN² + ND)`，分布只压缩为 `2D` 个标量
（完整高斯混合模型需要 `O(KD²)` 个参数），无矩阵分解、无迭代优化。

两个**默认关闭**的可选机制：变化幅度驱动的自适应迁移宽度缩放
（`clip(1+β|D_k|/mean|D_k|, 0.5, 2.0)`）与多源历史二阶矩融合（历史 3、权重下限 0.05）。
v3.0 修正了前者的 D_k 初始化语义（见 1.6）。

### 5.2 六个算法（仅 `respond` 不同，其余全部相同）

| 标识 | 变化时的行为 |
|---|---|
| `nsga2_static` | 原样带入上一种群，无任何响应 |
| `dnsga2_inspired` | 保留种群前缀 80% + 20% 均匀随机移民 |
| `ar1_predictive` | 逐维 OLS 一阶自回归外推（历史 5）+ 20% 移民；**方差不预测** |
| `gmm_inspired` | 高斯混合（≤3 分量，`reg_covar=1e-4`）采样 + **10%** 移民 |
| `pps` | 中心全步差分 `2μ_{k−1}−μ_{k−2}` + 整体平移 + 20% 移民 |
| `ldt_nsga2` | 半步趋势外推 + 高斯采样迁移 + 50/20/30 重组（**本文方法**） |
| `ldt_alpha1` | 同上，但 `α = 1.0`（机制对照变体） |

> 命名诚实性：`dnsga2_inspired` / `gmm_inspired` 是"受启发"的项目实现，**不是**原文复现；
> `pps` 是 PPS 的简化实现（无独立方差/流形预测项）。GMM 的移民比例是 10% 而非 20%（源码既有事实）。

### 5.3 九个测试问题

| 问题 | M | D | 变量边界 | 类型 | 说明 |
|---|---|---|---|---|---|
| FDA1 | 2 | 10 | x₁∈[0,1]，其余 [−1,1] | I | 解集平移，前沿时不变 |
| FDA3 | 2 | 10 | 同上 | II | 解集与前沿均变 |
| FDA4 | 3 | 12 | 全部 [0,1] | I | 三目标球面前沿，距离变量切片 `X[:, M−1:]` |
| DMOP1（变体） | 2 | 10 | [0,1] | — | 项目变体（系数 1、平方前沿） |
| DMOP2（变体） | 2 | 10 | [0,1] | — | 项目变体（含 +G 项、指数固定 2） |
| DMOP3（变体） | 2 | 10 | [0,1] | **静态** | 静态对照组，非动态问题 |
| **DMOP1A** | 2 | 10 | [0,1] | **III** | 权威 dMOP1：**解集静止、前沿旋转** |
| **DMOP2A** | 2 | 10 | [0,1] | **II** | 权威 dMOP2：解集与前沿均变 |
| **DMOP3A** | 2 | 10 | [0,1] | **I** | 权威 dMOP3：`r` 按环境轮换，解集移动 |

> dMOP1–dMOP3 的权威定义依据四份独立复现文献交叉确认（Goh & Tan 2009 原文受付费墙限制、
> 无开放获取版本，**本项目未能直接阅读原始公式**）。dMOP3 的 `G(t)` 存在 `sin` 与 `|sin|`
> 两种文献写法，本实现采用后者并已显式声明。

### 5.4 统计检验

- 两算法：配对 **Wilcoxon 符号秩检验**（种子完全对齐）+ **Holm** 逐步下降校正；
- 多算法：**Friedman** 检验 + **Nemenyi** 事后检验，`k=6, N=30, α=0.05` 时 `CD = 1.37668`；
- 主指标口径：**稳态 MIGD** —— 每个环境段最后一代 IGD 的段间平均，再对种子平均。
  另有"变化后 MIGD"（段首代）与"全代 MHV"两个辅助口径。

---

## 6. 实验结果（v3.0 口径不变）

> **v3.0 不重跑任何实验。** 本节数字与 v2.5 完全一致，因为它来自同一批冻结数据。

### 6.1 规模与状态

| 实验组 | 因子设计 | 运行记录数 | 失败 |
|---|---|---|---|
| v2.4 主比较 | FDA1/3/4 × τ{10,25,50} × 6 算法 × 30 种子 | 1620 | 0 |
| v2.4 补充 | DMOP1/2/3 × τ25 × 6 算法 × 30 种子 | 540 | 0 |
| v2.4 消融 | 3 问题 × 4 变体 × 30 种子 | 360 | 0 |
| v2.4 敏感性 | N / n_t / ρ / α 四组扫描 | 1260 | 0 |
| v2.5 阶段一 | FDA1/3/4 × τ{10,25,50} × 6 算法 × 30 种子（含 `ldt_alpha1`） | 1620 | 0 |
| v2.5 阶段二 | 权威 dMOP1A/2A/3A × τ25 × 6 算法 × 30 种子 | 540 | 0 |
| v2.5 阶段二附加 | 权威 dMOP × `ldt_nsga2`（α=0.5）× 30 种子 | 90 | 0 |
| **合计** | | **6030** | **0** |

### 6.2 核心结论

1. **步长归因（主结论）**：α=1.0 在 12 个配置中的 11 个上更优，
   与 PPS 差距缩到 `-8.5%~+3.4%`；**2 个配置反超**（FDA4 τ=10 / τ=25，p 分别 1.86e-9 / 5.72e-7）。
   → 原"LDT 显著劣于 PPS"主要来自外推步长，而非迁移范式。
2. **反例与自适应规则**：唯一例外是 **DMOP1A（Type III，解集静止）**，
   该问题 α=0.5 领先 α=1.0 达 **127%**（0.0251 vs 0.0570）。
   → 外推步长应随"解集是否在显著移动"自适应。
   诚实补充：该格上两者 MHV 几乎相同（0.5702 vs 0.5703），差异主要反映 IGD 口径下的收敛速度。
3. **组件贡献**：分布迁移（去掉后退化 13–30%）≫ 趋势外推（12–19%）≫ 旧种群保留（无一致增益）。
4. **计算开销**：GMM 单次响应耗时约为同配置 LDT 的 **4.4–11 倍**且性能不占优；
   LDT 与 AR(1)、PPS 同处 10⁻² s 量级。
5. **基线结论**：PPS 在 9 个主比较组合均值第一、LDT（α=0.5）稳居第二；
   LDT 显著优于 static/dnsga2/gmm；在 Friedman–Nemenyi 整体排序下与 PPS 于 4 个配置不可区分。

### 6.3 图表与汇总

- 论文插图：`results/figures_thesis/`（11 张）
- 汇总图表：`results/figures/`（箱线图、MIGD–τ 曲线、恢复曲线、消融图、敏感性图等）
- 论文级汇总：`results/summary_v25/`（`analysis_report.md` 可直接粘贴）
- v2.4 冻结汇总：`results/summary/`

---

## 7. 研究路线图：方案已规划、但尚未实现的部分

**必须明确：以下机制在 v3.0 中并未实现。** 优化方案的阶段 D
（"新算法机制"）属于研究建议，需要独立实验验证后才能写进论文结论。
v3.0 只完成了方案中**不影响搜索行为**的阶段 A/B/C，正是为了给阶段 D 铺好可追溯的实验底座。

| 方案条目 | 内容 | v3.0 状态 |
|---|---|---|
| 3.1 | **自适应趋势步长**：按"解集位移证据 vs 估计噪声"在线调节 α_k | ❌ **未实现**（下一优先项） |
| 3.2 | **相关性建模**：LDT-Cov（shrinkage 协方差）/ LDT-LowRank（PCA 子空间） | ❌ 未实现 |
| 3.3 | **迁移门控**：按历史预测误差调节迁移预算 ρ_k，并把候选验证评估计入总预算 | ❌ 未实现 |
| 3.4 | 多源融合状态口径：区分 `mu_current_raw` / `mu_current_fused` | ⚠️ 仅记录为已知问题（K13），**未修**；该机制默认关闭，不影响主结论 |
| 4.1 | 分离"已知变化时刻"与"未知变化在线检测"（检测延迟/误报/漏报） | ❌ 未实现 |
| 4.2 | 扩展变化模式（非等间隔、方向反转、突变/渐变、POF 分离、噪声） | ❌ 部分覆盖（Type I/II/III 已覆盖） |
| 4.3 | 统一参考前沿、多分辨率收敛检查、新增"变化后前 w 代 IGD 曲线面积"指标 | ⚠️ 参考前沿已固定并有测试；多分辨率检查与新指标未实现 |
| 5.3 | `configs/*.yaml` 唯一配置入口的取舍 | ⚠️ 仍标注为过时（K16）；当前唯一配置来源是 CLI 参数 |

**方案本身给出的边界提醒同样适用**：静态审阅不能证明程序在目标机器上必然跑通；
算法增强方向是研究建议，效果、耗时与统计显著性都需实验验证；
不要为了增加创新点同时堆叠太多机制，应优先用一条主线完成
"假设—设计—对照—消融—统计—局限"的闭环。

---

## 8. 文档在哪儿

| 想要什么 | 去哪儿 |
|---|---|
| **本科毕业论文（28 页）** | [`../文档/论文材料/LDT-NSGA-II本科毕业论文.docx`](../文档/论文材料/) |
| **硕士（实验）论文 v3.0** | [`../文档/论文材料/LDT-NSGA-II实验论文v3.0.docx`](../文档/论文材料/) |
| 开题报告 DOCX | [`../文档/论文材料/LDT-NSGA-II开题报告.docx`](../文档/论文材料/) |
| **答辩演示（含动画图表）** | [`../文档/论文材料/LDT-NSGA-II V2.5答辩.html`](../文档/论文材料/) |
| 改论文内容 | `../文档/论文材料/__源文件/chapterN_*.md`（本科版在 `_本科版/`），然后跑对应构建脚本 |
| 重做答辩演示 | `python -m scripts.export_defense_data` → `python ../文档/论文材料/__源文件/make_defense.py` |
| 零基础运行手册 | 论文附录 A（与本文第 4 节一致） |
| 全部核心源码 | 论文附录 D（由构建脚本从本仓库原文件逐字嵌入，不会与代码脱节） |
| 独立审计/复核报告 | `../文档/论文材料/_audit/` |
| 中期检查材料（V2.5 论文等） | `../文档/论文材料/参考归档/` |
| 项目优化方案 V2.0 | `../文档/论文材料/参考归档/LDT-NSGA-II项目优化与精进方案 V2.0.docx` |

**重建文档**：

```powershell
# 硕士（实验）论文 v3.0
python "..\文档\论文材料\_源文件\build_docx.py"

# 本科毕业论文（28 页；页数由脚本里传给 finalize_docx.py 的行距系数控制）
python "..\文档\论文材料\_源文件\_本科版\build_ug_docx.py"

# 开题报告（29 页）
python "..\文档\论文材料\_源文件\build_proposal_docx.py"

# 答辩演示（先从 results/ 导出真实数据，再生成自包含 HTML）
python -m scripts.export_defense_data
python "..\文档\论文材料\_源文件\make_defense.py"
```

> 两份学位论文共用同一套样式模板（`_源文件/reference-thesis.docx`）与同一条 finalize
> 流水线（分页、目录定位、页眉页码、Word 域更新），因此排版风格一致、结论完全同源。

---

## 9. 已知问题清单

论文以 **K1–K33** 完整披露。摘要如下（全部**不影响第 6 章结论**）：

| 编号 | 问题 | 状态 |
|---|---|---|
| K1 | AR1/GMM/PPS 的 `immigrant_ratio` 曾是死参数（移民被定容截断丢弃） | ✅ 已修复并重跑 1080 次 |
| K2 / K3 | DMOP3 实为静态；DMOP2 是项目变体非标准 dMOP2 | ✅ 已核实并如实标注 |
| K4–K9 | 坐标系混用、HV 回退 3 处 bug、无效字段、Friedman 并列秩、损坏 `main.py`、诊断列硬编码 | ✅ 已修复（0 次重跑） |
| K10 / K11 | FDA4 前沿实测 507 点；PlatEMO 无 DMOP*.m | ✅ 已核实 |
| K12 | HTML 报告非自包含 | ✅ 已修复 |
| **K13** | 多源融合时趋势项量纲口径不一致（`mu_eff` vs `mu_e`） | ⚠️ **如实保留**（默认未启用；方案 3.4 列为待办） |
| K14 | `run_supplementary.py` 默认只跑 DMOP1,DMOP2 | ✅ **v3.0 已修**（默认三问题 + 非零退出） |
| K15 | `build_report.py` 生成的 HTML 报告叙述文字是旧版硬编码 | ⚠️ 待修（表格数据正确） |
| K16 | `configs/*.yaml` 已过时且未被任何脚本读取 | ⚠️ 待清理（当前唯一配置来源是 CLI） |
| K17 | README 曾写 `sigma_min = 0.001*(u−l)`（真实尺度），代码是归一化 1e-3 | ✅ 本版已更正 |
| K18–K21 | 数值复核发现的 4 处不一致（含 GMM 移民实为 10%、"3780 次独立运行"措辞不准） | ✅ 已修正 |
| K22–K26 | 权威 dMOP1/2/3 定义错误、数据补齐、时变前沿判定、`get_true_pf` 重构 | ✅ 已修复并重跑 |
| **K27** | **逐代 CSV 的 `code_version` 是硬编码**，无法据此区分代码版本 | ✅ **v3.0 已修**（注入 + `run_id`）；历史戳保持原样 |
| **K28** | **统计分组键不含实验来源与扫描变量**，不同实验会被合并计数 | ✅ **v3.0 已修**（`experiment_group` + `config_hash` + `run_id`） |
| **K29** | **`build_report.py` 硬编码本机绝对路径**，干净克隆无法构建 | ✅ **v3.0 已修** |
| **K30** | **`requirements.txt` 缺 pandas**，分析路径在干净环境失败 | ✅ **v3.0 已修** |
| **K31** | **首次环境变化的 D_k 被记为 0.0**，把"不可测量"伪装成"无变化" | ✅ **v3.0 已修**（`Dk=None` + 每环境基线） |
| **K32** | **截断的运行被算作"已完成"**，续跑会追加出重复行 | ✅ **v3.0 已修**（`partial` 判定 + 重跑前清理） |
| **K33** | **无 CI 质量门** | ✅ **v3.0 已新增**（Linux+Windows × py3.9/3.11） |

---

## 10. 复现实验的完整命令

```powershell
cd LDT-NSGA-II

# ---- v2.4 四组实验（写入 results/raw，约数小时，可中断续跑）----
python -m experiments.run_comparison --problems FDA1,FDA3,FDA4 --taus 10,25,50 `
    --seeds 0-29 --out results/raw
python -m experiments.run_supplementary --problems DMOP1,DMOP2,DMOP3 --seeds 0-29 `
    --out results/raw --algos nsga2_static,dnsga2_inspired,ar1_predictive,gmm_inspired,pps,ldt_nsga2
python -m experiments.run_ablation
python -m experiments.run_sensitivity

# ---- v2.5 两阶段（写入 results/raw_v25，与 v2.4 严格隔离）----
python -m experiments.run_v25_alpha --problems FDA1,FDA3,FDA4 --taus 10,25,50 `
    --seeds 0-29 --out results/raw_v25 --phase comparison
python -m experiments.run_v25_alpha --problems DMOP1A,DMOP2A,DMOP3A --taus 25 `
    --seeds 0-29 --out results/raw_v25 --phase supplementary
python -m experiments.run_v25_alpha --problems DMOP1A,DMOP2A,DMOP3A --taus 25 `
    --seeds 0-29 --out results/raw_v25 --phase supplementary --algos ldt_nsga2

# ---- 分析与自检 ----
python -m scripts.check_dmop_authoritative     # 权威 dMOP 解析自检
python -m scripts.thesis_analysis              # 论文级汇总 → results/summary_v25/
python -m scripts.generate_dmop_section        # 生成 5.9.2 节结果表
python -m scripts.make_thesis_figures          # 论文插图
python -m scripts.export_defense_data          # 答辩演示数据
python analyze.py                              # v2.4 冻结汇总 + 图表
python build_report.py                         # 自包含 HTML 报告
python main.py                                 # 一键：测试 → 汇总 → HTML 报告
```

> **注意**：v3.0 起每个驱动脚本都会在结束时打印
> `[SUMMARY] group=... expected=... completed=... failed=... missing=... partial=... complete=...`，
> 并在 `complete=False` 时以非零码退出。CI 就是这个约定的消费者。

> **并行加速**：`run_v25_alpha.py` 可按 (问题, τ) 拆成 9 个独立进程并行；
> 在 32 逻辑核机器上把阶段一从约 5 小时压到约 0.5 小时。`run_v25_comparison.ps1`
> 是当时用的并行启动脚本（含日志目录与 PID 记录），可参考。

---

## 11. 环境与版本说明

| 项 | 值 |
|---|---|
| Python | 3.9.13（实测） |
| pymoo | 0.6.1.5（⚠️ 0.6.x 无 `algo.n_eval`，真实计数器在 `algo.evaluator.n_eval`，见 `_algo_n_eval()`） |
| numpy / scipy / scikit-learn | 2.0.2 / 1.13.1 / 1.6.1 |
| pandas / matplotlib / PyYAML / tqdm / pytest | 2.3.3 / 3.9.4 / 6.0.3 / 4.70.1 / 8.4.2 |
| CI 矩阵 | ubuntu-latest、windows-latest × Python 3.9、3.11 |
| **数据版本戳** | 历史 v2.4 数据的 `code_version` 标为 `v2.2` 或 `v2.4`；v2.5 数据标 `v2.5`；**v3.0 起由配置注入**。该字段在 v3.0 前不可信（K27），判断旧数据版本请依据目录 + `run_meta.jsonl` + git 历史 |
