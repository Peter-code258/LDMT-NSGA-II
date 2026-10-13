# LDT-NSGA-II

轻量级分布迁移动态多目标优化器（**L**ightweight **D**istribution-**T**ransfer NSGA-II），
统一基于 pymoo 的 NSGA-II 内核实现。

本仓库是学位论文《面向动态多目标优化的轻量级分布迁移 NSGA-II：算法设计、实现与实验研究》
的**配套研究代码 + 全部实验数据**。基于**同一批数据**产出**两份学位论文**：

| 交付物 | 规格 | 与本文数据的关系 |
|---|---|---|
| **本科毕业论文** | 28 页 · 7 表 · 2 图 | 精简重写，保留全部核心结论与已知问题清单 |
| **硕士学位论文（详细版）** | 177 页 · 51 表 · 11 图 · 1255 公式对象 | 含独立文献综述章与附录 A–E（附录 D 逐字嵌入本仓库源码） |
| **答辩演示 V2.5** | 15 页 · 9 张动画图表 | 数字全部由 `scripts/export_defense_data.py` 从 `results/` 导出 |

文档都在 [`../文档/`](../文档/) 目录下，本目录只放代码与数据。

- 当前代码版本：**v2.5**
- 已完成实验：**6030 次运行记录，0 失败**（v2.4 四组 3780 次 + v2.5 两组 2250 次）
- 单元测试：**63 个，全部通过**（`python -m pytest -q`）

---

## 目录

- [1. 这次更新了什么（v2.4 → v2.5）](#1-这次更新了什么v24--v25)
- [2. 项目文件夹地图：每个目录是干什么的](#2-项目文件夹地图每个目录是干什么的)
- [3. 快速上手：五分钟跑通](#3-快速上手五分钟跑通)
- [4. 算法与实验设计](#4-算法与实验设计)
- [5. 实验结果（v2.5 最新）](#5-实验结果v25-最新)
- [6. 文档在哪儿](#6-文档在哪儿)
- [7. 已知问题清单](#7-已知问题清单)
- [8. 复现实验的完整命令](#8-复现实验的完整命令)

---

## 1. 这次更新了什么（v2.4 → v2.5）

这是一次**有实质结论变化**的更新，不只是整理文件。按重要性排序：

### 1.1 新增 α=1.0 变体对照实验（最重要的结论更新）

v2.4 的结论是"简化版 PPS 全面第一，LDT 稳居第二且显著劣于 PPS"。v2.5 设计了
**"外推步长 × 迁移范式"解耦实验**来追问差距从哪来：

- 新增算法 `ldt_alpha1`（`algorithms/ldmt_nsga2.py` 中的 `LDTNSGA2Alpha1`），
  即把趋势系数从默认半步 `α=0.5` 改为全步 `α=1.0`，其余完全不变；
- 用**同协议、同批次**重跑 5 个基线 + `ldt_alpha1`，共 1620 次；
- 结果：**α=1.0 在 12 个"问题 × τ"配置中的 11 个上更优**，与 PPS 的差距从
  `+1.4%~+33.9%` 缩小到 `-8.5%~+3.4%`，其中 4 个配置统计不可区分、
  **2 个配置显著反超**（FDA4 τ=10：0.1068 vs 0.1167，p=1.86e-9）。

→ **结论修正**：v2.4 观察到的劣势主要来自**外推步长**，不是"矩迁移 vs 点预测"的范式差异。

### 1.2 新增并更正权威 dMOP1/2/3 实现（纠正了两处定义错误）

新增 `problems/dmop1a.py` / `dmop2a.py` / `dmop3a.py`，并新增独立自检脚本
`scripts/check_dmop_authoritative.py`。经四份独立复现文献交叉核对后，**更正了两处此前的错误**：

| 问题 | 原（错） | 现（对） |
|---|---|---|
| dMOP1 | Type I，`g` 含 `G(t)`，前沿 `√` 型 | **Type III**（解集静止、前沿随 `H(t)` 旋转），`g = 1+9Σx_i²`（不含 G(t)），前沿 `1-u^{H(t)}` |
| dMOP3 | 项目变体的欺骗性前沿 + `H(t)` 解集位置 | Type I，`f₁ = x_r`（r 按环境轮换），`g = 1+9Σ_{i≠r}(x_i−G)²`，前沿 `1−√f₁` |

受影响的 360 次运行作废并归档到 `results/raw_v25_badauthoritative/`，修正后重跑。
`_row_metrics` 的"按 t 取前沿"分支也一并补上了 `DMOP1A`（Type III 的前沿同样随 t 变）。

### 1.3 新增"解集静止时外推有害"的反例与自适应步长规则

补齐 `ldt_nsga2`（α=0.5）在权威 dMOP 上的 90 次运行后，发现**唯一例外**：
在 **DMOP1A（Type III，Pareto 解集静止、只有前沿旋转）** 上，
α=0.5 以 0.0251 对 0.0570 **大幅领先** α=1.0（p=1.86e-9，相对改进 127%）。

→ 机理：解集不动时中心增量 `μe−μprev` 几乎只含采样噪声，全步外推把噪声放大成系统性偏移。
→ 由此提出可操作规则：**外推步长应随"解集是否在显著移动"自适应**；12 个配置已把该规则两端都验证过。

### 1.4 新增脚本

| 文件 | 作用 |
|---|---|
| `experiments/run_v25_alpha.py` | v2.5 两阶段实验驱动；内置**拒绝写入冻结目录**的保护（防止污染 v2.4 数据） |
| `scripts/check_dmop_authoritative.py` | 权威 dMOP 自检：解析 PS 落在解析 PF 上的偏差、问题类型、`r` 切换策略 |
| `scripts/thesis_analysis.py` | **论文级**汇总与统计，从原始逐代 CSV 从零重算全部表格（独立于 `analyze.py`，可交叉验证） |
| `scripts/generate_dmop_section.py` | 生成论文 5.9.2 节的权威 dMOP 结果表（数字由脚本产出，不手抄） |

### 1.5 修复的缺陷（K14–K26）

- `get_true_pf` 的默认目标数判定由"按问题名硬编码"改为
  `M = 3 if name in ("FDA4","FDA5") else 2`（K26）；
- 论文写作期间的独立数值复核比对了 **269 项断言**，发现并修正 4 处不一致
  （K18–K21），其中包括 **GMM 基线的随机移民实际是 10% 而非 20%**；
- 补跑缺失数据（K24）、更正 `_row_metrics` 时变前沿判定（K25）。

> 完整清单（K1–K26，含每项的"已修复/已核实/如实披露"状态）见论文第 5.11 节。

### 1.6 目录结构整理

按"文档归 `文档/`、代码归 `LDT-NSGA-II/`"重新组织：

- 论文与开题报告的 Markdown、DOCX **构建系统**（原 `_build_thesis/`）移到
  `../文档/论文材料/_源文件/`；
- 生成的论文/开题报告 Markdown 移到 `../文档/论文材料/`；
- 删除失效的开发期脚本 `_check_pilot.py`、`_run_batch2.ps1`，以及可再生的
  `__pycache__/`、`.pytest_cache/`。

---

## 2. 项目文件夹地图：每个目录是干什么的

### 2.1 顶层：两个目录，职责分明

```text
轻量级分布迁移动态多目标优化器/
├── LDT-NSGA-II/     ← 代码 + 实验数据（本目录）
└── 文档/             ← 论文、开题报告、答辩、审计报告（日志已归入 _参考归档/）
```

### 2.2 `LDT-NSGA-II/`（代码与数据）

| 目录 / 文件 | 作用 |
|---|---|
| **`algorithms/`** | 6 个算法 + 驱动器 + 共享数学设施 |
| ├ `nsga2_core.py` | `make_nsga2` 算子工厂（SBX 0.9/η15、PM η20、去重）、`DynamicDriver` 环境分段驱动、pymoo 评估计数兼容 |
| ├ `ldmt_nsga2.py` | **LDT-NSGA-II 主算法**（归一化矩特征 → 趋势中心 → 高斯采样 → 50/20/30 重组）；文件末尾含 v2.5 的 `LDTNSGA2Alpha1` |
| ├ `distribution_features.py` | `normalize` / `denormalize` / `extract_features`（含 `SIGMA_MIN=1e-3`）/ 非支配集判定 |
| ├ `environment_manager.py` | 时间规则 + `ChangeDetector`（变化幅度 D_k） |
| ├ `transfer_strategy.py` | 可插拔策略层（默认实验未启用） |
| ├ `ar_nsga2.py` / `pps_nsga2.py` / `gmm_transfer_nsga2.py` | 三个预测/迁移基线（K1 修复后：先预留移民名额再拼接） |
| └ `dnsga2_response.py` | 静态基线 + DNSGA-II 型响应基线 |
| **`problems/`** | 9 个测试问题 + 解析参考前沿 |
| ├ `base.py` | `DynamicProblem`：时间规则 + 共享固定探针（10 点、seed 2026、真实边界） |
| ├ `fda1.py` / `fda3.py` / `fda4.py` | FDA 系列（与 PlatEMO 逐式核对一致） |
| ├ `dmop1.py` / `dmop2.py` / `dmop3.py` | 三个**项目变体**（DMOP3 为静态对照，不得称标准 dMOP） |
| ├ `dmop1a.py` / `dmop2a.py` / `dmop3a.py` | **权威** dMOP1/2/3（v2.5 新增，覆盖 Type III/II/I） |
| └ `true_fronts.py` | 9 个问题的解析 Pareto 前沿；按 (问题, M, 点数, t) 缓存；FDA4 去重后 507 点 |
| **`evaluation/`** | 指标与统计 |
| ├ `igd.py` / `hv.py` | IGD、HV（HV 含 K5 修复后的容斥回退） |
| ├ `statistics.py` | Holm 逐步下降校正 + 配对/独立 Wilcoxon |
| └ `friedman_test.py` | Friedman 卡方（平均秩）+ Nemenyi 双档 q 表与临界差 |
| **`experiments/`** | 实验驱动 |
| ├ `common.py` | `RunConfig`、算法/问题注册表、`ExperimentRunner`（断点续跑 + 运行级元数据） |
| ├ `run_comparison.py` | v2.4 主比较实验 |
| ├ `run_supplementary.py` | v2.4 DMOP 补充实验（⚠️ 默认参数只跑 DMOP1,DMOP2，见 K14） |
| ├ `run_ablation.py` | 消融实验（4 个变体） |
| ├ `run_sensitivity.py` | 敏感性实验（N / n_t / ρ / α 四组扫描） |
| └ `run_v25_alpha.py` | **v2.5 两阶段实验**（α=1.0 对照 + 权威 dMOP 补充） |
| **`scripts/`** | 分析与自检工具 |
| ├ `check_dmop_authoritative.py` | 权威 dMOP 解析自检 |
| ├ `thesis_analysis.py` | 论文级汇总与统计（从原始 CSV 从零重算） |
| ├ `generate_dmop_section.py` | 生成论文 5.9.2 节结果表 |
| ├ `export_defense_data.py` | **导出答辩演示所需的全部真实数据**（`→ 文档/论文材料/_源文件/defense_data.json`） |
| └ `make_thesis_figures.py` | 生成论文插图（`results/figures_thesis/`） |
| **`tests/`** | 9 个测试文件、**63 个单元测试** |
| **`results/`** | **实验数据**（详见下表） |
| `configs/` | YAML 配置镜像（⚠️ **已过时且未被任何脚本读取**，见 K16，勿作配置依据） |
| `ldmt-experiment-report/` | `build_report.py` 生成的 HTML 报告（叙述文字为旧版，见 K15） |
| `analyze.py` | v2.4 冻结汇总与图表（写 `results/summary/`、`results/figures/`） |
| `build_report.py` / `serve_report.py` | HTML 报告生成与本地预览 |
| `main.py` | 一键流水线（pytest → analyze → build_report → 打开报告） |
| `conftest.py` | pytest 配置（注册 `smoke` marker） |
| `requirements.txt` | 依赖清单 |
| `README.md` | 本文件 |

### 2.3 `results/` 各子目录（哪些数据参与论文统计）

| 目录 | 内容 | 参与论文统计？ |
|---|---|---|
| `raw/` | v2.4 四组实验逐代 CSV（126 文件）+ `run_meta.jsonl` | ✅ 第 5.2–5.6、5.8–5.9.1 节 |
| `raw_v25/` | v2.5 两阶段实验逐代 CSV（72 文件） | ✅ 第 5.7、5.9.2 节 |
| `raw_v25_badauthoritative/` | 首版**错误** dMOP 定义下的 12 个 CSV | ❌ 仅施工留痕 |
| `raw_legacy_v21/` | v2.1 时期旧数据 | ❌ |
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
| ├ `LDT-NSGA-II毕业论文（硕士学位论文）.docx` / `.md` | 硕士论文正式稿（**177 页** / 51 表 / 11 图 / 1255 公式对象） |
| ├ `LDT-NSGA-II开题报告.docx` / `.md` | 开题报告（29 页 / 12 表） |
| ├ `LDT-NSGA-II V2.5答辩.html` | **当前答辩演示**（15 页 / 9 张动画图表） |
| ├ `LDT-NSGA-II V2.4答辩.html` | v2.4 答辩幻灯片（历史参考） |
| ├ `_源文件/` | **文档构建系统**：分章节 Markdown 源 + 构建脚本（改论文从改这里开始） |
| ├ `_源文件/_本科版/` | 本科版的独立入口与分章节源 |
| ├ `_audit/` | 独立审计与复核报告（代码审计、数值复核、文献核查、dMOP 定义考证） |
| ├ `_preview/` | 渲染出来的页面预览 PNG 与 PDF（QA 用，可删） |
| └ `_参考归档/` | 用户提供的早期参考文档与 v2.4 时期草稿、开发日志 |

**改论文的正确姿势**：编辑 `文档/论文材料/_源文件/chapterN_*.md`，然后运行
`python 文档/论文材料/_源文件/build_docx.py` 重新生成 DOCX 与交付用 Markdown。
本科版改 `_源文件/_本科版/*.md`，跑 `_源文件/_本科版/build_ug_docx.py`。

---

## 3. 快速上手：五分钟跑通

```powershell
cd LDT-NSGA-II

# 1) 安装依赖（实测 pymoo 0.6.1.5 + Python 3.9.13）
python -m pip install -r requirements.txt

# 2) 63 个单元测试（约 5 秒）
python -m pytest -q

# 3) 权威 dMOP 解析自检（应输出 ALL CHECKS PASSED）
python -m scripts.check_dmop_authoritative

# 4) 单次冒烟（分钟级；务必写入独立目录！）
python -m experiments.run_comparison --problems FDA3 --taus 25 --seeds 0 `
    --out results/_smoke --n-gen 50 --limit-seeds 1
```

> ⚠️ **两条必须知道的坑**
> 1. 冒烟数据用的是缩小预算（50 代而非 300 代），一旦写进 `results/raw` 就会**污染正式统计**。
> 2. `run_supplementary.py` 自身默认只跑 `DMOP1,DMOP2`（360 次）；复现论文的 540 次
>    必须显式加 `--problems DMOP1,DMOP2,DMOP3`（已知问题 K14）。

---

## 4. 算法与实验设计

### 4.1 LDT-NSGA-II 核心

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

两个**默认关闭**的可选机制：变化幅度驱动的自适应迁移宽度缩放（`clip(1+β|D_k|/mean|D_k|, 0.5, 2.0)`）
与多源历史二阶矩融合（历史 3、权重下限 0.05）。

### 4.2 六个算法（仅 `respond` 不同，其余全部相同）

| 标识 | 变化时的行为 |
|---|---|
| `nsga2_static` | 原样带入上一种群，无任何响应 |
| `dnsga2_inspired` | 保留种群前缀 80% + 20% 均匀随机移民 |
| `ar1_predictive` | 逐维 OLS 一阶自回归外推（历史 5）+ 20% 移民；**方差不预测** |
| `gmm_inspired` | 高斯混合（≤3 分量，`reg_covar=1e-4`）采样 + **10%** 移民 |
| `pps` | 中心全步差分 `2μ_{k−1}−μ_{k−2}` + 整体平移 + 20% 移民 |
| `ldt_nsga2` | 半步趋势外推 + 高斯采样迁移 + 50/20/30 重组（**本文方法**） |
| `ldt_alpha1` | 同上，但 `α = 1.0`（v2.5 机制对照变体） |

> 命名诚实性：`dnsga2_inspired` / `gmm_inspired` 是"受启发"的项目实现，**不是**原文复现；
> `pps` 是 PPS 的简化实现（无独立方差/流形预测项）。GMM 的移民比例是 10% 而非 20%（源码既有事实）。

### 4.3 九个测试问题

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

### 4.4 统计检验

- 两算法：配对 **Wilcoxon 符号秩检验**（种子完全对齐）+ **Holm** 逐步下降校正；
- 多算法：**Friedman** 检验 + **Nemenyi** 事后检验，`k=6, N=30, α=0.05` 时 `CD = 1.37668`；
- 主指标口径：**稳态 MIGD** —— 每个环境段最后一代 IGD 的段间平均，再对种子平均。
  另有"变化后 MIGD"（段首代）与"全代 MHV"两个辅助口径。

---

## 5. 实验结果（v2.5 最新）

### 5.1 规模与状态

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

### 5.2 核心结论

1. **步长归因（v2.5 主结论）**：α=1.0 在 12 个配置中的 11 个上更优，
   与 PPS 差距缩到 `-8.5%~+3.4%`；**2 个配置反超**（FDA4 τ=10 / τ=25，p 分别 1.86e-9 / 5.72e-7）。
   → 原"LDT 显著劣于 PPS"主要来自外推步长，而非迁移范式。
2. **反例与自适应规则**：唯一例外是 **DMOP1A（Type III，解集静止）**，
   该问题 α=0.5 领先 α=1.0 达 **127%**（0.0251 vs 0.0570）。
   → 外推步长应随"解集是否在显著移动"自适应。
3. **组件贡献**：分布迁移（去掉后退化 13–30%）≫ 趋势外推（12–19%）≫ 旧种群保留（无一致增益）。
4. **计算开销**：GMM 单次响应耗时约为同配置 LDT 的 **4.4–11 倍**且性能不占优；
   LDT 与 AR(1)、PPS 同处 10⁻² s 量级。
5. **v2.4 基线结论（仍成立）**：K1 修复后 PPS 在 9 个主比较组合均值第一、
   LDT（α=0.5）稳居第二；LDT 显著优于 static/dnsga2/gmm；在 Friedman–Nemenyi
   整体排序下与 PPS 于 4 个配置不可区分。

### 5.3 图表与汇总

- 论文插图：`results/figures_thesis/`（11 张）
- 汇总图表：`results/figures/`（箱线图、MIGD–τ 曲线、恢复曲线、消融图、敏感性图等）
- 论文级汇总：`results/summary_v25/`（`analysis_report.md` 可直接粘贴）
- v2.4 冻结汇总：`results/summary/`

---

## 6. 文档在哪儿

| 想要什么 | 去哪儿 |
|---|---|
| **本科毕业论文（28 页）** | [`../文档/论文材料/LDT-NSGA-II本科毕业论文.docx`](../文档/论文材料/) |
| 硕士学位论文详细版（177 页） | [`../文档/论文材料/LDT-NSGA-II毕业论文（硕士学位论文）.docx`](../文档/论文材料/) |
| 开题报告 DOCX | [`../文档/论文材料/LDT-NSGA-II开题报告.docx`](../文档/论文材料/) |
| **答辩演示（含动画图表）** | [`../文档/论文材料/LDT-NSGA-II V2.5答辩.html`](../文档/论文材料/) |
| 改论文内容 | `../文档/论文材料/_源文件/chapterN_*.md`（本科版在 `_本科版/`），然后跑对应构建脚本 |
| 重做答辩演示 | `python -m scripts.export_defense_data` → `python ../文档/论文材料/_源文件/make_defense.py` |
| 零基础运行手册 | 论文附录 A（与本文第 3 节一致） |
| 全部核心源码 | 论文附录 D（由构建脚本从本仓库原文件逐字嵌入，不会与代码脱节） |
| 独立审计/复核报告 | `../文档/论文材料/_audit/` |
| 开发过程日志 | `../文档/论文材料/_参考归档/`（3 份 docx + 实现方案 md） |

**重建文档**：

```powershell
# 硕士论文（177 页）
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

## 7. 已知问题清单

论文第 5.11 节以 **K1–K26** 完整披露。摘要如下（全部**不影响第 5 章结论**）：

| 编号 | 问题 | 状态 |
|---|---|---|
| K1 | AR1/GMM/PPS 的 `immigrant_ratio` 曾是死参数（移民被定容截断丢弃） | ✅ 已修复并重跑 1080 次 |
| K2 / K3 | DMOP3 实为静态；DMOP2 是项目变体非标准 dMOP2 | ✅ 已核实并如实标注 |
| K4–K9 | 坐标系混用、HV 回退 3 处 bug、无效字段、Friedman 并列秩、损坏 `main.py`、诊断列硬编码 | ✅ 已修复（0 次重跑） |
| K10 / K11 | FDA4 前沿实测 507 点；PlatEMO 无 DMOP*.m | ✅ 已核实 |
| K12 | HTML 报告非自包含 | ✅ 已修复 |
| K13 | 多源融合时趋势项量纲口径不一致 | ⚠️ 如实保留（默认未启用，不影响结论） |
| **K14** | `run_supplementary.py` 默认只跑 DMOP1,DMOP2（无法复现 540 次） | ⚠️ **待修** |
| **K15** | `build_report.py` 生成的 HTML 报告叙述文字是 v1 硬编码（表格数据正确） | ⚠️ **待修** |
| **K16** | `configs/*.yaml` 已过时且未被任何脚本读取 | ⚠️ **待清理** |
| **K17** | `README` 曾写 `sigma_min = 0.001*(u−l)`（真实尺度），代码是归一化 1e-3 | ⚠️ 本版已更正 |
| K18–K21 | 数值复核发现的 4 处不一致（含 GMM 移民实为 10%、"3780 次独立运行"措辞不准） | ✅ 已修正 |
| K22–K26 | 权威 dMOP1/2/3 定义错误、数据补齐、时变前沿判定、`get_true_pf` 重构 | ✅ 已修复并重跑 |

---

## 8. 复现实验的完整命令

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
python analyze.py                              # v2.4 冻结汇总 + 图表
python main.py                                 # 一键：测试 → 汇总 → HTML 报告
```

> **并行加速**：`run_v25_alpha.py` 可按 (问题, τ) 拆成 9 个独立进程并行；
> 在 32 逻辑核机器上把阶段一从约 5 小时压到约 0.5 小时。`run_v25_comparison.ps1`
> 是当时用的并行启动脚本（含日志目录与 PID 记录），可参考。

---

## 环境与版本说明

| 项 | 值 |
|---|---|
| Python | 3.9.13（实测） |
| pymoo | 0.6.1.5（⚠️ 0.6.x 无 `algo.n_eval`，真实计数器在 `algo.evaluator.n_eval`，见 `_algo_n_eval()`） |
| numpy / scipy / scikit-learn | 2.0.2 / 1.13.1 / 1.6.1 |
| **数据版本戳** | 磁盘 v2.4 数据的 `code_version` 标为 **`v2.2`**（由 commit `50adf19c5c5a` 生成），代码当前为 `v2.4`；v2.5 数据标 `v2.5`。该字段仅作标记，不影响数值（见论文 4.1 节 / K20） |
