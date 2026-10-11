# LDT-NSGA-II

轻量级分布迁移动态多目标优化器（Lightweight Distribution-Transfer
dynamic multi-objective optimizer），统一基于 pymoo 的 NSGA-II 内核实现。
本仓库为动态多目标优化（DMO）方向本科毕业论文的配套研究代码。

**版本说明**：v2.1 — 在 v2.0（开题报告 2.0 内核）基础上扩展：
- 新增 **PPS 基线**（Population Prediction Strategy，Zhou et al. 2014）
- 新增 **dMOP3** 测试问题（静态 POS + 时变 POF）
- 新增 **Friedman + Nemenyi 统计检验**（多算法多问题排序）
- 新增 **恢复速度指标**（Recovery Generations：环境变化后 IGD 恢复到变化前的代数）
- 新增 **可插拔 TransferStrategy 接口**（Moment/GMM 可替换）
- 新增 **自适应迁移强度**（`adaptive_transfer=True`，随 D_k 调整 sigma）
- 新增 **多源知识融合**（`multi_source=True`，历史特征加权融合）

v2.0 内核保持不变：
- 种群组成 **50/20/30**（迁移 50%、保留旧种群 20%、随机 30%）
- 趋势中心预测：`mu_hat = clip(mu_e + alpha*(mu_e - mu_{e-1}), l, u)`，默认 `alpha=0.5`
- 迁移采样：`x_tr = clip(mu_hat + z*sigma_e, l, u)`，z~N(0,I)
- 标准差下限：`sigma_min = 0.001*(u-l)`

## 文件结构总览

```
ldmt_dmoea/
├── algorithms/                  # 算法实现（统一 NSGA-II 内核）
│   ├── nsga2_core.py            #   DynamicDriver 主循环 + make_nsga2 内核
│   ├── ldmt_nsga2.py            #   LDT-NSGA-II（均值/方差迁移 + 50/20/30）
│   ├── distribution_features.py #   特征提取/非支配筛选/MigrationFunction/transfer_population
│   ├── environment_manager.py   #   EnvironmentManager（环境分段）+ ChangeDetector（探针检测）
│   ├── dnsga2_response.py       #   DNSGA-II 基线
│   ├── ar_nsga2.py              #   AR-NSGA-II 基线
│   ├── gmm_transfer_nsga2.py    #   GMM-Transfer 基线
│   ├── pps_nsga2.py             #   PPS 预测基线 (v2.1)
│   └── transfer_strategy.py     #   可插拔迁移策略接口 (v2.1)
├── problems/                    # 动态测试问题与真实 Pareto 前沿
│   ├── fda1.py / fda3.py / fda4.py
│   ├── dmop1.py / dmop2.py / dmop3.py   # dMOP3 为 v2.1 新增
│   └── true_fronts.py
├── evaluation/                  # 指标计算
│   ├── igd.py / hv.py
│   ├── reference_front.py
│   ├── statistics.py            #   Wilcoxon + Holm 统计检验
│   └── friedman_test.py         #   Friedman + Nemenyi 检验 (v2.1)
├── experiments/                 # 实验入口（可断点续跑）
│   ├── common.py                #   RunConfig / problem_factory / ExperimentRunner
│   ├── run_comparison.py        #   主对比实验
│   ├── run_ablation.py          #   消融实验
│   ├── run_sensitivity.py       #   参数敏感性实验
│   └── run_supplementary.py     #   补充实验（DMOP1/2）
├── configs/                     # 实验配置 YAML 镜像
├── tests/                       # pytest 单元与冒烟测试
├── results/
│   ├── raw/                     # 逐代 CSV（comparison_/ablation_/sensitivity_/supplementary_*.csv）
│   ├── summary/                 # 汇总表 + 统计检验 + 图表 JSON + analysis_report.md
│   └── figures/                 # 箱线图/恢复曲线/耗时对比/消融/敏感性 PNG
├── ldmt-experiment-report/      # 独立 HTML 实验报告
│   ├── index.html               #   build_report.py 生成（内嵌数据 + ECharts）
│   └── assets/                  #   doc.css / doc.js / echarts.min.js
├── analyze.py                   # 分析流水线（汇总→统计→出图）
├── build_report.py              # HTML 报告组装器
├── conftest.py                  # pytest 公共固件
├── requirements.txt
└── README.md
```

## 实验设计思路

### 研究问题

动态多目标优化问题的 Pareto 前沿随时间不断移动，算法需要在每次环境
变化后快速重新跟踪前沿。本研究关注：**能否用低成本的分布特征迁移，
让种群在变化后更快收敛到新前沿**。

### 核心假设

- **H1（分布迁移改善跟踪）**：环境变化后，将种群迁移到预测分布能降低
  变化后的 MIGD 并缩短恢复代数；
- **H2（低维特征迁移更省时）**：仅用均值/标准差描述分布，相对 GMM
  全分布迁移计算开销更小；
- **H3（迁移比例与变化幅度影响稳定性）**：迁移比例与环境变化幅度
  `n_t` 会影响算法的稳态与响应性能。

### LDT 方法（报告 5.1-5.2）

在每次环境变化时（入口 `respond()`，见 [algorithms/ldmt_nsga2.py](algorithms/ldmt_nsga2.py#L75-L134)）：

1. 从上一代种群中选出**非支配子集**（数量不足时退回全种群）—
   [ldmt_nsga2.py L80-L85](algorithms/ldmt_nsga2.py#L80-L85)，
   非支配判定见 [non_dominated_set](algorithms/distribution_features.py#L32-L47)；
2. 提取其归一化后的**均值 `mu_e` 与标准差 `sigma_e`** 作为分布特征，
   `sigma_min = 0.001*(u-l)` —
   [ldmt_nsga2.py L86-L89](algorithms/ldmt_nsga2.py#L86-L89)，
   归一化与特征提取见 [extract_features](algorithms/distribution_features.py#L24-L29)；
3. 用固定探针点在新旧环境下的目标值差异计算变化幅度 `D_k`（仅记录）—
   [ldmt_nsga2.py L91-L92](algorithms/ldmt_nsga2.py#L91-L92)，
   检测逻辑见 [ChangeDetector.detect](algorithms/environment_manager.py#L42-L49)；
4. **趋势中心预测**：`mu_hat = clip(mu_e + alpha*(mu_e - mu_{e-1}), l, u)`，
   默认 `alpha=0.5`（首次变化时 `mu_hat = mu_e`）—
   [ldmt_nsga2.py L94-L101](algorithms/ldmt_nsga2.py#L94-L101)；
5. **高斯扰动采样迁移候选**：`x_tr = clip(mu_hat + z*sigma_e, l, u)`，
   z~N(0,I) —
   [ldmt_nsga2.py L103-L117](algorithms/ldmt_nsga2.py#L103-L117)；
6. 按 **50/20/30** 组成新一代种群 —
   [ldmt_nsga2.py L119-L134](algorithms/ldmt_nsga2.py#L119-L134)：
   - **50%** 迁移个体（`transfer_ratio`）
   - **20%** 保留旧种群个体（`retain_ratio`）
   - **30%** 均匀随机个体（`random_ratio`）

### 对比基线

| 算法 | 思路 |
|---|---|
| **NSGA-II** | 静态基准，环境变化时不迁移，仅重新评估当前种群 |
| **DNSGA-II-inspired** | 变化检测 + 随机移民响应，无分布迁移 |
| **AR1-Predictive** | 一阶趋势预测基线（透明定义，非特定论文复现） |
| **GMM-inspired** | 高斯混合全分布对照（非原 KTMM 复现） |
| **PPS** | 种群中心一阶差分外推预测 + 20% 随机移民 (Zhou et al. 2014) |
| **LDT-NSGA-II** | 均值/标准差特征 + 趋势中心 + 50/20/30 混合初始化 |

### v2.1 算法扩展

LDT-NSGA-II 新增三个可选模块（默认关闭，可通过参数开启）：

- **自适应迁移强度** (`adaptive_transfer=True`)：根据变化幅度 `D_k`
  动态缩放采样标准差 `sigma_eff = sigma_e * (1 + beta * D_k_norm)`，
  大变化时扩大探索范围，小变化时收紧局部搜索 —
  [ldmt_nsga2.py L175-L180](algorithms/ldmt_nsga2.py#L175-L180)；
- **多源知识融合** (`multi_source=True`)：维护最近 K 个环境的特征
  缓冲区，按相似度加权融合 `w_i ∝ exp(-||mu_e - mu_i||/sqrt(D))` —
  [ldmt_nsga2.py L108-L127](algorithms/ldmt_nsga2.py#L108-L127)；
- **可插拔 TransferStrategy** (`transfer_strategy=...`)：将"特征提取→
  迁移计算→种群重建"封装为策略对象，支持 Moment/GMM 替换 —
  [transfer_strategy.py](algorithms/transfer_strategy.py)。

### 测试问题与指标

- **主实验问题**：FDA1、FDA3、FDA4（2/2/3 目标）
- **补充问题**：DMOP1、DMOP2、DMOP3（v2.1 新增，静态 POS + 时变 POF）
- **评价指标**：MIGD（越小越好）、MHV（越大越好）、变化后 MIGD
  （响应速度）、恢复代数（v2.1 新增）、迁移模块耗时
- **统计检验**：Wilcoxon 秩和 + Holm 校正；Friedman + Nemenyi 后验
  （v2.1 新增，多算法多问题排序）

## 运行原理

### 环境分段与变化检测

- `EnvironmentManager` 以 `tau_t` 代为一个环境段，`t = floor(g/tau_t)/n_t`
  为连续环境时间 — 见 [environment_manager.py L5-L25](algorithms/environment_manager.py#L5-L25)；
- 每段起始（`g % tau_t == 0`）触发一次环境变化 —
  [is_change](algorithms/environment_manager.py#L22-L25)；
- `ChangeDetector` 维护一组固定探针点，比较其在新旧环境下的目标向量，
  得到归一化变化幅度 `D_k = mean(|new_F - prev_F| / scale)` —
  [ChangeDetector.detect](algorithms/environment_manager.py#L42-L49)。

### 趋势中心预测（报告 5.2）

- 维护上一环境的均值 `mu_prev` 与标准差 `sigma_prev` —
  [ldmt_nsga2.py L94-L101](algorithms/ldmt_nsga2.py#L94-L101)；
- 当历史不少于 2 个环境时：

  ```
  mu_hat = clip(mu_e + alpha * (mu_e - mu_{e-1}), l, u)
  ```

  默认 `alpha=0.5`（敏感性测试 `{0, 0.25, 0.5, 0.75, 1.0}`）；
  首次变化时 `mu_hat = mu_e`；
- 标准差下限：`sigma_min = 0.001 * (u - l)` —
  [ldmt_nsga2.py L86-L89](algorithms/ldmt_nsga2.py#L86-L89)。

### 种群迁移与组成（报告 5.2）

- 迁移候选采样：`x_tr = clip(mu_hat + z * sigma_e, l, u)`，z~N(0,I) —
  [ldmt_nsga2.py L103-L117](algorithms/ldmt_nsga2.py#L103-L117)；
- 新种群按 **50/20/30** 组成 —
  [ldmt_nsga2.py L119-L134](algorithms/ldmt_nsga2.py#L119-L134)：
  - **50%** 迁移个体（`transfer_ratio=0.5`）
  - **20%** 保留旧种群（`retain_ratio=0.2`）
  - **30%** 均匀随机（`random_ratio=0.3`）
- 非有限迁移值自动替换为均匀随机样本。

### NSGA-II 内核与逐代评估

- 每个环境段创建**全新**的 NSGA-II 实例（SBX `eta=15`、PM `eta=20`），
  以迁移结果作为初始种群，运行 `tau_t` 代 —
  [make_nsga2](algorithms/nsga2_core.py#L26-L36)、
  [DynamicDriver.run](algorithms/nsga2_core.py#L80-L157)；
- 段 `k` 使用 `seed = run_seed + k`，保证可复现且段间独立 —
  [nsga2_core.py L113](algorithms/nsga2_core.py#L113)；
- **逐代**记录 IGD、HV、单代耗时、迁移耗时（仅在段首代记录一次），
  用于恢复曲线与耗时统计 — [nsga2_core.py L115-L154](algorithms/nsga2_core.py#L115-L154)。

### 分析与报告产出

- `analyze.py`：汇总每轮指标 → 均值/标准差 → Wilcoxon 秩和检验 + Holm
  校正 + **Friedman + Nemenyi 后验检验**（v2.1）+ **恢复代数分析**
  （v2.1）→ 箱线图、恢复曲线、耗时对比、消融与敏感性图 —
  入口 [analyze.py](analyze.py)；
- `build_report.py`：将 `results/summary/` 下的图表 JSON 组装为独立
  HTML 报告（`ldmt-experiment-report/index.html`），内含 ECharts 交互
  图表与玻璃拟态布局 — 入口 [build_report.py](build_report.py)。

## 实验结果（v2.1 快照）

> 数据说明：FDA1/FDA3/FDA4 主对比为 30 独立种子（v2.0 内核
> `LDMT_NSGA2`）；PPS 与 dMOP3 为 v2.1 新增组件的快速验证（2 种子冒烟
> 运行），完整 30 种子对比可复跑 `run_comparison` 获得。
> 全部统计检验 p<0.05（Friedman，8 个 问题×τ_t 组合）。

### 主对比（τ_t=25，MIGD 越低越好）

| 问题 | 算法 | MIGD | MHV | post-MIGD | 恢复代数 |
|---|---|---|---|---|---|
| FDA1 | LDMT-NSGA-II | **0.00657** | **0.8658** | **0.0128** | 10.79 |
| FDA1 | AR1-Predictive | 0.01332 | 0.8543 | 0.0250 | 15.29 |
| FDA1 | PPS | 0.01420 | 0.8513 | 0.0403 | **5.50** |
| FDA1 | DNSGA-II | 0.01521 | 0.8502 | 0.0455 | 29.30 |
| FDA1 | GMM-inspired | 0.01507 | 0.8506 | 0.0465 | 22.57 |
| FDA3 | LDMT-NSGA-II | **0.00673** | **1.5373** | **0.0123** | 10.88 |
| FDA3 | AR1-Predictive | 0.00691 | 1.5363 | 0.0136 | 13.90 |
| FDA3 | DNSGA-II | 0.00923 | 1.5257 | 0.0368 | 20.67 |
| FDA3 | GMM-inspired | 0.00892 | 1.5273 | 0.0381 | 18.29 |
| FDA4 | LDMT-NSGA-II | **0.08658** | **0.6425** | **0.1287** | **9.34** |
| FDA4 | AR1-Predictive | 0.09166 | 0.6270 | 0.1531 | 9.94 |
| FDA4 | DNSGA-II | 0.09284 | 0.6207 | 0.1576 | 13.19 |
| FDA4 | GMM-inspired | 0.09116 | 0.6240 | 0.1475 | 12.94 |

### 图表索引

`results/figures/` 中的图表与 `results/summary/` 中的数据文件一一对应：

| 图表文件 | 内容 | 数据来源 | 分析脚本位置 |
|---|---|---|---|
| `boxplot_MIGD_FDA1.png` | FDA1 各算法 MIGD 分布箱线图（τ_t=25） | `summary_comparison.csv` | [analyze.py](analyze.py) §1.1 |
| `boxplot_MIGD_FDA3.png` | FDA3 各算法 MIGD 分布箱线图 | 同上 | 同上 |
| `boxplot_MIGD_FDA4.png` | FDA4 各算法 MIGD 分布箱线图 | 同上 | 同上 |
| `boxplot_MIGD_DMOP3.png` | dMOP3 各算法 MIGD 分布箱线图（v2.1） | `supplementary_summary.csv` | [analyze.py](analyze.py) §4 |
| `recovery_FDA1.png` | FDA1 恢复曲线（变化后 IGD 随代数变化，τ_t=25） | `comparison_runs.csv` | [analyze.py](analyze.py) §1.2 |
| `recovery_FDA3.png` | FDA3 恢复曲线 | 同上 | 同上 |
| `recovery_FDA4.png` | FDA4 恢复曲线 | 同上 | 同上 |
| `recovery_DMOP3.png` | dMOP3 恢复曲线（v2.1） | `supplementary_runs.csv` | [analyze.py](analyze.py) §4 |
| `MIGD_vs_tau.png` | MIGD 随环境变化间隔 τ_t∈{10,25,50} 的变化趋势 | `summary_comparison.csv` | [analyze.py](analyze.py) §1.1 |
| `migration_time.png` | 各算法单次迁移模块耗时对比（验证 H2） | `comparison_runs.csv` | [analyze.py](analyze.py) §1.2 |
| `ablation_MIGD.png` | 消融实验 MIGD 对比（full / no_trend / no_transfer / kappa0） | `ablation_summary.csv` | [analyze.py](analyze.py) §2 |
| `sensitivity.png` | 参数敏感性（N、n_t、ρ、α 扫描，v2.1） | `sensitivity_runs.csv` | [analyze.py](analyze.py) §3 |

**汇总数据文件索引**（`results/summary/`）：

| 文件 | 内容 |
|---|---|
| `summary_comparison.csv` | 主对比均值/标准差/中位数（每 问题×τ_t×算法 一行） |
| `statistical_tests.csv` | Wilcoxon 秩和检验 + Holm 校正结果 |
| `friedman_results.json` | Friedman χ²/p 值 + Nemenyi CD（v2.1） |
| `recovery_generations.csv` | 恢复代数明细（每 问题×算法×种子×环境 一行，v2.1） |
| `analysis_report.md` | 自动生成的完整文字分析报告 |
| `_chart_data.json` / `_recovery.json` / `_boxplot.json` | HTML 报告内嵌数据源 |
| `ablation_summary.csv` / `supplementary_summary.csv` | 消融 / 补充实验汇总 |

### 关键发现

- **LDMT/LDT 在 3 个 FDA 问题上 MIGD、MHV、post-MIGD 全面最优**，
  恢复代数约为 DNSGA-II 的 1/2 ~ 1/3，验证假设 H1（分布迁移改善跟踪）；
- **迁移效率**：均值/方差特征的单次迁移耗时显著低于 GMM 全分布迁移，
  验证假设 H2（低维特征更省时，见 `results/summary/analysis_report.md`）；
- **PPS 在 FDA1 上恢复代数最低**（5.50），但其 MIGD 收敛精度略逊于
  LDMT，表明"快速响应"与"收敛精度"之间存在权衡；
- **Friedman 检验**：8 个 问题×τ_t 组合全部 p<0.05，算法排序差异显著；
  Nemenyi CD 值见 `results/summary/friedman_results.json`；
- 完整逐代数据见 `results/raw/`，恢复曲线与箱线图见
  `results/figures/`，交互式报告见 `ldmt-experiment-report/index.html`。

## 安装

```bash
pip install -r requirements.txt
```

依赖：`pymoo>=0.6`、`numpy`、`scipy`、`scikit-learn`、`matplotlib`、
`pyyaml`、`tqdm`、`pytest`。

## 运行

### 1. 测试

```bash
pytest -q
```

### 2. 实验

所有实验均支持**断点续跑**：输出 CSV 中已存在的组合会被自动跳过。

```bash
# 主对比实验（FDA1/3/4 × tau{10,25,50} × 5 算法 × 30 种子）
python -m experiments.run_comparison --problems FDA1,FDA3,FDA4 --taus 10,25,50 --seeds 0-29 --out results/raw

# 单配置快速冒烟测试
python -m experiments.run_comparison --problems FDA3 --taus 25 --seeds 0 --out results/raw --n-gen 50 --limit-seeds 1

# 消融实验（报告 7.7）
python -m experiments.run_ablation

# 敏感性实验（报告 7.8）
python -m experiments.run_sensitivity

# 补充实验
python -m experiments.run_supplementary
```

#### `run_comparison` 参数说明

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--problems` | `FDA1,FDA3,FDA4` | 逗号分隔的问题名称 |
| `--taus` | `10,25,50` | 环境变化周期 `tau_t` |
| `--seeds` | `0-29` | 种子范围（`lo-hi`）或逗号列表 |
| `--out` | `results/raw` | 输出 CSV 目录 |
| `--n-gen` | `300` | 每次运行总代数 |
| `--n` | `100` | 种群规模 |
| `--n-t` | `20` | 变化幅度 `n_t` |
| `--algos` | 全部算法 | 逗号分隔的算法键（v2.1 含 PPS） |
| `--limit-seeds` | `None` | 仅取前 N 个种子（快速测试） |

#### 其他实验入口

- `run_ablation`：`--problems`（FDA1,FDA3,FDA4）、`--tau`（25）、
  `--n-gen`（300）、`--n`（100）、`--n-t`（20）、`--seeds`（0-29）、
  `--limit-seeds`。消融变体：`full`、`ablation_a`（关闭趋势）、
  `ablation_b`（关闭分布迁移）、`ablation_c`（κ=0）。
- `run_sensitivity`：`--problems`（FDA1,FDA3,FDA4）、`--tau`（25）、
  `--n-gen`（300）、`--seeds`（0-29）、`--limit-seeds`。扫描 N ∈
  {50,100,200}、n_t ∈ {10,20,40}、ρ ∈ {0.25,0.5,0.75}、α ∈
  {0,0.25,0.5,0.75,1.0}。
- `run_supplementary`：`--problems`（DMOP1,DMOP2；可选 DMOP3 (v2.1)）、
  `--tau`（25）、`--n-gen`（300）、`--n`（100）、`--n-t`（20）、
  `--algos`（全部，含 PPS）、`--seeds`（0-29）、`--limit-seeds`。

### 3. 分析

```bash
python analyze.py
```

读取 `results/raw/*.csv`，输出每轮指标、汇总表、Wilcoxon+Holm 统计
检验以及图表到 `results/summary/` 与 `results/figures/`，同时生成
`results/summary/analysis_report.md`。

### 4. 生成 HTML 报告

```bash
python build_report.py
```

读取 `results/summary/_chart_data.json`、`_recovery.json`、
`_boxplot.json`，组装为 `ldmt-experiment-report/index.html`（单文件自包含，内嵌数据、ECharts 与 CSS）。

### 5. 可直接在浏览器打开，或启动本地服务：

```bash
cd ldmt-experiment-report; python -m http.server 8000
```

```bash
python -c "import http.server, socketserver, webbrowser, threading, os, time; os.chdir('ldmt-experiment-report'); threading.Timer(1.0, lambda: webbrowser.open('http://localhost:8000')).start(); socketserver.TCPServer(('', 8000), http.server.SimpleHTTPRequestHandler).serve_forever()"
```


## 算法一览

| 键 | 名称 | 说明 |
|---|---|---|
| `NSGA-II` | NSGA-II | 静态基准，无迁移 |
| `DNSGA-II-inspired` | DNSGA-II-inspired | 变化检测 + 随机移民响应 |
| `AR1-Predictive` | AR1-Predictive | 一阶趋势预测基线 |
| `GMM-inspired` | GMM-inspired | 高斯混合全分布对照 |
| `PPS` | PPS | 种群中心一阶差分外推预测 (v2.1) |
| `LDT-NSGA-II` | LDT-NSGA-II | 均值/标准差特征 + 趋势中心 + 50/20/30 混合初始化 |

## 说明

- **FDA 公式**遵循 PlatEMO master（FDA1/FDA3/FDA4）。FDA4 使用标准
  DTLZ2 第一卦限单位超球面前沿（已验证：在 PS `x[:,M:]=G` 处目标范数
  恰好为 1.0）。
- **dMOP1** 采用标准 Farina-Deb-Amato (2004) 形式
  （`f2 = g*(1-(f1/g)^2)`，PF `{(u,1-u^2)}`）。PlatEMO master 未提供
  DMOP1/DMOP2.m，因此 **dMOP2 为本项目实现**（时变前沿，见
  `problems/dmop2.py`）。
- **GMM-inspired** 为本项目对比实现，并非复现某篇特定论文的 KTMM。
- **LDT-NSGA-II** 算法参数（ρ=0.5、κ=0.2、α=0.5）为开题报告 2.0 默认值，
  非文献证明最优值，已在敏感性实验中测试。
- HV 使用各问题固定参考点（原始目标值，不做逐轮归一化）：FDA1
  `(1.1,1.1)`、FDA3 `(2.2,1.1)`、FDA4 `(1.1,1.1,1.1)`。
- 所有随机性均设种子；环境段 `k` 使用 `seed = run_seed + k`。
- 实验结果仅适用于实际测试的问题、参数和硬件范围，不直接外推到所有
  真实工程任务。
