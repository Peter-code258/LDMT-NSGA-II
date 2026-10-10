# LDMT-NSGA-II

轻量级分布迁移动态多目标优化器（Lightweight Distribution-Migration
dynamic multi-objective optimizer），统一基于 pymoo 的 NSGA-II 内核实现。
本仓库为动态多目标优化（DMO）方向本科毕业论文的配套研究代码。

## 文件结构总览

```
ldmt_dmoea/
├── algorithms/                  # 算法实现（统一 NSGA-II 内核）
│   ├── nsga2_core.py            #   DynamicDriver 主循环 + make_nsga2 内核
│   ├── ldmt_nsga2.py            #   LDMT-NSGA-II（均值/方差迁移 + 70/20/10）
│   ├── distribution_features.py #   特征提取/非支配筛选/MigrationFunction/transfer_population
│   ├── environment_manager.py   #   EnvironmentManager（环境分段）+ ChangeDetector（探针检测）
│   ├── dnsga2_response.py       #   DNSGA-II 基线
│   ├── ar_nsga2.py              #   AR-NSGA-II 基线
│   └── gmm_transfer_nsga2.py    #   GMM-Transfer 基线
├── problems/                    # 动态测试问题与真实 Pareto 前沿
│   ├── fda1.py / fda3.py / fda4.py
│   ├── dmop1.py / dmop2.py
│   └── true_fronts.py
├── evaluation/                  # 指标计算
│   ├── igd.py / hv.py
│   ├── reference_front.py
│   └── statistics.py            #   Wilcoxon + Holm 统计检验
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

### LDMT 方法

在每次环境变化时（入口 `respond()`，见 [algorithms/ldmt_nsga2.py](algorithms/ldmt_nsga2.py#L68-L134)）：

1. 从上一代种群中选出**非支配子集**（数量不足时退回全种群）—
   [ldmt_nsga2.py L74-L79](algorithms/ldmt_nsga2.py#L74-L79)，
   非支配判定见 [non_dominated_set](algorithms/distribution_features.py#L32-L47)；
2. 提取其归一化后的**均值 `mu` 与标准差 `sigma`** 作为分布特征 —
   [ldmt_nsga2.py L82](algorithms/ldmt_nsga2.py#L82)，
   归一化与特征提取见 [extract_features](algorithms/distribution_features.py#L24-L29)；
3. 用固定探针点在新旧环境下的目标值差异计算变化幅度 `D_k` —
   [ldmt_nsga2.py L85-L86](algorithms/ldmt_nsga2.py#L85-L86)，
   检测逻辑见 [ChangeDetector.detect](algorithms/environment_manager.py#L42-L49)；
4. 用 AR(1) 风格外推预测下一环境的分布 `(mu_hat, sigma_hat)`，步长
   `lambda_k` 根据 `D_k` 的变化自适应裁剪到 `[0.5, 1.5]` —
   [ldmt_nsga2.py L88-L100](algorithms/ldmt_nsga2.py#L88-L100)，
   预测公式见 [MigrationFunction.predict](algorithms/distribution_features.py#L76-L87)；
5. 将上一代整体种群按 `(mu_prev, sigma_prev) -> (mu_hat, sigma_hat)`
   重新中心化与缩放，得到迁移种群 —
   [ldmt_nsga2.py L103-L107](algorithms/ldmt_nsga2.py#L103-L107)，
   迁移实现见 [transfer_population](algorithms/distribution_features.py#L90-L100)；
6. 按 **70/20/10** 组成新一代种群 —
   [ldmt_nsga2.py L109-L127](algorithms/ldmt_nsga2.py#L109-L127)：
   - **70%** 迁移个体
   - **20%** 局部变异个体（对迁移个体的随机维度加 `N(0, 0.05)` 噪声）
   - **10%** 均匀随机移民

### 对比基线

| 算法 | 思路 |
|---|---|
| **DNSGA-II** | 变化检测 + 随机移民响应，无分布迁移 |
| **AR-NSGA-II** | 仅用 AR(1) 预测均值，不做分布缩放迁移 |
| **GMM-Transfer** | 用高斯混合模型完整描述并迁移分布 |
| **LDMT-NSGA-II** | 均值/标准差特征 + AR(1) 迁移 + 70/20/10 组成 |

### 测试问题与指标

- **主实验问题**：FDA1、FDA3、FDA4（2/2/3 目标）
- **补充问题**：DMOP1、DMOP2
- **评价指标**：MIGD（越小越好）、MHV（越大越好）、变化后 MIGD
  （响应速度）、恢复曲线（IGD 随环境内代数变化）、迁移模块耗时

## 运行原理

### 环境分段与变化检测

- `EnvironmentManager` 以 `tau_t` 代为一个环境段，`t = floor(g/tau_t)/n_t`
  为连续环境时间 — 见 [environment_manager.py L5-L25](algorithms/environment_manager.py#L5-L25)；
- 每段起始（`g % tau_t == 0`）触发一次环境变化 —
  [is_change](algorithms/environment_manager.py#L22-L25)；
- `ChangeDetector` 维护一组固定探针点，比较其在新旧环境下的目标向量，
  得到归一化变化幅度 `D_k = mean(|new_F - prev_F| / scale)` —
  [ChangeDetector.detect](algorithms/environment_manager.py#L42-L49)。

### 分布迁移预测（MigrationFunction）

- 维护 `(mu_k, sigma_k, D_k)` 历史 —
  [MigrationFunction](algorithms/distribution_features.py#L50-L87)；
- 当历史不少于 2 个变化幅度时：

  ```
  lambda_k = clip(1.0 + alpha * (D_k - D_prev) / (D_prev + 1e-9), 0.5, 1.5)
  ```

  实现见 [distribution_features.py L78-L82](algorithms/distribution_features.py#L78-L82)；
  首次变化时 `lambda_k = 1.0`；
- 预测下一分布：

  ```
  mu_hat    = mu_prev  + lambda_k * (mu_prev  - mu_prev2)
  sigma_hat = max(sigma_min, sigma_prev + lambda_k * (sigma_prev - sigma_prev2))
  ```

  实现见 [distribution_features.py L83-L86](algorithms/distribution_features.py#L83-L86)。

### 种群迁移与组成

- `transfer_population` 将上一代种群 `Z` 按比例
  `sigma_hat / sigma_prev` 缩放并平移到 `mu_hat`，再裁剪回 `[0,1]` —
  [transfer_population](algorithms/distribution_features.py#L90-L100)；
- 迁移种群经 70/20/10 拆分后，局部变异部分对随机选取的 `D//3` 个维度
  加 `N(0, 0.05)` 噪声，移民部分在 `[xl, xu]` 内均匀采样 —
  [ldmt_nsga2.py L109-L127](algorithms/ldmt_nsga2.py#L109-L127)。

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
  校正 → 箱线图、恢复曲线、耗时对比、消融与敏感性图 —
  入口 [analyze.py](analyze.py)；
- `build_report.py`：将 `results/summary/` 下的图表 JSON 组装为独立
  HTML 报告（`ldmt-experiment-report/index.html`），内含 ECharts 交互
  图表与玻璃拟态布局 — 入口 [build_report.py](build_report.py)。

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
# 主对比实验（FDA1/3/4 × tau{10,25,50} × 4 算法 × 30 种子）
python -m experiments.run_comparison --problems FDA1,FDA3,FDA4 \
    --taus 10,25,50 --seeds 0-29 --out results/raw

# 单配置快速冒烟测试
python -m experiments.run_comparison --problems FDA3 --taus 25 \
    --seeds 0 --out results/raw --n-gen 50 --limit-seeds 1

# 消融 / 敏感性 / 补充实验
python -m experiments.run_ablation
python -m experiments.run_sensitivity
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
| `--algos` | 全部 4 种 | 逗号分隔的算法键 |
| `--limit-seeds` | `None` | 仅取前 N 个种子（快速测试） |

#### 其他实验入口

- `run_ablation`：`--problems`（FDA1,FDA3,FDA4）、`--tau`（25）、
  `--n-gen`（300）、`--n`（100）、`--n-t`（20）、`--seeds`（0-29）、
  `--limit-seeds`。消融变体：`full`、`no_feature_simplification`、
  `no_transfer`。
- `run_sensitivity`：`--problems`（FDA3）、`--tau`（25）、`--n-gen`
  （300）、`--seeds`（0-29）、`--limit-seeds`。扫描 N ∈ {50,100,150}
  与 n_t ∈ {10,20,40}。
- `run_supplementary`：`--problems`（DMOP1,DMOP2）、`--tau`（25）、
  `--n-gen`（300）、`--n`（100）、`--n-t`（20）、`--algos`（全部）、
  `--seeds`（0-29）、`--limit-seeds`。

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
`_boxplot.json`，组装为 `ldmt-experiment-report/index.html`（单文件
自包含，内嵌数据、ECharts 与 CSS）。可直接在浏览器打开，或启动
本地服务：

```bash
cd ldmt-experiment-report && python -m http.server 8000
```

## 算法一览

| 键 | 名称 | 说明 |
|---|---|---|
| `DNSGA2` | DNSGA-II | 变化检测 + 随机移民响应 |
| `AR_NSGA2` | AR-NSGA-II | AR(1) 均值预测，不做分布迁移 |
| `GMM_Transfer` | GMM-Transfer | 完整高斯混合分布迁移 |
| `LDMT_NSGA2` | LDMT-NSGA-II | 均值/方差特征 + AR(1) 迁移 + 70/20/10 组成 |

## 说明

- **FDA 公式**遵循 PlatEMO master（FDA1/FDA3/FDA4）。FDA4 使用标准
  DTLZ2 第一卦限单位超球面前沿（已验证：在 PS `x[:,M:]=G` 处目标范数
  恰好为 1.0）。
- **dMOP1** 采用标准 Farina-Deb-Amato (2004) 形式
  （`f2 = g*(1-(f1/g)^2)`，PF `{(u,1-u^2)}`）。PlatEMO master 未提供
  DMOP1/DMOP2.m，因此 **dMOP2 为本项目实现**（时变前沿，见
  `problems/dmop2.py`）。
- **GMM-Transfer** 为本项目对比实现，并非复现某篇特定论文的 KTMM。
- HV 使用各问题固定参考点（原始目标值，不做逐轮归一化）：FDA1
  `(1.1,1.1)`、FDA3 `(2.2,1.1)`、FDA4 `(1.1,1.1,1.1)`。
- 所有随机性均设种子；环境段 `k` 使用 `seed = run_seed + k`。
