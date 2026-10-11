# LDT-NSGA-II

轻量级分布迁移动态多目标优化器（Lightweight Distribution-Transfer
dynamic multi-objective optimizer），统一基于 pymoo 的 NSGA-II 内核实现。
本仓库为动态多目标优化（DMO）方向本科毕业论文的配套研究代码。

**版本说明：v2.4 — 已知问题清单（K2–K12）清理版（当前版本）**

相对 v2.3 的修正：**全部为代码质量与文档如实性修复，不改变任何实验
行为，0 次实验重跑**（分析流水线 analyze.py 与 build_report.py 已重跑
验证，Friedman 临界差 CD=1.37668 逐位不变）：

- **K2（DMOP3 身份更正）**：确认当前 DMOP3 为项目自定义**静态**欺骗
  问题（PS/PF 均不随时间变化），删除从未使用的死代码
  `G = |sin(0.5πt)|`；docstring 如实更正——dMOP 套件权威出处为
  **Goh & Tan (2009) IEEE TEVC 13(1):103–127**（非旧版误引的
  Farina 2004），且权威 dMOP3 是 Type I（PS 时变、PF 静态），本实现
  仅作"静态对照组"保留（Helbig 2012 博士论文 p.52–54 式 3.40）；
- **K3（DMOP2 出处标注）**：docstring 补充权威 dMOP2 定义
  （`g=1+9Σ(x_i−G)²`、`f₂=g(1−(f₁/g)^{H(t)})`、`H(t)=0.75sin(0.5πt)+1.25`，
  Type II）与本实现的差异说明，明确"项目变体、非标准 dMOP2"；
- **DMOP1 归属更正**：旧 docstring 误称"authoritative Farina-Deb-Amato
  (2004)"，实际为项目变体（系数 1 而非 9、平方前沿而非 sqrt 前沿），
  归属更正为 Goh & Tan (2009)；
- **K4（MomentTransferStrategy 坐标系）**：可插拔策略类旧实现混用
  归一化/真实坐标（真实尺度 sigma 下限 + 真实界 clip + 二次归一化）；
  现全程归一化空间、单次反归一化，与 LDTNSGA2 内联实现一致
  （该策略类非默认路径，不影响实验）；
- **K5（HV 精确回退三处 bug）**：`_hv_exact` 容斥实现修复——①奇数
  子集符号反转、②交集应取分量 max（误用 min）、③循环多算空集项；
  新增对拍测试验证与 pymoo HV 逐位一致（pymoo 在场时回退本就不执行）；
- **K6（无效字段）**：删除 `RunConfig.evals_per_env`（真实预算为
  N×τ_t×环境数），`code_version` 默认升至 v2.4；
- **K7（Friedman 并列秩）**：`_rank_per_seed` 改用
  `scipy.stats.rankdata` 平均秩（旧 argsort 法对并列给任意不同秩）；
  Nemenyi 的 `alpha` 参数现在真正选择 q 表（0.10/0.05 双档，Demšar
  2006 Table 5），不支持的 alpha 显式抛错；连续 MIGD 数据无并列，
  统计结果不变（CD=1.37668 已复核）；
- **K8（损坏文件）**：删除根目录 `main.py`（含未解决合并冲突标记的
  IDE 模板）；v2.4 起根目录 `main.py` 重建为**一键流水线入口**
  （测试 → analyze → build_report → 自动打开报告，见文件结构节）；
- **K9（诊断列真实性）**：static/DNSGA2/AR1/PPS 的
  `boundary_clip_ratio` 由硬编码 0.0 改为真实计算（该列仅作诊断，
  analyze.py 不消费）；顺带修复 ARNSGA2 接收但未存储 `n_probes` 的
  潜在 AttributeError；
- **K10（FDA4 参考前沿口径）**：实测确认角度网格 23×23=529 去重后为
  **507 点**，文档与测试统一该数字；
- **K11（PlatEMO 核查）**：经 GitHub git-tree API 递归核查 PlatEMO
  master（commit d25e65d）确认无 DMOP1/2/3.m，文档如实标注；
- **K12（报告自包含）**：`build_report.py` 改为构建时内联
  `assets/doc.css` 与 `assets/echarts.min.js`，生成的
  `ldmt-experiment-report/index.html` 现为**单文件自包含**（1.15 MB），
  可直接分发。

新增 `tests/test_v24_fixes.py` 15 个回归测试，测试总数 **63 个，
全部通过**。

---

<details>
<summary><b>v2.3 — K1 修复重跑版（历史版本）</b></summary>

四组实验共 **3780 次独立运行全部完成，0 失败**，输出 126 个逐代 CSV、
14 份汇总文件与 12 张 PNG 图表。相对 v2.2 的主要修正：

- **修复 K1（死参数）**：PPS/AR1/GMM 三个基线的 `immigrant_ratio` 原
  实现中追加的随机移民被定容截断全部丢弃（实际为"纯预测、无移民"）；
  现改为**先预留移民名额再拼接**（`algorithms/ar_nsga2.py`、
  `gmm_transfer_nsga2.py`、`pps_nsga2.py`），新增 3 个回归测试固化
  （测试总数 48 个）；v2.3 旧数据归档于
  `results/raw_legacy_v22_k1/`，不参与 v2.3 统计；
- **重跑受影响实验**：主比较与补充实验中三基线的全部 1080 次运行
  （36 个 CSV）在 K1 修复后重跑，LDT/static/DNSGA2 及消融/敏感性
  数据不受影响、沿用原结果；
- **结论更新（如实）**：K1 修复后 20% 随机移民注入使 PPS 的中心
  外推获得多样性补充，**PPS 在全部 9 个主比较组合与 DMOP1/DMOP2 上
  稳态 MIGD 第一**，LDT 退居第二——本仓库所有文档按新数据如实更新。

</details>

相对 v2.1 的主要修正（v2.2 引入）：

- 特征提取/迁移全程统一在归一化坐标系内完成，修复旧版归一化坐标混用问题；
- 修复 PPS 基线对已归一化均值二次归一化的坐标 bug；
- 修正 FDA4 实现逻辑（标准 DTLZ2 第一卦限单位超球面前沿）；
- 分离静态 NSGA-II 与 DNSGA-II 响应基线的行为；
- 修复统计检验流水线的运行时错误；
- 为采样标准差缩放增加上下限约束；
- 修复新版 pymoo 的 API 兼容性问题；
- 新增 6 个专项测试，测试总数 45 个，全部通过；v2.1 旧实验数据归档于
  `results/raw_legacy_v21/`，不参与 v2.2 统计。

**v2.1（历史版本）** 在 v2.0（开题报告 2.0 内核）基础上扩展：

- 新增 **PPS 基线**（Population Prediction Strategy，Zhou et al. 2014）
- 新增 **DMOP3** 测试问题
- 新增 **Friedman + Nemenyi 统计检验**（多算法多问题排序）
- 新增 **恢复速度指标**（Recovery Generations：环境变化后 IGD 恢复到变化前水平的代数）
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
LDT-NSGA-II/
├── algorithms/                  # 算法实现（统一 NSGA-II 内核）
│   ├── nsga2_core.py            #   DynamicDriver 主循环 + make_nsga2 内核
│   ├── ldmt_nsga2.py            #   LDT-NSGA-II（均值/方差迁移 + 50/20/30）
│   ├── distribution_features.py #   特征提取/非支配筛选/MigrationFunction/transfer_population
│   ├── environment_manager.py   #   EnvironmentManager（环境分段）+ ChangeDetector（探针检测）
│   ├── dnsga2_response.py       #   DNSGA-II 基线
│   ├── ar_nsga2.py              #   AR-NSGA-II 基线
│   ├── gmm_transfer_nsga2.py    #   GMM-Transfer 基线
│   ├── pps_nsga2.py             #   PPS 预测基线 (v2.1，v2.2 修复坐标 bug)
│   └── transfer_strategy.py     #   可插拔迁移策略接口 (v2.1)
├── problems/                    # 动态测试问题与真实 Pareto 前沿
│   ├── fda1.py / fda3.py / fda4.py
│   ├── dmop1.py / dmop2.py / dmop3.py
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
│   └── run_supplementary.py     #   补充实验（DMOP1/2/3）
├── configs/                     # 实验配置 YAML 镜像
├── documents/                   # 论文与答辩材料（documents/论文材料/）
│   ├── LDT-NSGA-II_毕业论文.docx
│   ├── LDT-NSGA-II V2.4毕业论文.md   # 论文工作稿（更新至 v2.3 数据、v2.4 修复状态）
│   └── LDT-NSGA-II V2.2答辩.html
├── tests/                       # pytest 单元与冒烟测试（63 个）
├── results/
│   ├── raw/                     # 逐代 CSV（comparison_/ablation_/sensitivity_/supplementary_*.csv）
│   ├── raw_legacy_v21/          # v2.1 旧数据归档（不参与统计）
│   ├── raw_legacy_v22_k1/       # v2.2 三基线旧数据归档（K1 修复前，不参与统计）
│   ├── summary/                 # 14 份汇总表/统计检验/图表 JSON/analysis_report.md
│   └── figures/                 # 12 张 PNG（箱线图/恢复曲线/耗时/消融/敏感性）
├── ldmt-experiment-report/      # HTML 实验报告（v2.4 起单文件自包含）
│   ├── index.html               #   build_report.py 生成（assets 已内联）
│   └── assets/                  #   doc.css / doc.js / echarts.min.js（内联源文件）
├── analyze.py                   # 分析流水线（汇总→统计→出图）
├── build_report.py              # HTML 报告组装器
├── main.py                      # 一键流水线入口（v2.4：测试→分析→报告→打开网页）
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

### LDT 方法

在每次环境变化时（入口 `respond()`，见 [algorithms/ldmt_nsga2.py](algorithms/ldmt_nsga2.py)）：

1. 从上一代种群中选出**非支配子集**（数量不足时退回全种群），
   非支配判定见 [non_dominated_set](algorithms/distribution_features.py)；
2. 在归一化空间提取其**均值 `mu_e` 与标准差 `sigma_e`** 作为分布特征，
   `sigma_min = 0.001*(u-l)`，归一化与特征提取见
   [extract_features](algorithms/distribution_features.py)；
3. 用固定探针点在新旧环境下的目标值差异计算变化幅度 `D_k`（仅记录），
   检测逻辑见 [ChangeDetector.detect](algorithms/environment_manager.py)；
4. **趋势中心预测**：`mu_hat = clip(mu_e + alpha*(mu_e - mu_{e-1}), 0, 1)`，
   默认 `alpha=0.5`（首次变化时 `mu_hat = mu_e`）；
5. **高斯扰动采样迁移候选**：`x_tr = clip(mu_hat + z*sigma_e, 0, 1)`，
   z~N(0,I)，再反归一化回真实决策空间；
6. 按 **50/20/30** 组成新一代种群（见 [algorithms/ldmt_nsga2.py](algorithms/ldmt_nsga2.py)）：
   - **50%** 迁移个体（`transfer_ratio`）
   - **20%** 保留旧种群个体（`retain_ratio`）
   - **30%** 均匀随机个体（`random_ratio`）

### 对比基线

| 算法 | 思路 |
|---|---|
| **NSGA-II**（static） | 静态基准，环境变化时不迁移，仅重新评估当前种群 |
| **DNSGA-II-inspired** | 变化检测 + 随机移民响应，无分布迁移 |
| **AR1-Predictive** | 一阶趋势预测基线（透明定义，非特定论文复现） |
| **GMM-inspired** | 高斯混合全分布对照（非原 KTMM 复现） |
| **PPS** | 简化版：种群中心一阶差分外推 + 整体平移 + 20% 随机移民（Zhou et al. 2014） |
| **LDT-NSGA-II** | 均值/标准差特征 + 趋势中心 + 50/20/30 混合初始化 |

> **K1 修复记录（v2.3）**：v2.2 及更早版本中，PPS/AR1/GMM 三个基线的
> `immigrant_ratio` 是死参数——"平移/采样后恰好 N 行 + 末尾 append
> 移民 + 定容截断前 N 行"的实现顺序导致随机移民被全部截断丢弃，三者
> 实际为"纯预测、无随机移民"行为。v2.3 已修复为**先预留 `n_imm` 个
> 名额、对预测个体切片后再拼接移民**（见
> [algorithms/ar_nsga2.py](algorithms/ar_nsga2.py)、
> [algorithms/gmm_transfer_nsga2.py](algorithms/gmm_transfer_nsga2.py)、
> [algorithms/pps_nsga2.py](algorithms/pps_nsga2.py)），并以
> `tests/test_k1_immigrant_fix.py` 3 个测试固化。基线行为变化属实质
> 修正，受影响的 1080 次运行已全部重跑，旧数据归档于
> `results/raw_legacy_v22_k1/`。

### v2.1 算法扩展

LDT-NSGA-II 另有三个可选模块（默认关闭，可通过参数开启）：

- **自适应迁移强度** (`adaptive_transfer=True`)：根据变化幅度 `D_k`
  动态缩放采样标准差，大变化时扩大探索范围，小变化时收紧局部搜索；
- **多源知识融合** (`multi_source=True`)：维护最近 K 个环境的特征
  缓冲区，按相似度加权融合 `w_i ∝ exp(-||mu_e - mu_i||/sqrt(D))`；
- **可插拔 TransferStrategy** (`transfer_strategy=...`)：将"特征提取→
  迁移计算→种群重建"封装为策略对象，支持 Moment/GMM 替换 —
  见 [transfer_strategy.py](algorithms/transfer_strategy.py)。

### 测试问题与指标

- **主实验问题**：FDA1、FDA3、FDA4（2/2/3 目标）
- **补充问题**：DMOP1、DMOP2、DMOP3（均为项目变体/自定义实现，非
  Goh & Tan 2009 原版，见 docstring；DMOP3 为静态对照组——K2 已在
  v2.4 如实标注并清理死代码）
- **主指标：稳态 MIGD**。口径以 [analyze.py](analyze.py) 的
  `per_run_metrics` 实现为准：对一次运行，取**每个环境段最后一代**
  （下一次变化前的末代种群）的 IGD 再对全部环境段取算术平均。论文与
  所有统计检验均使用该口径；辅助口径还有 post-change MIGD（每段首代
  IGD 平均）与全代 MIGD（所有代 IGD 平均）。
- **其他指标**：MHV（越大越好）、恢复代数、迁移模块耗时。
- **统计检验**：Wilcoxon 秩和 + Holm 校正；Friedman + Nemenyi 后验
  （多算法多问题排序，9 个"问题×τ"组合）。

## 运行原理

### 环境分段与变化检测

- `EnvironmentManager` 以 `tau_t` 代为一个环境段，`t = floor(g/tau_t)/n_t`
  为连续环境时间 — 见 [environment_manager.py](algorithms/environment_manager.py)；
- 每段起始（`g % tau_t == 0`）触发一次环境变化；
- `ChangeDetector` 维护一组固定探针点，比较其在新旧环境下的目标向量，
  得到归一化变化幅度 `D_k = mean(|new_F - prev_F| / scale)`。

### 趋势中心预测

- 维护上一环境的均值 `mu_prev` 与标准差 `sigma_prev`（均在归一化空间）；
- 当历史不少于 2 个环境时：

  ```
  mu_hat = clip(mu_e + alpha * (mu_e - mu_{e-1}), 0, 1)
  ```

  默认 `alpha=0.5`（敏感性测试 `{0, 0.25, 0.5, 0.75, 1.0}`）；
  首次变化时 `mu_hat = mu_e`；
- 标准差下限：`sigma_min = 0.001 * (u - l)`，v2.2 另对缩放设置上下限。

### 种群迁移与组成

- 迁移候选采样：归一化空间 `z_tr = clip(mu_hat + z * sigma_e, 0, 1)`，
  z~N(0,I)，再反归一化；
- 新种群按 **50/20/30** 组成：
  - **50%** 迁移个体（`transfer_ratio=0.5`）
  - **20%** 保留旧种群（`retain_ratio=0.2`）
  - **30%** 均匀随机（`random_ratio=0.3`）
- 非有限迁移值自动替换为均匀随机样本。

### NSGA-II 内核与逐代评估

- 每个环境段创建**全新**的 NSGA-II 实例（SBX `eta=15`、PM `eta=20`），
  以迁移结果作为初始种群，运行 `tau_t` 代 —
  见 [nsga2_core.py](algorithms/nsga2_core.py)；
- 段 `k` 使用 `seed = run_seed + k`，保证可复现且段间独立；
- **逐代**记录 IGD、HV、单代耗时、迁移耗时（仅在段首代记录一次），
  用于恢复曲线与耗时统计。

### 分析与报告产出

- [analyze.py](analyze.py)：汇总每轮指标 → 均值/标准差 → Wilcoxon 秩和
  检验 + Holm 校正 + Friedman + Nemenyi 后验检验 + 恢复代数分析 →
  箱线图、恢复曲线、耗时对比、消融与敏感性图；
- [build_report.py](build_report.py)：读取 `results/summary/` 下的图表
  JSON 组装 HTML 报告。v2.4 起（K12 修复）构建时内联
  `assets/doc.css` 与 `assets/echarts.min.js`，生成的 index.html 为
  **单文件自包含**，可直接分发。

## 实验结果（v2.3 数据 — K1 修复重跑版；v2.4 无实验行为变化）

> 数据状态：全部为 **30 个对齐种子**的完整结果（敏感性扫描项除外），
> 统一协议 N=100、n_gen=300、n_t=20。**3780 次运行 0 失败**：
> 主比较 1620 次（54 CSV）、补充 540 次（18 CSV）、消融 360 次
> （12 CSV）、敏感性 1260 次（42 CSV）。其中 PPS/AR1/GMM 三基线的
> 1080 次运行为 K1 修复后重跑（v2.2 对应数据已归档于
> `results/raw_legacy_v22_k1/`，不参与统计）。下表数字全部直接来自
> `results/summary/` 真实汇总。

### 主比较：稳态 MIGD（越低越好，30 种子均值）

| 问题 | τ | PPS | AR1 | **LDT** | GMM | static | DNSGA2 |
|---|--:|--:|--:|--:|--:|--:|--:|
| FDA1 | 10 | **0.01225** | 0.01729 | 0.01560 | 0.02027 | 0.02037 | 0.02109 |
| FDA1 | 25 | **0.00620** | 0.00651 | 0.00741 | 0.00843 | 0.00865 | 0.00865 |
| FDA1 | 50 | **0.00494** | 0.00499 | 0.00520 | 0.00526 | 0.00532 | 0.00536 |
| FDA3 | 10 | **0.01319** | 0.01838 | 0.01766 | 0.02204 | 0.02152 | 0.02273 |
| FDA3 | 25 | **0.00658** | 0.00696 | 0.00794 | 0.00900 | 0.00915 | 0.00923 |
| FDA3 | 50 | **0.00513** | 0.00518 | 0.00540 | 0.00551 | 0.00554 | 0.00553 |
| FDA4 | 10 | **0.11668** | 0.15071 | 0.12362 | 0.15951 | 0.21009 | 0.16279 |
| FDA4 | 25 | **0.08812** | 0.09872 | 0.09039 | 0.10110 | 0.10197 | 0.10272 |
| FDA4 | 50 | **0.07355** | 0.07463 | 0.07455 | 0.07567 | 0.07641 | 0.07586 |

### 核心结果（如实）：K1 修复后 PPS 全面第一，LDT 稳居第二

- **PPS 在全部 9 个"问题×τ"组合上稳态 MIGD 均为第一**——K1 修复后
  20% 随机移民注入生效，弥补了纯中心外推在多样性上的短板；
- **LDT 在全部 9 个组合上均为第二**：对 static/DNSGA2/GMM 保持显著
  优势（Wilcoxon-Holm 多数 p_adj≈9.3e-9），与 AR1 互有胜负，显著
  劣于 PPS；
- FDA4 τ=10 档：PPS 0.1167 vs LDT 0.1236（差距 5.9%），但 **LDT 的
  标准差 0.0047 小于 PPS 的 0.0075**，稳定性更好；Friedman 平均秩
  PPS 1.1、LDT 1.9，秩差 0.8 小于临界差 CD=1.377，Nemenyi 意义下
  两者不可区分；
- MHV 口径方向一致：FDA4 τ=10 时 PPS 0.5976、LDT 0.5707 居前二，
  第三名 DNSGA2 仅 0.4855。

### LDT 相对最佳基线（PPS）的差距与统计判读

| 问题 | τ=10 | τ=25 | τ=50 |
|---|--:|--:|--:|
| FDA1 | +27.4% | +19.6% | +5.3% |
| FDA3 | +33.9% | +20.7% | +5.3% |
| FDA4 | +6.0% | +2.6% | +1.4% |

（正=劣于 PPS，来自 `analysis_report.md` 1.1 节）

- **Wilcoxon-Holm 配对检验**：LDT vs PPS 在 9 个组合全部显著劣后
  （FDA4 τ50 最缓和，p_adj=0.0133，*；其余多为 p_adj≈9.3e-9）；
- **Friedman-Nemenyi 检验**：LDT vs PPS 在 FDA1 τ10（秩差 1.3）与
  FDA4 全部三档（0.8/0.73/1.0）上秩差小于 CD=1.377、不可区分；
  在 FDA1 τ25/50 与 FDA3 全部三档显著劣后；
- LDT vs AR1：FDA1 τ10 LDT 显著更好（p_adj=0.0013，**）；FDA3 τ10
  与 FDA4 τ50 无显著差异；其余档 AR1 显著更好；
- 9 个组合 Friedman 检验全部显著（χ² 范围 64.2–139.4，p<1e-9；
  JSON 中 `p_value=0.0` 为数值下溢）。

### 恢复速度与计算开销（τ=25，恢复代数越低越好）

| 问题 | PPS | AR1 | LDT | static | DNSGA2 | GMM |
|---|--:|--:|--:|--:|--:|--:|
| FDA1 | **10.29** | 16.02 | 21.29 | 21.86 | 21.98 | 22.01 |
| FDA3 | **13.97** | 17.03 | 21.64 | 22.13 | 22.13 | 22.11 |
| FDA4 | **9.65** | 16.05 | 18.50 | 20.64 | 20.67 | 20.69 |

单次迁移耗时（30 种子均值）：static 约 1e-4–9e-4 s，DNSGA2 约
3e-4–2.5e-3 s，LDT/AR1/PPS 约 0.006–0.05 s；GMM 约 0.006–0.62 s
（τ=10 档可达 0.6 s），**同配置配对**下为 LDT 的约 **4.4–11 倍**
（τ=10 约 10–11 倍、τ=25 约 8–9 倍、τ=50 约 4–7 倍）而其 MIGD
不占优——H2（低维特征更省时）成立；而"恢复最快"属于 PPS，不属于
LDT，这是对本方法不利但如实报告的结果。

### 消融实验（τ=25，稳态 MIGD）

| 问题 | full 完整 LDT | A：关趋势外推（α=0） | B：关分布迁移（改均匀随机） | C：κ=0（不保留旧个体） |
|---|--:|--:|--:|--:|
| FDA1 | **0.00741** | 0.00884 | 0.00965 | 0.00746 |
| FDA3 | **0.00794** | 0.00943 | 0.01027 | 0.00788 |
| FDA4 | **0.09039** | 0.10106 | 0.10223 | 0.09018 |

1. **分布迁移贡献最大**：去掉后退化 13–30%（变体 B 最差）；
2. **趋势外推有效**：去掉后退化 12–19%（变体 A）；
3. **旧种群保留（κ=20%）在 FDA 系列无可检测增益**：变体 C 与 full
   统计等价，FDA3/FDA4 均值甚至略优——默认 50/20/30 中的"20"缺乏
   本实验范围内的证据支持，需在多峰/非平稳场景再验证。

### 补充实验：DMOP1–3（τ=25，稳态 MIGD）

| 问题 | PPS | AR1 | LDT | GMM | static | DNSGA2 |
|---|--:|--:|--:|--:|--:|--:|
| DMOP1 | **0.00739** | 0.00789 | 0.00879 | 0.01019 | 0.01015 | 0.01039 |
| DMOP2 | **0.00688** | 0.00731 | 0.00822 | 0.00966 | 0.00972 | 0.00983 |
| DMOP3 | 0.19199 | 0.19200 | 0.19212 | 0.19220 | **0.19150** | **0.19150** |

DMOP1/DMOP2 与二目标 FDA 结论一致（PPS 最优，AR1 第二，LDT 第三）。
**DMOP3 六算法极差不足 0.8%**：K2 已在 v2.4 定性——该问题为项目
自定义**静态**欺骗问题（PS/PF 均不随时间变化），作为"迁移机制在
静态场景是否有害"的对照组保留；其上的算法比较不具动态优化区分度。
以 Goh & Tan (2009) 权威 dMOP3（Type I）重跑已列入展望。

### 参数敏感性（τ=25，跨问题/逐问题稳态 MIGD）

- 种群规模 N：50 → 0.0529，100（默认）→ 0.0352，200 → 0.0239；
- 环境严重度 n_t：10 → 0.0376，20（默认）→ 0.0352，40 → 0.0330；
- 迁移比例 ρ：0.25 →（0.00767/0.00826/0.09297），0.50（默认）→
  （0.00741/0.00794/0.09039），0.75 →（0.00735/0.00779/0.08890）
  （顺序为 FDA1/FDA3/FDA4）；
- 趋势步长 α：

  | α | FDA1 | FDA3 | FDA4 |
  |--:|--:|--:|--:|
  | 0 | 0.00884 | 0.00943 | 0.10106 |
  | 0.25 | 0.00826 | 0.00879 | 0.09568 |
  | 0.50（默认） | 0.00741 | 0.00794 | 0.09039 |
  | 0.75 | 0.00656 | 0.00702 | 0.08685 |
  | **1.00** | **0.00621** | **0.00664** | **0.08583** |

**不利发现**：α=1.0（全步外推）在三个问题上一致优于默认 α=0.5，
改善约 16%/16%/5%。论文已将"以 α=1.0 重跑主比较"列为第一优先展望，
以区分 LDT 与 PPS 的差距究竟来自外推步长还是矩迁移范式。

> 口径说明：N 与 n_t 扫描见 `sensitivity_runs.csv`；ρ/α 扫描数值由
> `results/raw/` 中对应原始 CSV 按"每环境段末代 IGD"重新聚合得到
> （α=0 与消融变体 A 交叉验证一致）。

### 图表索引（`results/figures/`，12 张 PNG）

| 图表文件 | 内容 | 数据来源 |
|---|---|---|
| `boxplot_MIGD_FDA1.png` | FDA1 各算法稳态 MIGD 分布箱线图（三档 τ） | `summary_comparison.csv` |
| `boxplot_MIGD_FDA3.png` | FDA3 各算法稳态 MIGD 分布箱线图 | 同上 |
| `boxplot_MIGD_FDA4.png` | FDA4 各算法稳态 MIGD 分布箱线图 | 同上 |
| `boxplot_MIGD_DMOP3.png` | DMOP3 各算法稳态 MIGD 分布（六算法重合，K2） | `supplementary_summary.csv` |
| `MIGD_vs_tau.png` | 稳态 MIGD 随 τ∈{10,25,50} 的变化趋势 | `summary_comparison.csv` |
| `migration_time.png` | 各算法单次迁移模块耗时对比（验证 H2） | `comparison_runs.csv` |
| `recovery_FDA1.png` / `recovery_FDA3.png` / `recovery_FDA4.png` | FDA 逐问题恢复曲线（τ=25，变化后 IGD 随代数变化） | `comparison_runs.csv` |
| `recovery_DMOP3.png` | DMOP3 恢复曲线 | `supplementary_runs.csv` |
| `ablation_MIGD.png` | 消融实验分组条形图（full / A / B / C） | `ablation_summary.csv` |
| `sensitivity.png` | 参数敏感性（N、n_t、ρ、α 扫描） | `sensitivity_runs.csv` + ρ/α 原始 CSV |

### 汇总数据文件索引（`results/summary/`，14 份）

| 文件 | 内容 |
|---|---|
| `summary_comparison.csv` | 主对比均值/标准差/中位数（每 问题×τ×算法 一行） |
| `comparison_runs.csv` | 主对比逐运行明细（30 种子） |
| `statistical_tests.csv` | Wilcoxon 秩和检验 + Holm 校正结果 |
| `friedman_results.json` | 9 组 Friedman χ²/平均秩 + Nemenyi CD 与显著对 |
| `recovery_generations.csv` | 恢复代数明细（每 问题×算法×种子×环境 一行） |
| `ablation_summary.csv` / `ablation_runs.csv` | 消融实验汇总 / 逐运行明细 |
| `supplementary_summary.csv` / `supplementary_runs.csv` | DMOP1/2/3 汇总 / 逐运行明细 |
| `sensitivity_runs.csv` | N、n_t 敏感性扫描明细 |
| `analysis_report.md` | 自动生成的完整文字分析报告 |
| `_chart_data.json` / `_recovery.json` / `_boxplot.json` | HTML 报告数据源 |

### 关键发现（有边界的双向结论）

- **K1 修复后的整体格局**：PPS 稳态 MIGD 在全部 9 个主比较组合与
  DMOP1/DMOP2 上均为最优——修复后的 20% 随机移民注入使其中心外推
  获得多样性补充；LDT 在全部组合上位居第二；
- **LDT 的可辩护位置**：对 static、DNSGA2、GMM 保持显著优势（多数
  p_adj≈9.3e-9）；与 AR1 互有胜负（FDA1 τ10 显著更好）；与 PPS 的
  差距随 τ 增大收敛到 5% 左右，且在 FDA4 全部三档与 FDA1 τ10 上
  Nemenyi 意义下与 PPS 不可区分；FDA4 上 LDT 的种子间标准差更小
  （τ=10：0.0047 vs PPS 0.0075），稳定性更优；
- **消融**：分布迁移贡献最大（去掉退化 13–30%），趋势外推次之
  （12–19%），旧种群保留无可检测增益；
- **对设计不利的发现同样如实报告**：α=1.0 一致优于默认 α=0.5；
  恢复速度 PPS 最快而非 LDT；DMOP3 与时间无关、六算法等价（K2）；
- **效率**：GMM 单次迁移耗时约为 LDT 的 4.4–11 倍且性能不占优，
  H2 成立；
- 完整逐代数据见 `results/raw/`，交互报告见
  [ldmt-experiment-report/index.html](ldmt-experiment-report/index.html)
  （v2.4 起单文件自包含，K12 已修复）。

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

v2.4 共 63 个测试，全部通过是数据发布前提。

### 2. 实验

所有实验均支持**断点续跑**：输出 CSV 中已存在的组合会被自动跳过。

```bash
# 主对比实验（FDA1/3/4 × tau{10,25,50} × 6 算法 × 30 种子 = 1620 次）
python -m experiments.run_comparison --problems FDA1,FDA3,FDA4 --taus 10,25,50 --seeds 0-29 --out results/raw

# 单配置快速冒烟测试
python -m experiments.run_comparison --problems FDA3 --taus 25 --seeds 0 --out results/raw --n-gen 50 --limit-seeds 1

# 消融实验（3 问题 × 4 变体 × 30 种子 = 360 次）
python -m experiments.run_ablation

# 敏感性实验（N / n_t / rho / alpha 扫描，共 1260 次）
python -m experiments.run_sensitivity

# 补充实验（DMOP1/2/3 × 6 算法 × 30 种子 = 540 次）
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
| `--n-t` | `20` | 变化严重度 `n_t` |
| `--algos` | 全部算法 | 逗号分隔的算法键（6 个，含 PPS） |
| `--limit-seeds` | `None` | 仅取前 N 个种子（快速测试） |

#### 其他实验入口

- `run_ablation`：问题 FDA1/FDA3/FDA4、τ=25、30 种子。消融变体：
  `full`、`ablation_a`（关趋势，α=0）、`ablation_b`（关分布迁移，
  迁移部分替换为均匀随机）、`ablation_c`（κ=0，不保留旧个体）；
- `run_sensitivity`：扫描 N ∈ {50,100,200}、n_t ∈ {10,20,40}、
  ρ ∈ {0.25,0.5,0.75}、α ∈ {0,0.25,0.5,0.75,1.0}；
- `run_supplementary`：问题 DMOP1/DMOP2/DMOP3、τ=25、6 个算法、
  30 种子。

### 3. 分析

```bash
python analyze.py
```

读取 `results/raw/*.csv`，输出逐运行指标、汇总表、Wilcoxon+Holm 与
Friedman+Nemenyi 统计检验、恢复代数以及图表到 `results/summary/` 与
`results/figures/`，同时生成 `results/summary/analysis_report.md`。

### 4. 生成 HTML 报告

```bash
python build_report.py
```

读取 `results/summary/_chart_data.json`、`_recovery.json`、
`_boxplot.json`，组装为 `ldmt-experiment-report/index.html`。
v2.4 起（K12 修复）CSS/JS/ECharts 全部内联，**单文件自包含**，
可直接拷贝分发。

### 5. 浏览报告

可直接在浏览器打开，或启动本地服务：

```bash
cd ldmt-experiment-report; python -m http.server 8000
```

## 算法一览

| 键 | 名称 | 说明 |
|---|---|---|
| `nsga2_static` | NSGA-II | 静态基准，无迁移 |
| `dnsga2_inspired` | DNSGA-II-inspired | 变化检测 + 随机移民响应 |
| `ar1_predictive` | AR1-Predictive | 一阶趋势预测基线（v2.3 修复随机移民注入） |
| `gmm_inspired` | GMM-inspired | 高斯混合全分布对照（v2.3 修复随机移民注入） |
| `pps` | PPS | 简化版种群中心一阶差分外推 + 整体平移 + 20% 随机移民（v2.3 修复移民注入） |
| `ldt_nsga2` | LDT-NSGA-II | 均值/标准差特征 + 趋势中心 + 50/20/30 混合初始化 |

## 说明与已知局限

- **FDA 公式**遵循 PlatEMO master（FDA1/FDA3/FDA4）。FDA4 使用标准
  DTLZ2 第一卦限单位超球面前沿（已验证：在 PS `x[:,M:]=G` 处目标范数
  恰好为 1.0）。
- **DMOP1/DMOP2/DMOP3** 均为**项目变体/自定义实现**，非
  Goh & Tan (2009, IEEE TEVC 13(1):103–127) 原版（v2.4 已更正归属：
  旧文档误引 Farina 2004；权威定义见 Helbig 2012 博士论文 p.52–54
  式 3.38–3.40）。PlatEMO master（commit d25e65d）经 git-tree API
  递归核查确认未提供 DMOP1/2/3.m（K11）。差异已逐条写入各文件
  docstring，结果按"补充实验"口径引用，不得称为标准 dMOP。
- **DMOP3（K2，v2.4 已定性）**：项目自定义静态欺骗问题，PS/PF 均不
  随时间变化，作为静态对照组保留；六算法稳态 MIGD 极差不足 0.8%，
  其上的比较不具动态优化区分度。
- **GMM-inspired** 为本项目对比实现，并非复现某篇特定论文的 KTMM。
- **LDT-NSGA-II** 算法参数（ρ=0.5、κ=0.2、α=0.5）为开题报告 2.0 默认值，
  非文献证明最优值；敏感性实验已显示 α=1.0 一致更优，重比较列入展望。
- HV 使用各问题固定参考点（原始目标值，不做逐轮归一化）：FDA1
  `(1.1,1.1)`、FDA3 `(2.2,1.1)`、FDA4 `(1.1,1.1,1.1)`。
- 所有随机性均设种子；环境段 `k` 使用 `seed = run_seed + k`；主比较
  与消融/补充实验使用同一套 30 个对齐种子。
- 实验结果仅适用于实际测试的问题、参数和硬件范围，不直接外推到所有
  真实工程任务。论文以 K1–K12 清单完整披露了全部已知实现问题：
  **K1 在 v2.3 修复并重跑；K2–K12 在 v2.4 全部处置**（修复或如实
  标注，无实验行为变化、0 次重跑，详见上文版本说明）。
