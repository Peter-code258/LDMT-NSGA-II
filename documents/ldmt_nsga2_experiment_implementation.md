# LDMT-NSGA-II 实验实现方案

## Context（背景与目标）

开题报告《渐变环境下基于分布迁移的轻量级动态多目标优化算法》要求：以 NSGA-II 为基础，设计轻量级种群分布特征迁移算法 LDMT-NSGA-II，在 FDA1/FDA3/FDA4（及 dMOP1/2 补充）上与 DNSGA-II、AR 预测、GMM 全分布迁移三类基线对比，并以 IGD/HV/单代时间/恢复曲线评价。本方案按报告附录 B 的目录结构、5.3 节的算法公式、5.5 节的实验配置落地代码与跑批，产出可复现的 CSV/图/统计报告。

用户已确认：① 完整 30 次/组（后台执行）；② 安装并使用 pymoo 作为 NSGA-II 内核。

## 关键架构决策

1. **统一 NSGA-II 内核**：所有 4 个算法共用 `pymoo.algorithms.moo.nsga2.NSGA2`（SBX 交叉 + 多项式变异 + 二元锦标赛选择 + 快速非支配排序 + 拥挤距离），仅响应机制不同——保证报告 5.5.2/5.7 要求的"同一套 NSGA-II 算子"。
2. **自定义动态驱动循环**：pymoo 无逐环境分段 FDA 管理，故按报告 5.7 在内核之上实现驱动器——每个环境段：设时间参数 t → 重评估种群/探针 → 触发响应生成新初始种群 → `minimize(problem, NSGA2(sampling=...), n_gen=τ_t)` 跑 τ_t 代。利用 pymoo 的 `sampling=X`（数组）自动重评估暖启动个体。
3. **FDA 自定义实现**：按 Farina 等(2004)原始公式实现为 `pymoo.core.problem.Problem` 子类，固定变量数/边界/时间公式，并核对解析真实 PF。
4. **统计与可视化**：scipy.stats.wilcoxon + Holm 校正；matplotlib 出图；全部指标从原始 CSV 生成。

## 文件清单（对应附录 B）

### 配置与依赖
- `requirements.txt` — pymoo(>=0.6)、numpy、scipy、scikit-learn、matplotlib、pyyaml、tqdm
- `configs/main_experiment.yaml` — 3问题×3τ_t×4算法×30次, n_t=20, N=100, 300代
- `configs/ablation.yaml` — FDA1/3/4, τ_t=25, 3变体(完整/移除特征简化/移除迁移函数), 30次
- `configs/sensitivity.yaml` — FDA3: N∈{50,100,150}; n_t∈{10,20,40}; 各30次

### 测试问题 `problems/`
- `fda1.py` — FDA1(M=2,n=20): PS 随 t 移动, PF 稳定(f1+f2=1)
- `fda3.py` — FDA3(M=2,n=20): PS 与 PF 同时变化
- `fda4.py` — FDA4(M=3,n=12): 三目标, PF 为平面
- `dmop1.py`, `dmop2.py` — 补充基准
- `true_fronts.py` — 各问题各环境的真实 PF 采样(解析或密集网格)

### 算法 `algorithms/`
- `nsga2_core.py` — `make_nsga2()` 工厂; `run_environment_segment()` 单段驱动; `DynamicDriver` 主循环
- `environment_manager.py` — 时间规则 t(g)=floor(g/τ_t)/n_t; 固定探针集检测变化; 归一化探针差分 D_k
- `distribution_features.py` — 边界归一化、均值 μ/逐维 σ 提取、迁移函数(中心/尺度外推)、迁移强度 λ_k(由 D_k 调节, 截断[0.5,1.5])
- `dnsga2_response.py` — DNSGA-II: 上一代种群(重评估) + 随机移民 + 变异响应
- `ar_nsga2.py` — AR(1) 中心预测: 近5环境 μ 序列拟合, 不稳则一阶差分回退
- `gmm_transfer_nsga2.py` — GMM 全分布迁移: 拟合 GMM 采样生成新初始种群(对照实现, 非原文复现)
- `ldmt_nsga2.py` — 所提算法: 70%迁移+20%局部变异+10%随机移民

### 评估 `evaluation/`
- `igd.py` — IGD: 真实 PF→解集最近距离均值
- `hv.py` — HV: 固定参考点(归一化目标尺度外留余量), pymoo 的 HV 指标
- `reference_front.py` — 各问题真实 PF 生成与缓存
- `statistics.py` — Wilcoxon 秩和 + Holm 校正; 均值/标准差/中位数/IQR; 汇总表

### 实验 `experiments/`
- `common.py` — 共享跑批框架(种子0-29、计时、CSV 写出、预热)
- `run_comparison.py` — 主对比(后台可执行)
- `run_ablation.py` — 消融
- `run_sensitivity.py` — 参数敏感性
- `run_supplementary.py` — dMOP1/2

### 测试 `tests/`
- `test_problem_dimensions.py` — 变量数/边界/目标维度
- `test_environment_change.py` — t 变化触发目标值改变、探针检测
- `test_indicators.py` — IGD/HV 对已知解集的正确性
- `test_smoke.py` — 各算法小规模(5代)可跑通

### 结果 `results/` — raw/summary/figures（运行时生成）
- 原始长表 CSV: algorithm/problem/seed/n_t/tau_t/N/env_idx/gen/IGD/HV/gen_runtime/migration_runtime/n_eval/code_version
- 汇总表 + 箱线图/恢复曲线/对比柱状图

## 核心算法实现要点（报告 5.3）

- **特征提取**：z=(x-l)/(u-l)；从非支配解(不足则全种群)算 μ_k、σ_k，成本 O(ND)。
- **迁移函数**：μ̂_k=μ_{k-1}+λ_k(μ_{k-1}-μ_{k-2})；σ̂_k=max(σ_min, σ_{k-1}+λ_k(σ_{k-1}-σ_{k-2}))；λ_k 由探针归一化变化量 D_k/D_{k-1} 调节并截断[0.5,1.5]，历史不足回退 λ=1。
- **个体映射**：ẑ=clip(μ̂_k + (σ̂_k/(σ_{k-1}+ε))⊙(z_i-μ_{k-1}),0,1) → 反归一化。
- **种群组成**：70%迁移 + 20%局部变异(对迁移个体加小高斯扰动) + 10%随机移民。比例待预实验冻结。

## 实验配置（报告 5.5）

| 实验 | 问题 | τ_t | n_t | N | 代数 | 算法 | 次数 |
|---|---|---|---|---|---|---|---|
| 主对比 | FDA1/3/4 | 10,25,50 | 20 | 100 | 300 | 4 | 30 |
| 消融 | FDA1/3/4 | 25 | 20 | 100 | 300 | 完整/去特征/去迁移 | 30 |
| 敏感-N | FDA3 | 25 | 20 | 50/100/150 | 300 | LDMT | 30 |
| 敏感-n_t | FDA3 | 25 | 10/20/40 | 100 | 300 | LDMT | 30 |
| 补充 | dMOP1/2 | 25 | 20 | 100 | 300 | 4 | 30 |

种子 0-29 固定；n_t=20；HV 参考点固定；先预热 1 次不计入。

## 执行步骤

1. 安装 pymoo 与依赖；建目录骨架。
2. 实现 FDA1/3/4 + 真实 PF；单元测试核对维度/边界/PF。
3. 实现 NSGA-II 内核驱动 + 环境管理器 + 探针检测；冒烟测试。
4. 实现 4 算法（DNSGA-II/AR/GMM/LDMT）；冒烟测试。
5. 实现 IGD/HV/统计；指标正确性测试。
6. 小预实验(5次×50代)冻结迁移比例等参数。
7. 后台执行主对比(30次)→消融→敏感性→补充，分批落 CSV。
8. 出图、统计检验、生成分析报告。

## 验证（Verification）

- **单元测试**：`pytest tests/` — 问题维度、环境变化触发、IGD/HV 对已知解集、各算法冒烟跑通。
- **预实验**：FDA3 上 5 次×50 代，确认 IGD 可算、HV 参考点固定、环境按 τ_t 触发、迁移模块计时非零。
- **数据校验**：每组运行数=预定值；IGD/HV 无 NaN/Inf；CSV 字段完整。
- **结果一致性**：所有图表源自原始 CSV；恢复曲线记录变化后第1/5代与下次变化前 IGD。

## 后台执行说明

主对比等大批次以 `run_in_background` 执行，分批写入 `results/raw/`；每批完成自动校验后再启动下一批，避免单次过长阻塞。预计总耗时数小时（用户已确认后台执行）。
