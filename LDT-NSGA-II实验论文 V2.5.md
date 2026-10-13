---
title: "面向动态多目标优化的轻量级分布迁移 NSGA-II：算法设计、实现与实验研究"
---

# 摘要

动态多目标优化问题（dynamic multi-objective optimization problem, DMOP）的目标函数或约束随时间变化，其 Pareto 最优解集（PS）与 Pareto 前沿（PF）随之移动，要求优化算法不仅能收敛，更要能**持续跟踪**变化后的最优区域。基于种群的方法在环境变化后往往因种群"记忆"过时而收敛变慢，因此**如何利用历史信息生成新环境的初始种群（warm-start）**，是动态多目标进化算法的核心问题。

本文在经典 NSGA-II 框架之上，设计并实现了一个轻量级的分布迁移算法 **LDT-NSGA-II**（Lightweight Distribution-Transfer NSGA-II）。其核心思想是：在归一化决策空间 $[0,1]^D$ 中提取上一环境非支配解集的**逐维均值 $\mu_e$ 与逐维标准差 $\sigma_e$** 作为分布特征，用趋势外推公式 $\hat\mu = \mu_e + \alpha(\mu_e - \mu_{prev})$ 预测新环境的分布中心，再按 **50% 迁移 + 20% 保留 + 30% 随机移民**的比例重组新环境的初始种群；另提供两个默认关闭的可选机制（基于变化幅度的自适应迁移宽度缩放、多源历史二阶矩融合）。单次环境响应的时间复杂度为 $O(MN^2 + ND)$，仅需 $2D$ 个分布参数，不引入概率模型拟合。

围绕该算法，本文构建了一套完全可复现的实验系统：5 个行为受控的对比基线（静态 NSGA-II、DNSGA-II 型响应、AR(1) 预测、GMM 全分布迁移、PPS 中心预测）、按原始文献实现的权威 dMOP1/dMOP2/dMOP3（覆盖 Type I/II/III 三种动态类型）、以及 6 个测试问题（FDA1/FDA3/FDA4 与 PlatEMO 源码逐式核对一致）。实验分两轮：v2.4 的四组实验（比较 1620、补充 540、消融 360、敏感性 1260，共 **3780 次运行记录**，0 失败）；v2.5 的关键对照实验（$\alpha=1.0$ 变体主比较 1620 次 + 权威 dMOP 补充 540 次）。统计推断采用配对 Wilcoxon 符号秩检验 + Holm 逐步校正，以及 Friedman + Nemenyi 检验（$k=6$、$N=30$、$CD=1.37668$）。

**主要结果（如实报告）**：

1. **v2.4（默认 $\alpha=0.5$）**：在全部 9 个"问题 × 变化周期 $\tau$"配置上，简化版 PPS（中心全步外推 + 20% 随机移民）稳态 MIGD 均值第一，LDT 稳居第二。LDT 在三个非预测基线上显著占优（多数 $p_{adj}=9.31\times10^{-9}$），在二目标规则平移问题（FDA1/FDA3/DMOP1/DMOP2）上显著劣于 PPS（$9/9$ 配置），但与 PPS 的差距随变化减慢从 27–34% 收敛到 5.3% 以内；在 Friedman–Nemenyi 整体排序意义下，二者在 FDA1 $\tau=10$ 与 FDA4 全部三档共 4 个配置上**不可区分**（秩差 $0.73$–$1.30 < CD$），且 LDT 在三目标 FDA4 上种子间方差更小（$\tau=10$：0.0047 对 0.0075）。消融实验显示"分布迁移"组件贡献最大（去掉后退化 13–30%）、"趋势外推"次之（12–19%）、"旧种群保留"无一致可检测增益。
2. **敏感性实验证伪了默认参数**：稳态 MIGD 随趋势系数 $\alpha$ 从 0 增大到 1.0 单调下降，$\alpha=1.0$ 在三个问题上一致最优（相对默认 $\alpha=0.5$ 改善 16%/16%/5%）。
3. **v2.5 关键对照实验（本文的核心机制发现）**：把 $\alpha$ 改为 1.0（全步外推）后，LDT 在 **11/12 个配置上相对 $\alpha=0.5$ 显著改善**，与 PPS 的差距从 $+1.4\%\sim+33.9\%$ 缩小到 $-8.5\%\sim+3.4\%$——其中 **4 个配置统计不可区分**（FDA1 全部三档、FDA4 $\tau=50$），**2 个配置显著反超**（FDA4 $\tau=10$：0.1068 对 0.1167，$p=1.86\times10^{-9}$；FDA4 $\tau=25$：0.0858 对 0.0881，$p=5.72\times10^{-7}$）。**因此，v2.4 中观察到的"LDT 显著劣于 PPS"主要来源于外推步长（半步 vs 全步），而非"矩迁移 vs 点预测"的范式差异。**
4. **v2.5 的第二个机制发现（一个明确的反例）**：唯一的例外是权威 **DMOP1A（Type III：Pareto 解集静止、仅前沿随 $t$ 旋转）**——在该问题上 $\alpha=0.5$ 以 0.0251 对 0.0570 **大幅领先** $\alpha=1.0$（$p=1.86\times10^{-9}$，相对改进 127%）。原因是解集不动时 $\mu_e-\mu_{prev}$ 几乎只含采样噪声，全步外推会把这个噪声放大成系统性偏移。**这一反例恰好验证了原设计中"半步更保守"的考虑**，并给出一个可操作的设计规则：**外推步长应随"解集是否在移动"自适应**——本文的数据已经把该规则的两个端点都验证过了。
5. **开销**：GMM 全分布迁移的单次响应耗时约为同配置 LDT 的 **4.4–11 倍**，且性能不占优；LDT 与 AR(1)、PPS 同处 $10^{-2}$ s 量级。
6. **诚实披露**：本文以 K1–K26 清单形式完整列出全部已知问题与不确定性——包括一个已被修复并重跑 1080 次运行的基线"死参数"缺陷（K1）、两处测试问题归属错误（K2/K3）、一处权威 dMOP 定义错误（v2.5 复核发现并更正：dMOP1 实为 Type III 而非 Type I）、版本戳漂移（K20）以及"3780 次运行"措辞的不精确（K21）。

本文的贡献可概括为三点：提出一种仅依赖相邻环境、无需训练模型的轻量级矩迁移机制并给出完整消融证据；设计并完成"**外推步长 × 迁移范式**解耦"的对照实验，把 LDT 与最优预测基线的差距重新归因，并发现"解集静止时外推有害"这一反例；构建口径统一、带 63 个单元测试与规范统计检验、并对全部已知问题如实披露的可复现实验系统。

**关键词**：动态多目标优化；进化算法；NSGA-II；分布迁移；预测策略；反向世代距离；Wilcoxon 检验；Friedman 检验

# Abstract

Dynamic multi-objective optimization problems (DMOPs) involve time-varying objectives whose Pareto-optimal set (PS) and Pareto front (PF) move over time, requiring an optimizer not only to converge but also to *track* the moving optima. Population-based methods tend to lag after an environmental change because the population's "memory" becomes stale; **how to generate a good initial population (warm start) for the new environment from historical information** is therefore the central question of dynamic multi-objective evolutionary algorithms.

Building on the classical NSGA-II framework, this thesis designs and implements **LDT-NSGA-II** (Lightweight Distribution-Transfer NSGA-II). In the normalized decision space $[0,1]^D$, the per-dimension mean $\mu_e$ and standard deviation $\sigma_e$ of the previous environment's non-dominated set are extracted as distribution features; a trend-extrapolation rule $\hat\mu = \mu_e + \alpha(\mu_e - \mu_{prev})$ predicts the new distribution centre; and the new initial population is composed of **50% transferred, 20% retained and 30% random-immigrant** individuals. Two optional mechanisms (change-magnitude-driven adaptive transfer-width scaling, and multi-source second-moment fusion) are provided but disabled by default. A single response costs $O(MN^2 + ND)$ and carries only $2D$ distribution parameters, with no probabilistic model fitting.

The algorithm is embedded in a fully reproducible experimental system: five behaviourally controlled baselines (static NSGA-II, DNSGA-II-style response, AR(1) prediction, GMM full-distribution transfer, PPS centre prediction); the authoritative dMOP1/dMOP2/dMOP3 implemented from the original literature (covering dynamic Types I, II and III); and six benchmark problems (FDA1/FDA3/FDA4 verified formula-by-formula against the PlatEMO sources). Experiments were run in two rounds: the v2.4 suite (comparison 1620, supplementary 540, ablation 360, sensitivity 1260 = **3,780 run records**, zero failures) and the v2.5 controlled study (an $\alpha=1.0$ variant comparison of 1620 runs plus 540 authoritative-dMOP runs). Statistical inference uses the paired Wilcoxon signed-rank test with Holm correction, and the Friedman + Nemenyi test ($k=6$, $N=30$, $CD=1.37668$).

**Main results (reported honestly).** (1) With the default $\alpha=0.5$, the simplified PPS baseline attains the best steady-state MIGD in all nine problem-by-$\tau$ configurations while LDT-NSGA-II ranks second everywhere: significantly better than the static, DNSGA-II-style and GMM baselines (most Holm-adjusted $p_{adj}=9.31\times10^{-9}$), yet significantly worse than PPS on the two-objective regular-shift problems (9/9 configurations), with the gap shrinking from 27–34% to within 5.3% as changes slow down. Under the Friedman–Nemenyi overall ranking the two are nonetheless **indistinguishable in four configurations** (FDA1 $\tau=10$ and all three FDA4 settings; rank gaps 0.73–1.30 $< CD$), and LDT has the smaller across-seed variance on the three-objective FDA4 ($\tau=10$: 0.0047 vs 0.0075). Ablation shows that distribution transfer is the largest contributor (13–30% degradation when removed), trend extrapolation the second (12–19%), and old-population retention has no consistent detectable benefit. (2) The sensitivity study falsifies the default parameter: steady-state MIGD decreases monotonically as the trend coefficient $\alpha$ grows from 0 to 1.0, and $\alpha=1.0$ is uniformly best on the three FDA problems (16%/16%/5% better than the default). (3) **The v2.5 controlled study is the central mechanistic finding of this thesis**: with full-step extrapolation ($\alpha=1.0$), LDT improves in **11 of 12 configurations**, and the gap to PPS shrinks from $+1.4\%\sim+33.9\%$ to $-8.5\%\sim+3.4\%$ — statistically indistinguishable in four configurations (all three FDA1 settings, FDA4 $\tau=50$) and **significantly better in two** (FDA4 $\tau=10$: 0.1068 vs 0.1167, $p=1.86\times10^{-9}$; FDA4 $\tau=25$: 0.0858 vs 0.0881, $p=5.72\times10^{-7}$). **The v2.4 observation that LDT is significantly worse than PPS therefore stems mainly from the extrapolation step size (half-step vs full-step) rather than from the "moment transfer vs point prediction" paradigm difference.** (4) **The second mechanistic finding is a clear counter-example**: the only exception is the authoritative **DMOP1A (Type III: the Pareto-optimal *set* is static while the *front* rotates with $t$)**, where $\alpha=0.5$ beats $\alpha=1.0$ by a wide margin (0.0251 vs 0.0570, $p=1.86\times10^{-9}$, a 127% relative improvement). The reason is that when the solution set does not move, $\mu_e-\mu_{prev}$ contains little more than sampling noise, and a full-step extrapolation amplifies that noise into a systematic offset. **This counter-example vindicates the original "conservative half-step" rationale and yields an actionable design rule: the extrapolation step should adapt to whether the solution set is actually moving — this thesis has validated both endpoints of that rule.** (5) GMM full-distribution transfer costs 4.4–11 times as much per response as LDT on matched configurations while performing no better. (6) Known issues and uncertainties are disclosed in full as a K1–K26 list, including a repaired baseline "dead-parameter" defect for which 1080 runs were re-executed (K1), two problem-attribution errors (K2/K3), an authoritative dMOP definition error found and corrected in v2.5 (dMOP1 is Type III, not Type I), a version-stamp drift (K20), and the imprecision of the "3,780 independent runs" wording (K21).

**Keywords**: dynamic multi-objective optimization; evolutionary algorithm; NSGA-II; distribution transfer; prediction strategy; inverted generational distance; Wilcoxon test; Friedman test


\newpage

# 第 1 章 绪论

## 1.1 研究背景

在工程设计、生产调度、能源管理、交通控制、投资组合等大量实际问题中，决策者往往需要同时优化多个彼此冲突的目标：既要降低成本又要提高质量，既要缩短工期又要节约资源。这类问题称为多目标优化问题（multi-objective optimization problem, MOP）。由于目标之间相互冲突，通常不存在使所有目标同时达到最优的单一解，而是存在一组由 Pareto 支配关系定义的非支配解——Pareto 最优解集（PS），及其在目标空间中的像——Pareto 最优前沿（PF）。决策者最终从这组折中解中依据偏好做出选择。

经典多目标进化算法如 NSGA-II [1] 通过非支配排序、拥挤度距离与精英保留机制，在**静态环境**中取得了很好的求解效果，并在此后二十余年间成为多目标优化的基准内核。然而现实世界很少是一成不变的：原材料价格波动、设备状态漂移、交通流量随时间变化、天气与负荷持续改变，都会使目标函数或约束随时间变化。此时问题由静态 MOP 扩展为**动态多目标优化问题**（dynamic multi-objective optimization problem, DMOP）：

$$\min_{x \in \Omega} F(x, t) = (f_1(x,t), f_2(x,t), \ldots, f_M(x,t))^{\mathsf T}$$

在 DMOP 中，$PS(t)$ 与 $PF(t)$ 随环境（时间）而移动、形变甚至改变拓扑结构。Farina、Deb 与 Amato 于 2004 年给出了 DMOP 的经典定义并构造了 FDA 系列动态测试问题 [2]，将动态变化区分为解集变化而前沿不变、前沿本身变化、两者均变化等不同类型。对一个动态优化算法而言，它不仅要能在某一时刻逼近当前前沿，还要在环境变化后**快速跟踪**并重新收敛到新的前沿——这对响应速度、多样性保持与计算开销同时提出了要求。

![图 1-1 动态多目标优化问题示意](../../LDT-NSGA-II/results/figures_thesis/fig1_1_dmop_concept.png)

图 1-1：动态多目标优化问题示意（以 FDA1 为例）。左/右分别为相邻两个时刻的 Pareto 前沿与种群位置，箭头表示前沿随时间的迁移方向；旧种群滞后于新前沿，环境变化后需要重新跟踪。

## 1.2 问题的提出

基于种群的方法在环境突变后往往因种群"记忆"过时而收敛变慢。由此，近二十年的研究集中于同一个问题：

> **环境变化后，如何生成一个好的初始种群（warm-start）？**

现有策略大致分为三条主线（详见第 2 章）：

1. **响应式（diversity-based）**：检测到变化后，向种群注入随机移民或变异个体以恢复多样性。代表工作是 Deb 等对 NSGA-II 的动态改造（DNSGA-II），提出随机移民与变异移民两种变体 [3]。该类方法无需历史样本，但对新环境没有任何方向信息。
2. **预测式（prediction-based）**：用历史种群统计量（如中心点序列）建立时序模型，外推新环境的种群位置。代表工作是 Zhou 等提出的种群预测策略（PPS）[4]。
3. **迁移/分布学习式（transfer / distribution-learning）**：显式建模历史种群的**分布**（而非仅中心），并把分布知识迁移到新环境。代表方向包括以高斯混合模型刻画动态最优解分布的知识迁移方法 [14]。

这三类方法在"精度—开销"谱系上各占一端：响应式几乎零开销但恢复慢；预测式在规则变化下精度高但需要保存多代历史并拟合模型，且预测假设失效时会引入系统性偏差；分布迁移式表达力强但概率模型的参数估计成本与"黑箱度"都更高。

由此提炼出本课题关注的三个缺口：

- **G1（轻量化）**：能否设计一种只依赖**相邻（或少量）环境**、不需要训练时序模型或概率模型、单次响应计算量可忽略的迁移机制？
- **G2（可审查性与组件归因）**：在"复杂模型 vs 简单规则"的取舍中，简化到何种程度仍能保留预测/迁移的收益？各组件（趋势外推、分布迁移、旧解保留、随机移民）各自的边际贡献能否被独立量化？
- **G3（可复现性）**：能否建立一套口径统一、带单元测试与规范统计检验的实验系统，使结论可被第三方复现与证伪？

## 1.3 研究意义

**理论意义。** 其一，把"渐变环境"从定性描述转化为可度量的操作化定义（第 2 章定义 2.1），有助于深化对"环境变化幅度—种群分布迁移—收敛性能"三者关系的理解。其二，将分布特征迁移与经典 NSGA-II 框架正交地结合，探索一种不依赖重型预测模型、不依赖大量历史环境样本的轻量级响应机制，可丰富动态多目标优化的方法体系。其三，通过**把外推步长与迁移范式分离**的对照实验（第 5 章），为"预测式方法的优势究竟来自哪里"这一机制性问题提供实验证据。

**工程价值。** 轻量级方法实现简单、参数较少、计算与存储开销小，便于在算力受限或需要在线运行的实际系统中部署；基于相邻环境分布迁移的快速恢复能力，有助于缩短环境变化后的重新收敛时间，提高优化系统在动态工况下的稳定性与可用性。

**必须客观说明的是**：本课题的验证主要在 FDA 与 dMOP 等标准动态测试问题上进行；从测试问题到具体工程场景的迁移仍需结合实际问题进一步验证，本文不夸大其适用范围（适用边界见 5.9 节与第 7 章）。

## 1.4 本文的主要工作与创新点

本文完成的主要工作如下。

1. **算法设计**：提出并实现 LDT-NSGA-II（Lightweight Distribution-Transfer NSGA-II）。其内核是在归一化决策空间 $[0,1]^D$ 上提取上一环境非支配解集的**逐维均值 $\mu_e$ 与逐维标准差 $\sigma_e$** 作为分布特征，用趋势外推公式 $\hat\mu = \mu_e + \alpha(\mu_e - \mu_{prev})$ 预测新环境的分布中心（默认 $\alpha = 0.5$），再按 **50% 迁移 + 20% 保留 + 30% 随机移民** 的比例重组新环境的初始种群。另提供两个默认关闭的可选机制：基于变化幅度的自适应迁移宽度缩放（上限 $[0.5, 2.0]$）与多源历史融合（混合二阶矩方差 + 相似度权重下限 0.05）。单次响应的复杂度为 $O(MN^2 + ND)$，不引入概率模型拟合（第 3 章）。
2. **受控基线与测试问题**：实现 5 个行为明确的对比基线（静态 NSGA-II、DNSGA-II 型响应、AR(1) 预测、GMM 全分布迁移、PPS 中心预测），其行为由单元测试锁定；实现 6 个测试问题，其中 FDA1/FDA3/FDA4 与权威平台 PlatEMO 的 MATLAB 源码**逐式核对一致**；并按原始文献实现**权威的 dMOP1/dMOP2/dMOP3**（v2.5 新增），覆盖 **Type I、Type II、Type III 三种动态类型**（这是本文测试集相对既有材料的重要补强），以解析自检脚本验证"解析 Pareto 集落在解析 Pareto 前沿上"（最大偏差 0.000e+00）。
3. **完全可复现实验系统**：四组实验（比较 1620 次、补充 540 次、消融 360 次、敏感性 1260 次）共 **3780 次独立运行记录、0 失败**，全部支持断点续跑；每次运行落盘逐代 18 字段 CSV 与运行级元数据（git 提交、库版本、逐环境评估数）。新增 v2.5 两阶段实验（$\alpha = 1.0$ 对照 1620 次 + 权威 dMOP 补充 630 次，共 2250 次），全部结果由独立脚本在本论文写作期间重新计算（第 4 章）。
4. **诚实的评估与组件归因**：63 个单元测试；统计推断采用配对 Wilcoxon 符号秩检验 + Holm 逐步校正，以及 Friedman + Nemenyi 检验；对**当前实现的真实行为**逐一如实记录于"已知问题清单"（K1–K26）。消融实验分离出各组件贡献（分布迁移 > 趋势外推 ≫ 旧解保留），敏感性实验揭示默认半步系数 $\alpha = 0.5$ 并非最优；v2.5 的 $\alpha = 1.0$ 对照实验进一步回答了"差距来自外推步长还是迁移范式"这一关键问题，并**在权威 dMOP1（Type III，解集静止）上发现了全步外推反而有害的反例**，据此提出自适应步长规则（第 5 章、第 6 章）。

本文的**创新点**可概括为三点：

- **C1**：提出一种在统一归一化空间中仅用"逐维均值 + 逐维标准差"两个矩特征、仅依赖相邻环境、无需训练模型的轻量级分布迁移机制，并给出完整的组件化设计与消融证据；
- **C2**：设计并完成**"外推步长 × 迁移范式"解耦的对照实验**，把 LDT 与最优预测基线（PPS）的差距归因从"范式差异"重新定位到"步长差异"，并以同一随机种子下的配对检验给出统计证据；
- **C3**：构建口径统一、带断点续跑与运行级元数据、含 63 个单元测试与规范统计检验的可复现实验系统，并以"已知问题清单"的形式完整披露实现缺陷与不确定项，使结论可被第三方复现与证伪。

## 1.5 论文组织结构

全文共 7 章，组织结构如下：

- **第 1 章 绪论**：研究背景、问题的提出、研究意义、主要工作与创新点、组织结构。
- **第 2 章 动态多目标优化研究现状与文献综述**：DMOP 的形式化与分类、渐变环境的操作化定义、三条响应策略主线的综述、测试问题与实验平台、统计检验规范，以及研究缺口提炼。
- **第 3 章 LDT-NSGA-II 算法与测试问题设计**：算法的完整设计（归一化特征、趋势外推、三路重组、两个可选机制、共享探针与变化幅度）、复杂度分析、五个基线的设计、六个测试问题（含权威 dMOP）的公式与解析 Pareto 集/前沿推导。
- **第 4 章 实验系统与可复现性设计**：软硬件环境、代码结构、四组实验的因子设计与预算核算、数据管线与产物、测试体系、断点续跑与元数据、以及"从零复现"的完整步骤。
- **第 5 章 实验结果与分析**：实验设置与数据状态；FDA1/FDA3/FDA4 上的主比较结果；跨问题综合分析（Friedman 平均秩、恢复速度、计算开销）；消融实验；DMOP 补充实验与权威 dMOP 实验；敏感性实验；v2.5 的 $\alpha = 1.0$ 对照实验与全部 12 个配置的步长对照；结果解读（事实层与机制假说层分离）；已知问题与不确定性清单（K1–K26）。
- **第 6 章 结论与展望**：结论、适用边界、局限性与未来工作。
- **附录 A–E**：零基础运行手册、术语表、代码文件索引与章节映射、完整可运行源代码、以及实验数据文件清单。

参考文献、致谢与学位论文独创性声明与使用授权说明按学校模板置于文末。


\newpage

# 第 2 章 动态多目标优化研究现状与文献综述

> 撰写说明：本章全部文献著录信息均以网络检索核实为准（检索时间 2026 年 10 月），无法核实的字段在文中显式标注。引用编号与文末参考文献表一致。

## 2.1 文献检索与筛选方法

为保证综述的覆盖面与可追溯性，本章按如下流程筛选文献：

1. **数据库与检索式**：以 IEEE Xplore、Web of Science、ACM Digital Library、CNKI 为主要来源；英文检索式采用 `("dynamic multi-objective" OR "dynamic multiobjective") AND ("evolutionary" OR "optimization")`，中文检索式采用 `动态多目标 AND (进化 OR 优化)`。时间范围以 2002—2026 年为主，并对 2002 年前的奠基性工作（如 NSGA、SPEA）做追溯引用。
2. **纳入标准**：(a) 提出或系统评估动态多目标优化算法；(b) 提出或修订动态测试问题与性能指标；(c) 提出多算法比较的统计检验规范；(d) 被本课题直接使用的软件平台与工具。
3. **排除标准**：仅处理静态多目标优化的文献（除作为内核来源引用外）；无法获取全文或权威著录信息不完整的条目（此类若必须引用，则在参考文献表中显式标注"未能核实"）。
4. **筛选结果**：共纳入 40 余条，按"基础理论—响应/预测/迁移方法—测试基准与指标—实验方法论—软件平台"五类组织。

> **诚实声明**：本课题为工程型研究，作者不声称对领域做了穷尽式系统综述。本章的目的是把本课题所用到的每一条技术路线、每一个测试问题、每一种统计口径的原始出处落到可核查的文献上；对于"本文未采纳但同属主流"的方向（如强化学习型 DMOEA、多群体协同、约束动态优化），仅做定位性说明，不做深入评述。

## 2.2 多目标优化与 NSGA-II 基础

多目标优化问题（multi-objective optimization problem, MOP）的一般形式为

$$\min_{x \in \Omega} F(x) = (f_1(x), f_2(x), \ldots, f_M(x))^{\mathsf T}, \quad \Omega \subseteq \mathbb{R}^D$$

其中 $M \ge 2$ 为目标数。当各目标相互冲突时，通常不存在使所有目标同时达到最优的单一解，因此解的比较由 Pareto 支配关系给出：若 $F(a) \le F(b)$ 逐分量成立且至少存在一个分量严格小于，则称 $a$ 支配 $b$。不被 $\Omega$ 中任何解支配的解构成 Pareto 最优解集（Pareto-optimal set, PS），其在目标空间的像称为 Pareto 最优前沿（Pareto front, PF）。

| 概念 | 符号 | 定义 | 项目内实现 |
|---|---|---|---|
| Pareto 支配 | $a \prec b$ | $\forall j:\ f_j(a) \le f_j(b)$ 且 $\exists j:\ f_j(a) < f_j(b)$ | `algorithms/distribution_features.py::non_dominated_set`（带 $10^{-8}$ 容差） |
| Pareto 最优解集 | $PS$ | $\{x \in \Omega \mid \nexists x' \in \Omega: x' \prec x\}$ | 同上（种群层面的近似） |
| Pareto 最优前沿 | $PF$ | $\{F(x) \mid x \in PS\}$ | `problems/true_fronts.py::get_true_pf` |

围绕 MOP 的求解，进化计算因其种群并行性（一次运行即可逼近整条前沿）成为主流范式。Srinivas 与 Deb 于 1994 年提出第一代非支配排序遗传算法 NSGA [15]；Zitzler 与 Thiele 于 1999 年提出强度 Pareto 进化算法 SPEA 并系统引入超体积（hypervolume, HV）指标 [9]；Deb 等 2002 年提出的 NSGA-II [1] 以三项机制奠定了此后二十年的工程标准：

1. **快速非支配排序**：以 $O(MN^2)$ 的代价将种群按支配层次分层（第一层即非支配解），取代 NSGA 中 $O(MN^3)$ 的排序；
2. **拥挤度距离**：在同层内以目标空间中的邻居间距度量解的疏密，作为非支配排序的次级判据，从而在不引入额外参数的前提下维持分布多样性；
3. **精英保留**：父代与子代合并后按"层号升序、层内拥挤度降序"择优进入下一代，保证已找到的非支配解不丢失。

NSGA-II 的算子配置（模拟二进制交叉 SBX [1] 与多项式变异 PM、非支配排序与拥挤度的选择）被本课题完整沿用作为所有算法的统一搜索内核，见第 4 章与附录 D 的可运行代码。需要强调的是：**NSGA-II 及其变体在本课题中只承担"单环境内的静态优化器"角色**；DMOP 的全部难点在于环境变化后如何重新初始化种群，这正是本章后续各节的焦点。

## 2.3 动态多目标优化问题的形式化与分类

### 2.3.1 形式化定义

动态多目标优化问题（dynamic multi-objective optimization problem, DMOP）在 MOP 的基础上引入时间变量 $t$：

$$\min_{x \in \Omega} F(x, t) = (f_1(x,t), f_2(x,t), \ldots, f_M(x,t))^{\mathsf T}$$

相应地，最优解集与最优前沿都成为时间的函数：$PS(t)$ 与 $PF(t)$。Farina、Deb 与 Amato 于 2004 年系统给出了 DMOP 的定义、FDA 系列测试问题与近似质量度量 [2]，将动态性区分为以下几类：

| 类型 | $PS(t)$ | $PF(t)$ | 典型问题 | 跟踪难点 |
|---|---|---|---|---|
| I | 随时间变化 | 不变 | FDA1、权威 dMOP3 | 前沿形状已知，但最优决策变量的位置在漂移 |
| II | 随时间变化 | 随时间变化 | FDA3、权威 dMOP2 | 需同时跟踪解集与前沿，收敛判据本身也在变 |
| III | 不随时间变化 | 随时间变化 | 较少见 | 前沿形变而解集静止 |
| IV | 均不随时间变化 | 均不随时间变化 | 静态问题的动态包装 | 无动态性（本课题的 DMOP3 项目变体即属此类，见 3.4.4 节） |

### 2.3.2 时间模型与变化形态

动态性在算法层面必须被离散化为可执行的时间模型。本课题采用与 PlatEMO [6] 及 FDA 原始定义 [2] 一致的**离散阶跃时间规则**：

$$t = \frac{1}{n_t}\left\lfloor \frac{g}{\tau_t} \right\rfloor$$

其中 $g$ 为 0 基代数，$\tau_t$ 为每个环境持续的代数（变化周期），$n_t$ 为一个运行期内不同环境状态的数目。该模型下，环境在 $g = k\tau_t$ 处发生**阶跃**，$t$ 的取值集合为 $\{0, 1/n_t, \ldots, (n_t-1)/n_t\}$。

值得单独讨论的是"突变"与"渐变"。文献中常把 $\tau_t$ 小、每次变化幅度大的情形称为"突变（abrupt）"，把 $\tau_t$ 大、单步幅度小的情形称为"渐变（gradual）"或"连续变化（continuous）"。本课题关注的**渐变**并非"时间连续并在积分意义上变化"，而是如下可度量的性质：

> **定义 2.1（渐变环境，可操作化定义）**：设相邻两个环境的 Pareto 最优解集为 $PS(t_k)$ 与 $PS(t_{k+1})$，定义单步变化幅度
> $$d_k = \frac{1}{|\mathcal S|}\sum_{x \in \mathcal S} \left\| x - \Pi_{PS(t_{k+1})}(x) \right\|_2$$
> 其中 $\mathcal S \subset PS(t_k)$ 为均匀采样点集，$\Pi$ 为到集合的投影。若对全部 $k$ 有 $d_k$ 远小于决策空间直径，且方向 $\mathrm{dir}_k$ 在相邻步之间具有一致性，则称该环境序列为渐变环境。

在本课题所用的 FDA/dMOP 基准中，$G(t) = \sin(\tfrac{\pi t}{2})$ 在 $t$ 的小步长下近似线性，故 PS 每步做近似匀速平移；$n_t$ 越大，单步 $d_k$ 越小、变化越平滑。这正是第 5 章把 $n_t$ 作为敏感性因子之一的原因，也是"基于相邻环境分布迁移"这一设计点的适用前提。

### 2.3.3 评价指标

动态场景下的评价必须同时刻画"收敛质量"与"跟踪速度"，常用指标如下。

- **反向世代距离 IGD**（inverted generational distance）：设 $P^*$ 为 PF 上均匀采样的参考点集，$P$ 为算法所得解集，则
  $$\mathrm{IGD}(P, P^*) = \frac{1}{|P^*|}\sum_{p \in P^*} \min_{q \in P} \|p - q\|_2$$
  IGD 同时惩罚收敛不足与分布不均，值越小越好。其通行引用来源见 Bosman 与 Thierens [8]；本课题实现见 `evaluation/igd.py`（欧氏距离由 SciPy 计算）。
- **超体积 HV**：解集与固定参考点围成的目标空间体积，具有严格的 Pareto 单调性，不依赖真实前沿 [9]。本课题实现封装 pymoo 的精确 HV。
- **动态扩展指标 MIGD**：对一次运行，按环境分段聚合 IGD。本课题采用的口径为"**稳态 MIGD**"——取每个环境段最后一代的 IGD，对全部环境段求平均；另报告"变化后 MIGD"（每段第 1 代）与"全代 MIGD"两个辅助口径。三个口径的确切代码定义见 5.1.3 节。
- **恢复曲线与恢复代数**：以环境内代数 $j$ 为横轴、跨种子平均 IGD 为纵轴，刻画变化后需要多少代回到稳态水平，其精确定义见 5.1.3 节。

## 2.4 环境变化响应策略的三条主线

环境变化后如何生成新环境的初始种群，是 DMOP 求解的核心问题。文献中的策略可归纳为三条主线（图 2-1）。

### 2.4.1 多样性保持与响应式方法

该类方法不做任何预测，仅在检测到变化后向种群注入新个体以恢复多样性。Deb 等在将 NSGA-II 改造用于水火电调度动态问题时提出了 DNSGA-II [3]，给出两种变体：**随机移民**（以 $[x_l, x_u]$ 上的均匀分布替换种群中一定比例的个体）与**变异移民**（对现有解施加高斯扰动）。该类方法的优点是实现简单、不需要历史样本、对任意变化形态都不会"预测反被预测误"；缺点是没有方向信息，丢失了上一环境已经获得的大部分收敛成果，在变化幅度较大时恢复缓慢。本课题将 DNSGA-II 型响应作为受控基线之一。

### 2.4.2 预测式方法

该类方法假设环境变化具有可学习的时间规律，通过历史最优解或其统计量建立时序模型，外推新环境的最优解位置。代表性工作包括：

- **前向预测策略**（forward-looking prediction，FPS）一类方法，用历史最优解的时间序列做自回归/线性外推，作为新环境的初始种子；
- **种群预测策略 PPS**（population prediction strategy）[4]：由 Zhou、Jin 与 Zhang 提出，其核心思想是**把种群表示为"中心 + 流形"两部分**，对中心序列做自回归外推、对流形做相应调整，再按外推结果生成新环境种群。PPS 是当前最具代表性的预测式 DMOEA 之一，也是本课题最重要的参照基线。

预测式方法在变化规律明显、变化平稳时优势显著，但其代价是：需要保存多个历史环境样本、需要拟合时序模型（计算与存储开销增加），并且**当预测假设不成立时会引入系统性偏差**。本课题第 5 章的实验正是围绕"外推步长"这一具体机制对该类方法做了解剖（5.5.4 节、5.6.4 节）。

### 2.4.3 知识迁移与分布学习式方法

该类方法把"复用过去的搜索经验"显式建模为跨环境的知识迁移，可进一步分为两个子层：

- **个体层面的迁移**：直接迁移历史最优个体或其线性组合/映射结果。Jiang 等提出的迁移式动态多目标进化算法（transfer-based DMOEA）通过构造非支配解在不同环境之间的映射来复用历史解；
- **分布层面的迁移**：对最优种群的**整体分布**（均值、协方差、概率密度、低维流形）建模，再迁移到新环境。Zou 等提出的 **KTMM**（knowledge transfer with mixture model）[14] 用高斯混合模型刻画动态最优解的分布并实现跨环境迁移，代表了该方向较新的进展。

分布层面的方法对个体对应关系不敏感、在渐变环境下平滑性较好，但需要估计概率模型的参数，开销与"黑箱度"都更高。本课题的 **LDT-NSGA-II 正是落在预测式与分布迁移式之间的一个轻量级折中**：只取"逐维均值 + 逐维标准差"两个矩特征（而非完整混合模型），只加一条确定性趋势外推（而非拟合时序模型），从而在恢复性能与计算/存储开销之间取平衡。

![图 2-1 动态多目标优化方法分类与本文定位](../../LDT-NSGA-II/results/figures_thesis/fig2_1_taxonomy.png)

图 2-1：动态多目标优化方法分类与本文定位。响应式、预测式、分布/知识迁移型三条主线及各自的代表算法；LDT-NSGA-II（高亮）位于预测式与迁移型之间的轻量级折中位置。

## 2.5 测试问题与实验平台

### 2.5.1 测试问题族

| 问题族 | 提出者与出处 | 特点 | 本课题使用情况 |
|---|---|---|---|
| FDA1–FDA5 | Farina、Deb、Amato 2004 [2] | 覆盖 Type I/II，时间参数 $G(t)=\sin(\tfrac{\pi t}{2})$ 族 | 使用 FDA1、FDA3、FDA4；FDA1/FDA3/FDA4 与 PlatEMO 源码逐式核对 |
| dMOP1–dMOP3 | Goh 与 Tan 2009 [5] | 面向共进化范式设计，含 Type I/II | 使用权威定义实现（见 3.4 节）及其对照 |
| JY 系列等 | 动态多目标优化的系统构造基准 | 可调变化频率与变化严重程度 | 列入展望（第 7 章） |

需要指出：**PlatEMO 平台并不包含 dMOP1–dMOP3 的实现**（本课题已通过 GitHub git-tree API 对 PlatEMO master 全树递归检索确认，见 5.8 节 K11），因此 dMOP 族的实现必须直接依据原始文献。这一问题在本课题的早期版本中确实导致过归属错误，第 5 章以"已知问题清单"的形式做了完整披露。

### 2.5.2 软件平台

- **PlatEMO** [6]：MATLAB 侧的进化多目标优化事实标准平台，集成大量问题、算法与指标。本课题以 PlatEMO 的 FDA 问题源码作为**逐式对照基准**。
- **pymoo** [7]：Python 侧的多目标优化框架，提供 NSGA-II 实现与 SBX/PM 算子。本课题的算法内核基于 pymoo 0.6.1.5。

### 2.5.3 一个被忽视的工程前提：索引与切片约定

在把 MATLAB 实现（PlatEMO，1 基索引）逐式转写为 Python（0 基索引）时，**切片边界是最容易出错的位置**。以 FDA4 为例，PlatEMO 中距离变量写作 `PopDec(:, M:end)`，对应 Python 应为 `X[:, M-1:]`（共 $D-(M-1) = 10$ 个变量），而非 `X[:, M:]`（9 个变量）。本课题在早期版本中曾出现该偏差，并在 v2.2 的源码审查中修正、以 40 点、容差 $10^{-12}$ 的单元测试固化（`tests/test_v22_fixes.py`）。这一细节在文献中极少被讨论，但对复现结果的影响是决定性的，故在此单列说明。

## 2.6 实验方法论与统计检验规范

多算法比较若只报均值，很容易把随机波动误读为算法优劣。机器学习领域的比较方法论由 Demšar 系统整理 [13]，本课题完整采用其规范：

1. **两算法比较**：当两算法的重复实验使用**完全对齐的随机种子**时，使用配对 Wilcoxon 符号秩检验 [10]（双侧）；若种子无法配对，退化为 Wilcoxon 秩和检验。选择非参数检验的原因是多目标优化指标（如 MIGD）的分布通常不满足正态性假设。
2. **多重比较校正**：同一组实验中把 LDT 与 5 个基线各做一次检验，共 5 次比较，必须控制族错误率（FWER）。本课题采用 **Holm 逐步下降校正** [11]：把 $m$ 个原始 $p$ 值升序排列为 $p_{(0)} \le \cdots \le p_{(m-1)}$，则
   $$\tilde p_{(r)} = \min\left(1,\ \max_{j \le r} (m-j)\, p_{(j)}\right)$$
3. **多算法多问题整体检验**：对每个数据集构造"重复 × 算法"的秩矩阵，先用 **Friedman 检验** [12] 判断"所有算法等价"这一原假设，若拒绝再用 **Nemenyi 事后检验** [13] 定位差异来源。Nemenyi 的临界差为
   $$CD = q_\alpha \sqrt{\frac{k(k+1)}{6N}}$$
   其中 $k$ 为算法数、$N$ 为数据集数、$q_\alpha$ 按 Demšar Table 5 查表。对 $k=6$、$N=30$、$\alpha=0.05$，$q_{0.05}=2.850$，故 $CD = 2.850\sqrt{42/180} \approx 1.377$——平均秩差小于 $CD$ 的两个算法在整体排序意义上不可区分。
4. **报告规范**：同时给出均值、标准差、中位数与四分位距；不使用"只报最优一格"的挑选式报告；对不显著的结果如实写出 $p$ 值而非隐去。

需要指出的是，**逐问题的 Wilcoxon 检验与跨问题的 Friedman–Nemenyi 检验回答的是不同问题**：前者问"在某一个（问题，$\tau$）配置下两算法是否可区分"，后者问"在若干配置的整体排序意义上两算法是否可区分"。前者的统计功效更高、更严格，后者更接近"综合表现是否相当"的工程判读。本课题两类检验同时报告，并在 5.5 节与 5.9 节明确说明两者的宽严差异——这是本文如实呈现结论的重要技术前提。

## 2.7 研究现状评述与研究缺口

综合本章文献，可得如下认识：

1. **响应式方法**（DNSGA-II 类 [3]）简单鲁棒、无预测风险，但缺乏方向信息，恢复慢；
2. **预测式方法**（PPS 类 [4]）在规则变化下收敛最快，但依赖历史样本与模型拟合，存在开销大、假设失效即偏差累积的风险；
3. **分布/知识迁移式方法**（KTMM 类 [14]）表达力强、对个体对应关系不敏感，但概率模型的参数估计成本高、可解释性与可审查性较弱；
4. **实验方法论**已相当成熟 [13]，但大量文献仍存在"只报均值、不做配对检验""参考前沿点数与口径不一致""复现信息不足"等问题。

由此提炼出本课题的三个研究缺口：

- **G1（轻量化缺口）**：是否存在一种只依赖**相邻（或少量）环境**、不需要训练时序模型或概率模型、单次响应计算量可忽略的迁移机制？
- **G2（可审查性缺口）**：在"黑箱模型 vs 简单规则"的取舍中，简化到何种程度仍能保留预测/迁移的收益？各组件各自的边际贡献能否被独立量化？
- **G3（可复现性缺口）**：能否建立一套口径统一（问题定义逐式可核、时间模型与预算一致、指标定义写死在代码中）、带单元测试与规范统计检验的实验系统，使上述结论可被第三方复现与证伪？

第 3 章的方法设计与第 4 章的实验系统分别针对 G1/G2 与 G3 给出答案。


\newpage

# 第 3 章 LDT-NSGA-II 算法与测试问题设计

本章给出算法的完整设计。**所有公式、常量与分支条件均与项目源码逐行对应**，读者可按文件与行号直接核对；标识为"代码位置"的条目均已在 v2.5 阶段由独立代码审计复核（见 5.8 节）。

## 3.1 设计目标与总体框架

### 3.1.1 设计目标

针对 2.7 节提出的三个研究缺口，LDT-NSGA-II（Lightweight Distribution-Transfer NSGA-II，轻量级分布迁移 NSGA-II）的设计目标如下：

- **T1（轻量）**：单次环境响应的计算量应与种群规模 $N$、决策维数 $D$ 成线性关系，不引入概率模型拟合、矩阵求逆或时序模型训练；
- **T2（仅依赖相邻环境）**：预测与迁移只使用"当前环境末种群 + 上一环境中心"这一最小历史窗口，不要求保存多个完整历史种群；
- **T3（可逐行审查）**：每个环节（特征提取、中心预测、采样、重组）都是独立可测的确定性函数，可被单元测试锁定；
- **T4（对预测失准稳健）**：必须保留一部分随机个体作为多样性兜底，使分布预测错误时不至于"全种群被锁死在错误区域"。

### 3.1.2 与 NSGA-II 的关系：策略与内核分离

LDT-NSGA-II **不是**对 NSGA-II 选择机制的修改，而是把 DMOP 求解显式拆成两个正交部分：

| 部分 | 职责 | 本课题实现 |
|---|---|---|
| **内核**（静态优化器） | 在**单一固定环境**内逼近 $PS(t_k)$、$PF(t_k)$ | 未改动的 pymoo NSGA-II（SBX/PM/快速非支配排序/拥挤度/精英保留） |
| **策略**（环境响应器） | 环境切换时生成新环境的**初始种群**（warm-start） | LDT 分布迁移模块，接口 `respond(prev_X, prev_F, problem, env_index, rng)` |

这一拆分带来两个重要的方法论后果：

1. **公平性**：所有对比算法共用同一内核、同一遗传算子、同一探针集、同一预算。算法之间**唯一的差异就是 `respond` 的实现**，因此性能差异可以归因到响应策略本身；
2. **可比性**：迁移种群不是"运行中注入"，而是作为新环境 NSGA-II 的**初始采样**。每段预算严格等于 $\tau_t \times N$ 次目标函数评价，不存在"某些算法因为注入个体而多算了评估次数"的不公平。

### 3.1.3 环境分段驱动器

系统运行骨架是环境分段驱动器 `DynamicDriver`（代码位置：`algorithms/nsga2_core.py`，第 65–183 行）：

1. 把 $n_{gen}$ 代划分为 $n_{envs} = \max(1, \lfloor n_{gen}/\tau_t \rfloor)$ 个环境段（默认 $300/25 = 12$ 段），每段恰好 $\tau_t$ 代；
2. 第 $k$ 段开始前调用 `problem.set_generation(k * tau_t)` 执行时间阶跃；
3. 当 $k > 0$ 时调用 `algorithm.respond(prev_X, prev_F, problem, k, rng)` 得到 warm-start 初始种群，随后为该段**新建一个全新的 pymoo NSGA2 实例**（`make_nsga2()`），以 warm-start 种群为 `sampling`，`seed = run_seed + k`，跑满 $\tau_t$ 代；
4. 逐代记录 IGD / HV / 各类耗时 / 累计评估数等 18 个字段；段末种群作为下一段的 `prev_X` / `prev_F`。

"每段新建 NSGA-II 实例"这一设计值得说明：它保证环境切换后父代种群不携带上一环境的内部状态（如已用代数、算子自适应状态），使"每段预算完全相同"这一公平性条件成立；而 `seed = run_seed + k` 又保证段与段之间不因同种子而产生相同的随机序列。

![图 3-1 LDT-NSGA-II 算法框架](../../LDT-NSGA-II/results/figures_thesis/fig3_1_framework.png)

图 3-1：LDT-NSGA-II 算法框架。主流程沿用 NSGA-II 内核；虚线框为本文新增的**分布迁移模块**（环境变化触发：非支配筛选 → 归一化特征提取 → 趋势中心预测 → 高斯采样迁移 → 50/20/30 种群重组），迁移结果作为新环境 NSGA-II 的初始种群。

## 3.2 归一化坐标系与分布特征提取

### 3.2.1 为什么必须固定坐标系

设第 $i$ 维决策变量的真实边界为 $x_l^{(i)}$ 与 $x_u^{(i)}$。DMOP 基准中不同维度的尺度可以差异很大：FDA1 的第一维在 $[0,1]$ 上而其余维在 $[-1,1]$ 上；dMOP1 的第一维在 $[0,1]$、其余维在 $[-1,1]$。若直接在真实空间计算"逐维标准差"，则宽度为 1 的维度与宽度为 2 的维度会被赋予相差一倍的重要性；更严重的是，若把"归一化后的均值"与"真实边界"混用于后续采样，采样结果的物理含义就是错的。

因此本算法的**第一条设计约束**是：**全部分布计算（特征提取、中心预测、融合、协方差缩放、采样、截断）都在归一化坐标系 $[0,1]^D$ 中完成，反归一化只执行一次**。

### 3.2.2 归一化与非支配筛选

设第 $k-1$ 环境结束时的种群为 $X \in \mathbb{R}^{N \times D}$，逐维归一化（代码位置：`algorithms/distribution_features.py::normalize`）：

$$Z = \frac{X - x_l^{\mathsf T}}{x_u^{\mathsf T} - x_l^{\mathsf T}} \in [0,1]^{N \times D}$$

随后按目标空间的 Pareto 支配关系（判定带 $10^{-8}$ 容差）取出非支配子集 $ND \subseteq Z$，即

$$ND = \{z \in Z \mid \nexists z' \in Z:\ F(z') \prec F(z)\}$$

**回退规则**：若 $|ND| < \max(10, N/4)$（代码位置：`ldmt_nsga2.py` 第 145–146 行），说明非支配层过小、其矩特征不能代表种群分布，此时**退回使用全种群** $Z$。对默认 $N = 100$，该阈值为 $25$。该规则的设计意图是：既要过滤掉被支配的差解，又不能在非支配层退化时把特征建立在个别解上。

### 3.2.3 逐维矩特征

对 $ND$ 计算逐维均值和逐维总体标准差（`ddof=0`）（`distribution_features.py::extract_features`）：

$$\mu_e = \frac{1}{|ND|}\sum_{z \in ND} z \in \mathbb{R}^D, \qquad
\sigma_e = \sqrt{\frac{1}{|ND|}\sum_{z \in ND}(z - \mu_e)^2} \in \mathbb{R}^D$$

**标准差下限**：$\sigma_e \leftarrow \max(\sigma_e,\ \sigma_{\min})$，其中 $\sigma_{\min} = 10^{-3}$（归一化单位，约等于跨度的 0.1%）（`ldmt_nsga2.py` 第 149 行）。设此下限的原因是防止**协方差塌缩**：当种群已经高度收敛（$\sigma_e \to 0$）时，高斯采样会退化为"重复采样中心点"，迁移部分完全丧失多样性；一旦预测中心有偏，整个迁移子群将同时落在错误位置。

$\mu_e$ 与 $\sigma_e$ 构成所谓的**分布特征**——前者刻画"种群中心在哪"，后者刻画"种群散布多宽"。这就是本文所说"轻量级"的具体含义：分布被压缩为 $2D$ 个标量，而完整高斯混合模型需要 $O(K D^2)$ 个参数。

## 3.3 趋势中心预测

### 3.3.1 半步一阶差分外推

设上一环境（第 $k-2$ 段末）的归一化中心为 $\mu_{prev}$，则第 $k-1$ 环境的预测中心为（`ldmt_nsga2.py` 第 159–164 行）：

$$\hat\mu = \mathrm{clip}\big(\mu_{\text{eff}} + \alpha\,(\mu_{\text{eff}} - \mu_{prev}),\ 0,\ 1\big), \qquad \alpha = 0.5 \text{（默认 \texttt{trend\_alpha}）}$$

其中 $\mu_{\text{eff}}$ 是融合后的当前中心（3.5.2 节）：当多源融合关闭（默认）时 $\mu_{\text{eff}} \equiv \mu_e$，故上式退化为"用当前中心减历史中心"。其含义是**沿"最近一次中心位移"的方向外推半步**。若种群中心做匀速平移，即 $\mu_{k} = \mu_{k-1} + \delta$，则真实的新中心为 $\mu_{k-1} + \delta = \mu_e + \delta$，而 $\mu_e - \mu_{prev} = \delta$，故 $\alpha = 1$ 恰好命中。默认取 $\alpha = 0.5$ 是**保守折中**：在外推不足（$\alpha$ 偏小）与过冲（$\alpha$ 偏大）之间取中。第 5 章的敏感性实验将如实显示这一保守选择在本组基准上并非最优（5.6.4 节），v2.5 的对照实验则进一步把它与"迁移范式"的贡献分离开来（5.7 节）。

**首次变化无历史**：当 $\mu_{prev}$ 不存在（即 $k = 1$）时，取 $\hat\mu = \mathrm{clip}(\mu_e, 0, 1)$，退化为"直接沿用上一环境中心"。

### 3.3.2 截断的必要性

`clip(·, 0, 1)` 不是可有可无的数值细节：由于外推方向由历史位移决定，当 PS 运动方向发生反转或到达边界附近时，外推结果可能落在 $[0,1]^D$ 之外。截断保证迁移候选始终是可行域内的点，从而不会在反归一化后产生越界解。

### 3.3.3 记录供下一步使用的特征

第 $k$ 环境响应完成后，算法把自己这一环境的**未融合**特征 $\mu_e$ 保存为下一轮的 $\mu_{prev}$（`ldmt_nsga2.py` 第 179 行，保存的是 `mu_e.copy()` 而非 `mu_eff`）。这一点在多源融合启用时会造成**趋势项口径不一致**：趋势项的第一个量是融合后的 $\mu_{\text{eff}}$，而减去的 $\mu_{prev}$ 是未融合的历史中心。由于 $\mu_{\text{eff}} - \mu_{prev} = (\mu_{\text{eff}} - \mu_e) + (\mu_e - \mu_{prev})$，两者之差恰为融合校正项，因此该写法等价于"在趋势项上叠加一个融合引起的位移"，并**不构成数值错误**；但它是否为设计意图，作者**不确定**——这是本课题如实披露的一处设计不确定（见 5.8 节 K13）。默认实验未开启多源融合，故不影响第 5 章的全部结论。

## 3.4 种群重组：50% 迁移 + 20% 保留 + 30% 随机移民

新环境的初始种群由三部分拼接而成，比例参数在构造时被断言约束为 $\rho + \kappa + r = 1$（`ldmt_nsga2.py` 第 47–51 行），默认 $\rho = 0.5,\ \kappa = 0.2,\ r = 0.3$。各部分个体数取 $n_{\text{part}} = \mathrm{round}(\text{ratio} \cdot N)$。

### 3.4.1 迁移个体（$\rho = 0.5$）

在归一化空间以 $\hat\mu$ 为中心、$\sigma_e$ 为逐维尺度做对角高斯采样（`ldmt_nsga2.py` 第 187–197 行）：

$$z \sim \mathcal{N}(0, I_D), \qquad x_{\text{norm}} = \mathrm{clip}\big(\hat\mu + z \odot \sigma_e,\ 0,\ 1\big)$$

再**反归一化恰好一次**映射回真实决策空间：

$$x = x_l^{\mathsf T} + x_{\text{norm}} \odot \big(x_u^{\mathsf T} - x_l^{\mathsf T}\big)$$

若出现非有限值行，则以真实边界上的均匀随机解替换（防御性分支，正常路径不会触发）。

这一部分是**算法性能的主要来源**：它把"上一环境的最优解分布"按预测方向整体搬移到新环境，同时保留了一阶矩（中心）与二阶矩（散布）信息。注意协方差被假设为对角（逐维独立），这是一个有意的简化——完整协方差需要 $O(D^2)$ 参数并在小样本下易出现奇异，而基准问题的 PS 运动本身近似沿坐标轴方向的整体平移。

### 3.4.2 保留个体（$\kappa = 0.2$）

从上一种群中随机抽取 $\kappa N$ 个个体（`rng.permutation`，`ldmt_nsga2.py` 第 199–205 行）。其设计意图是保留上一环境已经验证过的好解：若新环境变化很小，这些个体本身就已接近新 $PS$。

第 5 章的消融实验将如实报告：**在本组基准上该组件无可检测增益**（5.6 节）。

### 3.4.3 随机移民（$r = 0.3$）

在真实边界内均匀采样 $rN$ 个个体（`ldmt_nsga2.py` 第 207–212 行）：

$$x \sim \mathcal{U}(x_l, x_u)$$

这是对应设计目标 T4 的**多样性兜底**：当方向预测完全错误时，随机移民保证种群中始终存在覆盖整个决策空间的个体，使 NSGA-II 仍有能力在段内找到正确区域，而不会全种群被锁死在错误位置。

### 3.4.4 定容与配比

三部分拼接后由 `_resize` 精确补齐或截断到 $N$（`ldmt_nsga2.py` 第 228–235 行）：若总行数超过 $N$ 则截断**前 $N$ 行**（顺序为迁移 → 保留 → 随机）；若不足则用真实边界上的均匀随机行**逐行追加**补齐。由于 $\mathrm{round}$ 的舍入误差，三部分之和可能差 1–2 行，`_resize` 负责吸收这一误差。

![图 3-2 分布迁移策略示意](../../LDT-NSGA-II/results/figures_thesis/fig3_2_transfer.png)

图 3-2：分布迁移策略示意（归一化决策空间二维投影）。上一环境种群以均值 $\mu_e$ 与协方差椭圆刻画；趋势中心预测 $\hat\mu=\mu_e+\alpha(\mu_e-\mu_{prev})$ 沿中心增量方向外推（默认半步 $\alpha=0.5$）；迁移个体由 $\hat\mu$ 为中心、$\sigma_e$ 为宽度的高斯采样生成，再按 50/20/30 与保留个体、随机移民重组。

## 3.5 两个可选机制（默认关闭）

### 3.5.1 自适应迁移宽度缩放

设计直觉是"环境变化越剧烈，迁移分布应撒得越宽"。用共享探针的变化幅度 $D_k$（定义见 3.6 节）作为信号，并以运行至今的累计均值 $\overline{|D|}$ 归一化（`ldmt_nsga2.py` 第 166–176 行）：

$$\text{scale}_k = \mathrm{clip}\left(1 + \beta \frac{|D_k|}{\overline{|D|}},\ 0.5,\ 2.0\right), \qquad \sigma_{\text{eff}} = \text{scale}_k \cdot \sigma_e$$

其中 $\beta$ 为 `adaptive_beta`（默认 1.0），上下限为类常量 `SIGMA_SCALE_MIN = 0.5` 与 `SIGMA_SCALE_MAX = 2.0`，$\overline{|D|}$ 为**包含当前 $|D_k|$ 在内的运行累计均值**。归一化使缩放量对问题的绝对目标尺度不敏感；上下限则防止剧烈变化时分布"爆炸"（上限）或因 $D_k \approx 0$ 而分布塌缩（下限）。

严格地说，这是对 $\sigma$ 的**线性缩放**（方差按 $\text{scale}_k^2$ 放大），而非对完整协方差矩阵的缩放——因为算法本身只维护逐维标准差（对角近似）。该机制由 `adaptive_transfer=True` 启用，**默认关闭**，第 7 章将其列为待补实验之一。

### 3.5.2 多源历史融合

设计直觉是"单步差分对噪声敏感，多个历史环境的分布可以互相印证"。算法维护最近 `history_size = 3` 个环境的 $(\mu_i, \sigma_i)$ 历史（`ldmt_nsga2.py` 第 96–125 行，第 181–184 行）：

**相似度权重**：以当前环境的中心 $\mu_e$ 为参照，

$$w_i = \exp\left(-\frac{\|\mu_e - \mu_i\|_2}{\sqrt{D}}\right), \qquad D = \dim(\mu)$$

归一化因子 $\sqrt{D}$ 使权重对维数不敏感。权重低于 `FUSION_WEIGHT_FLOOR = 0.05` 的历史源被**剔除**（视为与当前环境分布差异过大、不宜迁移）；当前环境自身的权重固定为 1.0，然后与保留的历史权重一并归一化。

**混合分布的均值与二阶矩方差**：设归一化后权重为 $\tilde w_i$，

$$\mu_f = \sum_i \tilde w_i \mu_i, \qquad
\sigma_f^2 = \sum_i \tilde w_i\left(\sigma_i^2 + \|\mu_i - \mu_f\|^2\right)$$

第二个式子中的 $\|\mu_i - \mu_f\|^2$ 项是**混合分布二阶矩**的正确形式：多个来源的中心不重合时，总方差必须计入"中心之间的散布"，而不能只对标准差做朴素加权平均。这是 v2.2 源码审查中修正的一处缺陷（K4 相关）。该机制由 `multi_source=True` 启用，**默认关闭**。

### 3.5.3 可插拔策略层

`algorithms/transfer_strategy.py` 定义了 `TransferStrategy` 接口与两个实现（`MomentTransferStrategy`、`GMMTransferStrategy`），使 `respond` 的整体行为可以被替换而不改动主算法（`LDTNSGA2(transfer_strategy=...)`）。该层在 v2.4 修复了归一化/真实坐标混用问题；**默认实验未使用该层**，它只提供扩展点。

## 3.6 共享固定探针与变化幅度度量

### 3.6.1 共享固定探针

每个问题实例在构造时携带一组**固定**探针点（代码位置：`problems/base.py` 第 20–36 行）：`PROBE_SIZE = 10` 个点、由 `PROBE_SEED = 2026` 的独立随机流在**真实边界** $[x_l, x_u]$ 内均匀抽取，逐维独立且与该问题实例的维数一致。

三个设计要点：

1. **固定**：探针不随运行、种子、算法变化。同一问题下所有算法、所有种子看到的是同一组探针点，因此 $D_k$ 在算法之间可直接比较；
2. **真实边界**：探针在真实决策空间抽取，而不是"把种群行裁剪进 $[0,1]$"。否则对于边界为 $[-1,1]$ 的维度，$D_k$ 会被裁剪操作畸变；
3. **共享**：由问题对象持有（`problem.probe`），算法通过 `_ensure_probe` 取用，避免各算法各自造探针。

### 3.6.2 变化幅度 $D_k$

`ChangeDetector.detect`（`algorithms/environment_manager.py`）按如下公式计算相邻环境在探针点上的目标值平均相对偏移。设 $F^{\text{prev}}, F^{\text{new}} \in \mathbb{R}^{P \times M}$ 分别为上一环境与当前环境在 $P$ 个探针点上的目标值矩阵，则先逐目标维取归一化尺度

$$s_j = \max_{p} \big|F^{\text{prev}}_{pj}\big|, \qquad s_j \leftarrow 1 \text{ 若 } s_j \le 10^{-12}$$

再对全部 $P \times M$ 个元素取相对偏移的算术平均：

$$D_k = \frac{1}{PM}\sum_{p=1}^{P}\sum_{j=1}^{M}\frac{\big|F^{\text{new}}_{pj} - F^{\text{prev}}_{pj}\big|}{s_j}$$

即"新环境对探针点的目标值相对上一环境的平均相对偏移"。逐目标归一化（除以该目标上一环境的最大绝对值）使量纲与尺度差异很大的目标可以相加；$s_j \le 10^{-12}$ 时的兜底赋值避免除零。

**必须强调**：$D_k$ **不是变化触发器**。何时变化由驱动器的周期 $\tau_t$ 硬性决定（2.3.2 节的时间模型），算法不需要、也不被允许自行判定"是否发生变化"；$D_k$ 只作为**诊断量**（写入 `respond` 返回的 `info`）与自适应信号（3.5.1 节）。首次变化时无上一环境参照，$D_k = 0$。项目中的 `EnvironmentManager.is_change(g)` 与 `_detect_change` 的阈值判定在生产路径中均未被调用（仅被单元测试引用）。

| 项目 | 取值 | 代码位置 |
|---|---|---|
| 探针数 $P$ | 10 | `problems/base.py::PROBE_SIZE` |
| 探针随机种子 | 2026 | `problems/base.py::PROBE_SEED` |
| 探针分布 | 真实边界 $[x_l,x_u]$ 上均匀 | `problems/base.py` 第 33–36 行 |
| $\sigma$ 下限 | $10^{-3}$（归一化单位） | `ldmt_nsga2.py` 第 149 行 |
| 趋势系数 $\alpha$ | 0.5 | `LDTNSGA2.__init__` |
| 迁移/保留/随机比 | 0.5 / 0.2 / 0.3 | `LDTNSGA2.__init__` |
| 非支配子集回退阈值 | $\max(10, N/4)$（$N=100$ 时为 25） | `ldmt_nsga2.py` 第 145 行 |
| 自适应缩放区间 | $[0.5, 2.0]$ | `SIGMA_SCALE_MIN/MAX` |
| 自适应系数 $\beta$ | 1.0 | `adaptive_beta` |
| 融合历史容量 | 3 | `history_size` |
| 融合权重下限 | 0.05 | `FUSION_WEIGHT_FLOOR` |

表 3-1：LDT-NSGA-II 的全部数值常量（与源码逐项对应）。

## 3.7 计算复杂度与"轻量"的量化

设种群规模 $N$、决策维数 $D$、目标数 $M$、历史容量 $H = 3$。

| 环节 | 时间复杂度 | 说明 |
|---|---|---|
| 归一化 / 反归一化 | $O(ND)$ | 逐元素运算 |
| 非支配筛选 | $O(MN^2)$ | 逐一比较；与 NSGA-II 内部排序同阶，不增加数量级 |
| 矩特征提取 | $O(ND)$ | 均值与标准差各一遍扫描 |
| 趋势中心预测 | $O(D)$ | 一次向量加法与截断 |
| 多源融合（可选） | $O(HD)$ | 权重、均值、二阶矩各一遍 |
| 高斯采样迁移个体 | $O(\rho N D)$ | 一次标准正态采样 |
| 拼接与定容 | $O(ND)$ | |
| **单次响应合计** | $\mathbf{O(MN^2 + ND)}$ | 与 $N$、$D$ 线性相关的部分占主导；无矩阵分解、无迭代优化 |

对比之下，GMM 型分布迁移需要在每个环境切换时对非支配子集拟合高斯混合模型（期望最大化迭代，每次迭代 $O(NKD^2)$），其单次响应的计算量与内存占用都显著更高。第 5 章将以实测的 `migration_runtime_s` 给出两者的比例关系（5.5.6 节）。

## 3.8 五个对比基线的设计

全部基线与 LDT 共用同一驱动器、同一 NSGA-II 内核、同一算子配置、同一探针集与同一预算，**仅在 `respond` 上不同**。下表如实描述其真实行为（v2.5 由独立审计复核）。

| 算法标识 | 类型 | 环境变化时的 `respond` 行为 | 迁移耗时量级 | 出处 |
|---|---|---|---|---|
| `nsga2_static` | 静态基线 | 原样带入上一种群（`prev_X.copy()`），**无任何响应** | $10^{-4}$ s | NSGA-II [1] |
| `dnsga2_inspired` | 响应式 | 保留种群**前缀** $0.8N$ 行 + 注入 $0.2N$ 个均匀随机移民 | $10^{-4}\sim10^{-3}$ s | DNSGA-II [3]（项目实现，非原文复现） |
| `ar1_predictive` | 预测式 | 对非支配子集中心做逐维 OLS 一阶自回归外推（历史容量 5），预测个体整体平移 + **20%** 均匀随机移民；**方差不预测** | $10^{-2}$ s | 中心外推一族 [4] |
| `gmm_inspired` | 分布迁移 | 对非支配子集拟合高斯混合（`n_components = min(3, max(1, N/10))`，`reg_covar = 1e-4`，协方差类型依次尝试 diag→full→spherical），从混合模型采样 $N$ 行后保留前 $N - n_{\text{imm}}$ 行 + **10%** 移民 | $10^{-1}$ s | 分布迁移一族 [14]（项目实现，非 KTMM 复现） |
| `pps` | 预测式 | 中心历史容量 3，$\hat\mu = \mathrm{clip}(2\mu_{k-1} - \mu_{k-2}, 0, 1)$（**全步**一阶差分），种群整体平移 + **20%** 均匀随机移民；无独立方差预测项 | $10^{-2}$ s | PPS [4]（项目简化实现，非完整版） |
| **`ldt_nsga2`（本文）** | 迁移 + 预测 | 半步趋势外推 + 高斯采样迁移 + 50/20/30 三路重组 | $10^{-2}$ s | 本文 |

表 3-2：六个算法的 `respond` 行为对照。

需要明确三处**命名与归属的诚实声明**，以及一处**参数不对称**：

1. `dnsga2_inspired` 保留的是种群**前缀**而非随机抽样（Deb 原文为随机替换），因此它不是 DNSGA-II 的逐字复现，只沿用"注入随机移民恢复多样性"的思想；
2. `gmm_inspired` 是"用高斯混合模型做分布迁移"的项目实现，**不是** KTMM（或同类混合模型知识迁移方法）[14] 的复现；
3. `pps` 是 PPS [4] 的**简化**项目实现：只保留"中心一阶差分全步外推 + 整体平移 + 随机移民"，**没有** PPS 原文中"中心预测 + 形/方差预测"的完整结构；
4. **移民比例不一致**：`ar1_predictive` 与 `pps` 的默认 `immigrant_ratio = 0.2`，而 `gmm_inspired` 的默认值为 **0.1**——即 GMM 基线的随机移民比例只有另两个预测基线的一半。这是源码中的既有事实（`gmm_transfer_nsga2.py` 第 21 行），本文如实报告而不做统一，因为在同一实验中同时改变移民比例会使基线偏离其原始实现更远；但它意味着从第 5 章读出的 GMM 性能里，混杂了"移民更少"这一因素。

这些声明直接决定第 5 章结论的适用边界：本文报告的是"轻量矩迁移 vs 三种简化策略实现"的受控比较，**不能**被解读为对 PPS、DNSGA-II 或任何混合模型方法原始版本的性能评定。

## 3.9 测试问题集合

### 3.9.1 FDA 系列（与 PlatEMO 逐式核对）

FDA1、FDA3、FDA4 的实现逐式对照 PlatEMO master 的 MATLAB 源码，核对方式为抓取 `PlatEMO/Problems/Multi-objective optimization/FDA/` 下的源文件逐式比对。

**FDA1**（`problems/fda1.py`；$M=2$，$D=10$；$x_1 \in [0,1]$，$x_i \in [-1,1]$）：

$$G(t) = \sin\!\left(\frac{\pi t}{2}\right),\quad
g = 1 + \sum_{i=2}^{10}(x_i - G)^2,\quad
f_1 = x_1,\quad
f_2 = g\left(1 - \sqrt{f_1/g}\,\right)$$

$PS(t)$：$x_i = G(t)\ (i \ge 2)$，随 $t$ 平移；$PF = \{(u, 1-\sqrt u) : u \in [0,1]\}$ **时不变**。属 **Type I**。

**FDA3**（`problems/fda3.py`；$M=2$，$D=10$；边界同 FDA1）：

$$k(t) = 10^{\,2\sin(\pi t/2)},\quad
f_1 = x_1^{\,k(t)},\quad
G = \left|\sin\!\left(\frac{\pi t}{2}\right)\right|,$$
$$g = 1 + G + \sum_{i=2}^{10}(x_i - G)^2,\quad
f_2 = g\left(1 - \sqrt{f_1/g}\,\right)$$

指数 $k(t)$ 只对 $f_1$ 做单调重参数化，改变前沿上解的分布密度而不改变前沿集合本身；$PF(t)$ 随 $t$ 变化，属 **Type II**。

**FDA4**（`problems/fda4.py`；$M=3$，$D = M+9 = 12$；全部 $x_i \in [0,1]$）。前 $M-1$ 个变量为角度变量，**其余全部**（0 基下标 $M-1$ 起，共 $D-(M-1)=10$ 个）为距离变量：

$$g = \sum_{i=M}^{D}(x_i - G)^2, \qquad \theta_j = x_j \cdot \frac{\pi}{2}\ (j = 1,\ldots,M-1),$$
$$f_1 = (1+g)\cos\theta_1\cos\theta_2,\quad
f_2 = (1+g)\cos\theta_1\sin\theta_2,\quad
f_3 = (1+g)\sin\theta_1$$

$PS(t)$：距离变量 $= G(t)$；$PF$ 为第一卦限单位球面上一段，**时不变**，属 **Type I**。PlatEMO 中距离变量切片写作 `PopDec(:, M:end)`，对应 Python 的 `X[:, M-1:]`——这一索引对应关系是本课题 v2.2 修正的关键缺陷（2.5.3 节）。

### 3.9.2 权威 dMOP 系列（v2.5 新增，按原始文献实现）

dMOP1–dMOP3 源于 Goh 与 Tan 的共进化范式研究 [5]；其定义在文献中的复现版本另见 Helbig 与 Engelbrecht 的章节 [18] 等（详见本节末的可追溯性声明）。**PlatEMO 不包含 dMOP1–dMOP3 的实现**（本课题经 GitHub git-tree API 对 master 全树递归检索确认，见 5.11 节 K11），因此实现必须直接依据文献。v2.5 新增以下三个问题类，它们覆盖了本文测试集中**全部三种动态类型**：Type I（DMOP3A）、Type II（DMOP2A）、Type III（DMOP1A）。

> **坐标边界声明（易错点）**：权威 dMOP1/2/3 的**全部**决策变量都在 $[0,1]$ 内（$D=10$），这与 FDA 系列的约定**不同**——FDA1/FDA3 的 $x_1 \in [0,1]$ 而 $x_{i\ge2} \in [-1,1]$。本文早期草稿曾误用 FDA 的边界，v2.5 已更正。把 FDA 的边界套到 dMOP 上会改变可行域、进而改变 Pareto 集与实验结果。

**DMOP1A**（`problems/dmop1a.py`；$M=2$，$D=10$；全部 $x_i \in [0,1]$）：

$$H(t) = 0.75\sin\!\left(\frac{\pi t}{2}\right) + 1.25 \in [0.5,\ 2.0],$$
$$g = 1 + 9\sum_{i=2}^{10} x_i^2 \quad(\text{注意：不含 } G(t)\text{ 项}),\qquad
f_1 = x_1,\qquad
f_2 = g\left(1 - \left(\frac{f_1}{g}\right)^{H(t)}\right)$$

由 $g \ge 1$ 且等号仅在 $x_i = 0\ (i \ge 2)$ 时成立，得

$$PS = \{(x_1, 0, \ldots, 0) : x_1 \in [0,1]\}\quad\text{（与 } t \text{ 无关）}$$

而在 $PS$ 上 $g = 1$、$f_1 = x_1$，故

$$PF(t) = \left\{\left(u,\ 1 - u^{H(t)}\right) : u \in [0,1]\right\}\quad\text{（随 } t \text{ 变化）}$$

即**解集静止、前沿随时间旋转**（$H>1$ 凸、$H=1$ 退化为直线 $1-u$、$H<1$ 凹）。这是 **Type III** 问题，也是本文所用六个测试问题中**唯一"解集已经正确、而前沿形状仍在变"的配置**——它对算法的要求与 FDA1 恰好相反：不需要移动解集，而需要重新调整解集在前沿上的分布密度。

> **重要更正说明（v2.5）**：本课题早期草稿曾把权威 dMOP1 记为 **Type I**、$g$ 含 $G(t)$、前沿为 $\sqrt{\cdot}$ 型。经文献复核，该写法**错误**——它实际混合了 dMOP2 的 $g$ 与 dMOP3 的前沿。正确的权威 dMOP1 是上述 **Type III** 形式；本文按复核结果重新实现并在此显式更正。

**DMOP2A**（`problems/dmop2a.py`；$M=2$，$D=10$；全部 $x_i \in [0,1]$）：

$$G(t) = \sin\!\left(\frac{\pi t}{2}\right),\qquad
H(t) = 0.75\sin\!\left(\frac{\pi t}{2}\right) + 1.25 \in [0.5,\ 2.0],$$
$$g = 1 + 9\sum_{i=2}^{10}(x_i - G)^2,\quad
f_1 = x_1,\quad
f_2 = g\left(1 - \left(\frac{f_1}{g}\right)^{H(t)}\right)$$

$PS(t)$：$x_i = G(t)$（随 $t$ 平移）；在 $PS$ 上 $g = 1$ 且 $f_1 = x_1 = u$，故

$$PF(t) = \left\{\left(u,\ 1 - u^{H(t)}\right) : u \in [0,1]\right\}$$

前沿形状随 $H(t)$ 变化：$H = 2$ 时退化为凸抛物线 $1 - u^2$，$H \to 0.5$ 时向 $1-\sqrt u$ 弯曲。**$PS$ 与 $PF$ 同时随时间变化**，属 **Type II**。注意 dMOP1A 与 dMOP2A 共享同一条前沿族 $1 - u^{H(t)}$，区别只在解集是否移动——这使二者构成一对干净的对照：同样是"前沿旋转"，一个解集静止、一个解集平移。

**DMOP3A**（`problems/dmop3a.py`；$M=2$，$D=10$；全部 $x_i \in [0,1]$）：

$$G(t) = \left|\sin\!\left(\frac{\pi t}{2}\right)\right|,\qquad
r \in \{1, \ldots, D\}\ \text{（"位置相关"变量下标，按环境轮换）},$$
$$g = 1 + 9\sum_{i \ne r}(x_i - G(t))^2,\qquad
f_1 = x_r,\qquad
f_2 = g\left(1 - \sqrt{f_1/g}\,\right)$$

$PS(t) = \{x : x_r \in [0,1],\ x_i = G(t)\ (i \ne r)\}$——**解集随 $t$ 移动**；在 $PS$ 上 $g = 1$ 且 $f_1 = x_r$，故

$$PF = \left\{\left(u,\ 1 - \sqrt u\right) : u \in [0,1]\right\}\quad\text{（与 } t \text{ 无关）}$$

属 **Type I**。该问题的独特难点是**下标 $r$ 的切换**：当 $r$ 改变时，原本"位置相关"的变量变成"距离相关"、必须同时跳到 $G(t)$，会造成严重的多样性损失——这正是对 warm-start 机制的强压力测试。本实现按**环境**（而非按代）以轮转方式切换 $r$（`set_generation` 由 $\lfloor g/\tau_t\rfloor$ 推出 $r$），保证同一环境内目标函数不变、环境切换时 $r$ 恰好切换一次。

> **一处如实声明的文献分歧**：dMOP3 的 $G(t)$ 在两份复现文献中不一致——Helbig 与 Engelbrecht 的同行评议章节写作 $G(t)=\sin(\pi t/2)$，而 Brock 大学 2021 年的硕士论文写作 $G(t)=\lvert\sin(\pi t/2)\rvert$。由于原始论文无开放获取版本、无法据原文裁定，本实现采用**绝对值变体**（与 FDA3/FDA4/FDA5 的 $\lvert\sin(\cdot)\rvert$ 惯例一致，且该文献给出了自洽的 $PS$/$PF$ 推导）。两种选择下**问题的类型与 Pareto 前沿完全相同**，只有解集的运动轨迹不同。

**dMOP 定义的可追溯性声明**：Goh 与 Tan 2009 的原文受付费墙限制且无开放获取版本，本课题**未能直接阅读原始公式**。上述定义依据四份独立引用该文的复现文献交叉确认：Helbig 与 Engelbrecht 的 CSIR 同行评议章节（式 16–18）、Joćko 的 Brock 大学 2021 硕士论文（式 2.21–2.27）、Adekoya 的 Pretoria 大学 2019 博士论文（式 2.13）、MDPI *Mathematics* 2022, 10(15):1587 表 2–3。其中 MDPI 那一份把 dMOP2 的 $g$ 系数印为 1（其余三份均为 9），且该文在别处亦有明显笔误（如把 FDA2 的 $H(t)$ 印成 $0.75+0.7\sin$），故本课题采用**系数 9**并将该分歧如实记录于 `_audit/dmop_definitions.md`。

> **自检证据**：上述"把决策向量放在解析导出的 $PS$ 上，其目标值是否落在解析 $PF$ 上"的一致性已由脚本 `scripts/check_dmop_authoritative.py` 自动验证（$t \in \{0, 0.05, 0.2, 0.45, 0.8\}$、每个前沿 200 个采样点；DMOP3A 另对全部 10 个可能的 $r$ 逐一验证），最大偏差均为 **0.000e+00**。脚本同时验证了三项**问题类型**与文献一致：DMOP1A 的 $PS$ 静止（$\lvert\Delta PS\rvert = 0$）而 $PF$ 变化（$\lvert\Delta PF\rvert = 0.1205$）→ Type III；DMOP2A 的 $PS$ 与 $PF$ 都变化（$0.6494$ / $0.1205$）→ Type II；DMOP3A 的 $PS$ 变化（$0.6494$）而 $PF$ 静止（$0$）→ Type I；以及 dMOP3A 的 $r$ 在每个环境内恒定、跨环境共取到 10 个不同值。

**版本更正记录**：v2.5 期间，权威 dMOP1 与 dMOP3 的实现曾按错误的公式生成过一批实验数据（dMOP1 误作 Type I + $\sqrt{\cdot}$ 前沿、dMOP3 误用了项目变体的欺骗性前沿与 $H(t)$ 解集位置）。该批 12 个 CSV 已被移出统计目录、归档于 `results/raw_v25_badauthoritative/`，并在**修正后重新运行**。DMOP2A 的实现经复核为逐式正确，无需重跑。本论文第 5 章报告的全部 dMOP 数据均为修正后的版本。

### 3.9.3 三个项目变体（保留并如实标注）

早期的项目实现保留了 `DMOP1`、`DMOP2`、`DMOP3` 三个类（`problems/dmop1.py`、`dmop2.py`、`dmop3.py`），它们是**项目自定义变体**，与权威定义不同：

| 类 | 与权威定义的差异 | v2.5 后的定位 |
|---|---|---|
| `DMOP1` | $g$ 的系数为 1（非 9）；前沿指数为 2（非 $\sqrt{\cdot}$），得 $PF = \{(u, 1-u^2)\}$ | 项目变体，v2.4 补充实验使用；引用须注明"非标准 dMOP1" |
| `DMOP2` | $g$ 系数 1、含额外 $+G$ 项、前沿指数固定为 2（无 $H(t)$） | 项目变体；引用须注明"非标准 dMOP2" |
| `DMOP3` | 目标函数、$PS$、$PF$ **均与时间无关**（静态欺骗问题）；死代码 $G(t)$ 已于 v2.4 删除 | **静态对照组**（检验迁移机制在无变化场景是否有害），其上算法比较不具动态区分度 |

这三条差异是本课题 v2.4 源码审查的核心发现，被完整记入第 5.8 节的已知问题清单（K2、K3）。v2.5 新增权威实现后，第 5 章将**同时报告**两族问题的结果，并把差异的来源讲清楚。

### 3.9.4 参考前沿生成

参考前沿由 `problems/true_fronts.py` 生成，默认 500 点，按 (问题, $M$, 点数, $t$) 缓存：

- **二目标问题**：在 $u \sim \mathrm{linspace}(0,1)$ 上均匀取点代入解析式；其中前沿随 $t$ 变化的有 **FDA3、DMOP2（项目变体）以及 DMOP1A、DMOP2A（权威变体）**——注意 DMOP1A 属于"解集静止、前沿随时变化"的 Type III，其参考前沿同样必须按当前 $t$ 生成，这是实现中最容易漏掉的一处（若误用 $t=0$ 的前沿，IGD 会被系统性高估）；
- **FDA4（三目标）**：采用 $[0,\pi/2]$ 上的**确定性角度网格**（每维 23 个点，$23 \times 23 = 529$ 个组合）做球面参数化 $(\cos A\cos B,\ \cos A\sin B,\ \sin A)$，去重后实测 **507 点**。

必须如实声明两点：(a) 这是**项目自定义的确定性采样**，等角度间距而非等面积，与 PlatEMO 的 `UniformPoint` 规则不同；(b) 因此 FDA4 的 IGD 绝对值与"恰好 500 点"的表述略有出入。由于全部算法使用同一参考前沿，该差异不影响算法之间的比较，但会影响 IGD 的绝对数值与其他文献的横向对比——这一点在 5.1.3 节与 5.9 节再次说明。

## 3.10 本章小结

本章给出 LDT-NSGA-II 的完整设计：归一化空间中的逐维矩特征（均值 + 标准差 + $10^{-3}$ 下限）、半步一阶差分趋势外推（$\alpha = 0.5$）、50/20/30 的三路种群重组，以及自适应协方差缩放与多源二阶矩融合两个默认关闭的可选机制。单次环境响应的时间复杂度为 $O(MN^2 + ND)$，与种群规模和维数线性相关，不引入概率模型拟合。同时给出五个行为受控的对比基线、六个测试问题（FDA1/FDA3/FDA4 逐式核对 PlatEMO，DMOP1A/DMOP2A/DMOP3A 按原始文献实现并通过自检脚本验证），以及三个如实标注的项目变体。第 4 章转向实验系统：如何把这些设计变成一条可复现、可审计的实验流水线。


\newpage

# 第 4 章 实验系统与可复现性设计

本章说明实验是如何组织的，以及结论如何被第三方复现。设计原则是：**任何进入论文的数字都必须能由仓库内的脚本从原始数据重新算出来**，且口径（时间模型、预算、指标、统计方法）写在代码里而不是写在论文里。

## 4.1 软硬件环境与依赖

| 项 | 值 | 来源 |
|---|---|---|
| 操作系统 | Windows 10/11（开发与实验环境） | 运行元数据 `run_meta.jsonl` 的 `platform` 字段 |
| Python | 3.9.13 | 同上（`environment.python`） |
| 关键依赖 | pymoo 0.6.1.5、numpy 2.0.2、scipy 1.13.1、scikit-learn 1.6.1、pytest 8.4.2 | 实测 |
| 声明依赖 | `pymoo>=0.6`、`numpy`、`scipy`、`scikit-learn`、`matplotlib`、`pyyaml`、`tqdm`、`pytest`（`requirements.txt` 全文） | 仓库 |
| 并行度 | v2.5 两阶段实验按 (问题, $\tau$) 独立进程并行执行（最多 9 进程） | `run_v25_comparison.ps1` |

**版本兼容声明（必须知悉）**：pymoo 0.6.x **没有** `algo.n_eval` 属性，真实的评估计数器位于 `algo.evaluator.n_eval`。项目以兼容封装 `_algo_n_eval()`（`algorithms/nsga2_core.py`）读取，并在两处取后备。升级 pymoo 大版本时必须重新验证该封装，否则预算审计会静默失效。

**代码版本与数据版本戳的一致性声明**：`DynamicDriver.code_version` 当前取值为 `"v2.4"`，但**磁盘上 v2.4 全部 126 个原始 CSV 与 3240 条运行元数据记录标的是 `"v2.2"`**（由 commit `50adf19c5c5a` 生成）。原因是 `run_comparison.py` 仍显式传入 `code_version="v2.2"`。这属于"版本戳漂移"，不影响数值结果（代码逻辑在 v2.2→v2.4 之间未被改变行为），但本论文如实说明，并在引用数据时同时给出**数据生成提交**与**代码分析提交**。v2.5 两阶段实验的数据统一标为 `"v2.5"`。

## 4.2 代码结构

```text
LDT-NSGA-II/
├── algorithms/                    # 6 个算法 + 驱动器 + 共享数学设施
│   ├── nsga2_core.py              # make_nsga2 工厂 + DynamicDriver + n_eval 兼容
│   ├── ldmt_nsga2.py              # LDT-NSGA-II 主算法 + LDTNSGA2Alpha1（v2.5）
│   ├── distribution_features.py   # normalize / denormalize / extract_features / 非支配集
│   ├── environment_manager.py     # EnvironmentManager + ChangeDetector（D_k）
│   ├── transfer_strategy.py       # 可插拔策略层（默认实验未用）
│   ├── ar_nsga2.py                # AR(1) 预测基线
│   ├── gmm_transfer_nsga2.py      # GMM 分布迁移基线
│   ├── pps_nsga2.py               # PPS 中心预测基线
│   └── dnsga2_response.py         # 静态 + 响应式基线
├── problems/                      # 9 个测试问题 + 参考前沿
│   ├── base.py                    # DynamicProblem：时间规则 + 共享固定探针
│   ├── fda1.py / fda3.py / fda4.py
│   ├── dmop1.py / dmop2.py / dmop3.py      # 项目变体（DMOP3 为静态对照）
│   ├── dmop1a.py / dmop2a.py / dmop3a.py   # 权威 dMOP1/2/3（v2.5 新增）
│   └── true_fronts.py             # 解析参考前沿（按 t 缓存）
├── evaluation/                    # igd / hv / statistics / friedman_test / reference_front
├── experiments/
│   ├── common.py                  # RunConfig、算法/问题注册表、ExperimentRunner（断点续跑+元数据）
│   ├── run_comparison.py          # 主比较实验
│   ├── run_supplementary.py       # DMOP 补充实验
│   ├── run_ablation.py            # 消融实验
│   ├── run_sensitivity.py         # 敏感性实验
│   └── run_v25_alpha.py           # v2.5 两阶段实验驱动（α=1.0 对照 + 权威 dMOP）
├── scripts/
│   ├── check_dmop_authoritative.py  # 权威 dMOP 解析自检（v2.5）
│   ├── thesis_analysis.py           # 论文级汇总与统计（v2.5，独立于 analyze.py）
│   └── make_thesis_figures.py       # 论文插图生成
├── tests/                         # 9 个测试文件、63 个单元测试
├── main.py                        # 一键流水线入口
├── analyze.py                     # v2.4 统计汇总 + 图表 + JSON 导出
├── build_report.py / serve_report.py
└── results/
    ├── raw/                       # v2.4 原始逐代 CSV（126 文件，冻结，不修改）
    ├── raw_v25/                   # v2.5 原始逐代 CSV（独立目录，禁止与 raw 混写）
    ├── summary/                   # v2.4 汇总（冻结）
    ├── summary_v25/               # v2.5 汇总（由 scripts/thesis_analysis.py 生成）
    ├── figures/                   # 汇总图表
    └── figures_thesis/            # 论文插图
```

**目录隔离原则**：v2.5 的新实验被强制写入 `results/raw_v25/`。`analyze.py` 用 glob 读取 `results/raw/comparison_*.csv` 等模式，若把新算法写进 `results/raw` 会静默污染冻结的 v2.4 汇总。为此 `run_v25_alpha.py` 内置了拒绝逻辑：当 `--out` 指向 `results/raw` 或任何 `raw_legacy*` 目录时直接退出并报错。

## 4.3 实验协议与因子设计

**固定协议（四组实验共同）**：种群规模 $N = 100$、总代数 $n_{gen} = 300$、每环境代数 $\tau_t \in \{10, 25, 50\}$、环境数 $n_t = 20$、独立随机种子 $0\text{–}29$（30 个**对齐**种子）、环境段数 $n_{envs} = \max(1, \lfloor n_{gen}/\tau_t\rfloor)$（$\tau_t = 25$ 时为 12 段）。

**预算核算**：每个环境段的函数评价次数严格等于 $N \times \tau_t$（实测 $\tau_t = 10/25/50$ 分别对应 1000/2500/5000 次），单次运行总评价次数 $= N \times n_{gen} = 30\,000$，且**对六个算法完全相同**。这一等值是第 3 章"策略与内核分离"设计的直接后果，也是算法间可比性的基础。

| 实验组 | 脚本 | 因子设计 | 运行次数 | 原始 CSV |
|---|---|---|---|---|
| 主比较 comparison | `run_comparison.py` | 3 问题（FDA1/FDA3/FDA4）× $\tau \in \{10,25,50\}$ × 6 算法 × 30 种子 | 1620 | 54 |
| 补充 supplementary | `run_supplementary.py` | 3 问题（DMOP1/2/3）× $\tau = 25$ × 6 算法 × 30 种子 | 540 | 18 |
| 消融 ablation | `run_ablation.py` | 3 问题 × $\tau = 25$ × 4 变体（full / A 关趋势外推 / B 关分布迁移 / C $\kappa=0$）× 30 种子 | 360 | 12 |
| 敏感性 sensitivity | `run_sensitivity.py` | $\tau = 25$；$N \in \{50,100,200\}$、$n_t \in \{10,20,40\}$、$\rho \in \{0.25,0.5,0.75\}$、$\alpha \in \{0,0.25,0.5,0.75,1.0\}$ | 1260 | 42 |
| **v2.4 合计** | | | **3780** | **126** |

表 4-1：v2.4 四组实验的因子设计与规模（已由独立脚本按磁盘数据复核）。

敏感性的 1260 次由四项研究构成，逐项为：$N$ 三档 × 3 问题 × 30 种子 = 270；$n_t$ 三档 × 3 × 30 = 270；$\rho$ 三档 × 3 × 30 = 270；$\alpha$ 五档 × 3 × 30 = 450；合计 270+270+270+450 = **1260**。

**v2.5 新增两阶段实验**（写入 `results/raw_v25/`）：

| 阶段 | 脚本 | 因子设计 | 运行次数 |
|---|---|---|---|
| 阶段一：$\alpha = 1.0$ 对照 | `run_v25_alpha.py --phase comparison` | 3 问题（FDA1/FDA3/FDA4）× $\tau \in \{10,25,50\}$ × 6 算法（5 基线 + `ldt_alpha1`）× 30 种子 | 1620 |
| 阶段二：权威 dMOP 补充 | `run_v25_alpha.py --phase supplementary` | 3 问题（DMOP1A/DMOP2A/DMOP3A）× $\tau = 25$ × 6 算法 × 30 种子 | 540 |
| 阶段二附加：$\alpha=0.5$ 对照 | `run_v25_alpha.py --phase supplementary --algos ldt_nsga2` | 3 问题 × $\tau=25$ × 1 算法（$\alpha=0.5$）× 30 种子 | 90 |
| **v2.5 合计** | | | **2250** |

表 4-2：v2.5 两阶段实验的规模。**阶段一全部重跑五个基线**（而不是复用 v2.4 的基线结果），以保证与 `ldt_alpha1` 完全同批、同代码、同环境生成；因此 $\alpha = 1.0$ 与 $\alpha = 0.5$ 的对照在"基线数值"上也做了一次独立的重复验证（第 5 章 5.7 节）。**阶段二附加的 90 次**用于补齐 LDT 的 $\alpha=0.5$ 变体在权威 dMOP 问题上的读数——该变体不在 `run_v25_alpha.py` 的默认算法列表中（默认列表用 `ldt_alpha1` 替代了 `ldt_nsga2`），但若不补齐就无法回答"在这一族问题上外推步长是否同样重要"，因此单独补跑；其结果见 5.9.2 节与表 5-19。

至此全部分析数据为 $3780 + 2250 = 6030$ 次独立运行。

## 4.4 断点续跑与运行级元数据

**断点续跑**：`ExperimentRunner._existing_seeds(path)` 读取目标 CSV 中已经出现的 `seed` 集合，命中即跳过对应的 (algorithm, problem, $\tau$, seed) 组合。`run_ablation.py` 与 `run_sensitivity.py` 有等价的内联逻辑。这一机制使"全天候长跑"成为可能，也意味着**任何一步中断后重跑同一命令都是安全的、幂等的**。

**运行级元数据**：每次运行向 `results/raw/run_meta.jsonl` 追加一行 JSON，字段包括算法规范名与 CLI 名、问题、$\tau_t$、$n_t$、$N$、$n_{gen}$、种子、`n_eval_total` 与**逐环境评估数** `n_eval_per_env`、耗时、**git 提交号（前 12 位）**、`code_version`、状态，以及 Python/numpy/pymoo/sklearn 版本与平台信息。失败记录单独写入 `failed_runs.jsonl`。

**失败判定**：运行器对每一行的 IGD 与 HV 做有限性校验，出现 `NaN/Inf` 即抛异常并记为失败。v2.4 全部 3780 次运行中 `results/raw/failed_runs.jsonl` **不存在**，即 0 失败；v2.5 阶段一、阶段二亦为 0 失败（见 5.1.2 节的数据状态核验）。

## 4.5 数据管线与产物

**逐代原始数据**：`results/raw/*.csv` 与 `results/raw_v25/*.csv`，每行一代，字段恰为 18 项（`ROW_FIELDS`）：

```text
algorithm, problem, seed, population_size, total_generations, n_t, tau_t,
environment_index, generation, IGD, HV, generation_runtime_s,
migration_runtime_s, n_evaluations, code_version, lambda_k,
boundary_clip_ratio, n_nd
```

其中 `migration_runtime_s` 只在**变化后的第一代**（`gen_in_env == 0` 且 $k > 0$）写入实测迁移耗时，其余代记 0——因此对该列按运行求和恰好等于整轮的真实迁移总时间，不会因重复累加而放大。

**论文级汇总**：`scripts/thesis_analysis.py` 从上述两个目录读入全部逐代数据，按 5.1.3 节写死的口径计算每次运行的稳态 MIGD / MHV、变化后 MIGD、全代 MHV、恢复代数与迁移耗时，输出到 `results/summary_v25/`：

| 产物 | 内容 |
|---|---|
| `per_run_metrics.csv` | 每次运行一行的全部汇总指标 |
| `summary_v24.csv` / `summary_v25.csv` | 每（问题, $\tau$, 算法）的均值、标准差、中位数、MHV、变化后 MIGD、迁移耗时、运行数 |
| `wilcoxon_holm_v24.csv` / `wilcoxon_holm_v25.csv` | 配对 Wilcoxon + Holm 校正的逐配置检验结果 |
| `friedman_cells_*.csv` / `friedman_overall_*.csv` / `nemenyi_gaps_*.csv` | Friedman 检验、平均秩、Nemenyi 临界差对比 |
| `alpha1_vs_alpha05.csv` | $\alpha = 1.0$ 与 $\alpha = 0.5$ 的配对对照 |
| `recovery_summary.csv` | 恢复代数汇总 |
| `integrity.json` / `integrity_runs_by_prefix.csv` | 数据完整性核验（行数、文件数、运行数、版本戳） |
| `analysis_report.md` | 上述全部表格的 Markdown 汇总 |

**关键设计决定**：`thesis_analysis.py` **不写 `results/summary/`**，也不复用 `analyze.py` 的中间产物，而是从原始逐代 CSV 从零重算。这样论文中的每一个数字都有一条"原始数据 → 口径代码 → 论文表格"的完整链条；同时 `results/summary/` 中由 `analyze.py` 生成的 v2.4 冻结汇总保持原样，可用于交叉验证。第 5.1.3 节报告了两种独立实现的互检结果。

**图表**：`results/figures/` 由 `analyze.py` 生成（逐问题 MIGD 箱线图、MIGD–$\tau$ 曲线、消融分组条形图、敏感性图、迁移耗时图、各问题恢复曲线）；`results/figures_thesis/` 由 `scripts/make_thesis_figures.py` 生成论文专用插图（问题示意、方法分类、算法框架、迁移示意、跨问题 MIGD、收敛曲线、恢复轨迹、迁移耗时、消融、$N$ 敏感性、Friedman–Nemenyi 临界差图）。

## 4.6 测试体系（63 个单元测试）

运行方式：`python -m pytest -q`。测试总数经 `--collect-only` 实测为 **63**，全部通过是发布实验数据的先决条件。

| 测试文件 | 数量 | 覆盖内容 |
|---|---|---|
| `test_imports.py` | 5 | 全模块可导入、6 个 canonical id 注册表、问题实例探针形状与边界 |
| `test_v22_fixes.py` | 10 | FDA4 与 PlatEMO 逐点对照（`atol=1e-12`）、Holm 已知值、静态/响应式基线行为、归一化往返、自适应 $\sigma$ 上限、共享探针 |
| `test_k1_immigrant_fix.py` | 3 | K1 修复回归：AR1/GMM/PPS 移民注入后初始种群确含移民行且定容正确 |
| `test_v21_features.py` | 9 | PPS 冒烟、DMOP3、Friedman–Nemenyi、可插拔策略、多源融合 |
| `test_smoke.py` | 8 | 6 算法端到端冒烟（FDA1；LDT 另测 FDA3/FDA4） |
| `test_environment_change.py` | 4 | 变化判定、FDA1/FDA4 跨环境目标变化、FDA4 单位球面前沿 |
| `test_problem_dimensions.py` | 4 | 6 个问题的维数与目标数 |
| `test_indicators.py` | 5 | IGD 自距离为 0、按 $t$ 取前沿、HV 已知集合与参考点 |
| `test_v24_fixes.py` | 15 | HV 回退对拍 pymoo、Friedman 并列秩与 Nemenyi 双表、MomentTransfer 归一化一致性、DMOP 静/动态行为、四算法真实 `boundary_clip_ratio`、FDA4 前沿 507 点 |
| **合计** | **63** | |

表 4-3：测试体系分工（逐文件计数经实测核对，与代码仓库一致）。

测试体系覆盖四类风险：（a）**公式转写风险**——用 PlatEMO 源码逐点对照锁定 FDA4；（b）**参数失效风险**——用回归测试锁定 K1 类"死参数"缺陷；（c）**口径漂移风险**——用已知值测试锁定 Holm 校正与 Nemenyi 查表；（d）**接口契约风险**——用冒烟测试锁定 6 个算法都能端到端跑通。

v2.5 新增的权威 dMOP 实现采用**独立自检脚本**而非单元测试来验证：`scripts/check_dmop_authoritative.py` 把决策向量放在解析导出的 $PS(t)$ 上，检查其目标值是否落在解析 $PF(t)$ 上（最大偏差 0.000e+00），并检查三者与项目变体确实不同、$PF$ 的时变/时不变性质符合问题类型。

## 4.7 从零复现：逐步命令

以下命令均在项目根目录下执行（Windows PowerShell；`python` 需为装有 4.1 节依赖的解释器）。

**第 1 步：安装依赖并跑测试**

```powershell
pip install -r requirements.txt
python -m pytest -q                     # 预期 63 passed
```

**第 2 步：单次冒烟（分钟级，务必用独立输出目录）**

```powershell
python -m experiments.run_comparison --problems FDA3 --taus 25 --seeds 0 `
    --out results/_smoke --n-gen 50 --limit-seeds 1
```

⚠️ 冒烟数据使用的是缩小预算（50 代），**绝不可与正式数据混放**；`--out results/_smoke` 保证隔离。

**第 3 步：v2.4 四组实验（全量约数小时，可中断续跑）**

```powershell
python -m experiments.run_comparison --problems FDA1,FDA3,FDA4 --taus 10,25,50 `
    --seeds 0-29 --out results/raw
python -m experiments.run_supplementary --problems DMOP1,DMOP2,DMOP3 --seeds 0-29 `
    --out results/raw --algos nsga2_static,dnsga2_inspired,ar1_predictive,gmm_inspired,pps,ldt_nsga2
python -m experiments.run_ablation
python -m experiments.run_sensitivity
```

⚠️ **已知的可复现性陷阱**：`run_supplementary.py` 自身默认只跑 `DMOP1,DMOP2`（=360 次），复现论文的 540 次必须显式加 `--problems DMOP1,DMOP2,DMOP3`；`main.py` 未覆盖该默认值。该问题记入已知问题清单 K14。

**第 4 步：v2.5 两阶段实验**

```powershell
# 阶段一 α=1.0 对照（1620 次）；按 (问题,τ) 拆成 9 个进程并行可缩短到约 0.5 小时
python -m experiments.run_v25_alpha --problems FDA1,FDA3,FDA4 --taus 10,25,50 `
    --seeds 0-29 --out results/raw_v25 --phase comparison
# 阶段二 权威 dMOP 补充（540 次）
python -m experiments.run_v25_alpha --problems DMOP1A,DMOP2A,DMOP3A --taus 25 `
    --seeds 0-29 --out results/raw_v25 --phase supplementary
```

**第 5 步：权威 dMOP 自检 + 论文级统计**

```powershell
python -m scripts.check_dmop_authoritative    # 解析 PS/PF 一致性自检
python -m scripts.thesis_analysis             # 生成 results/summary_v25/ 全部论文表格
python analyze.py                             # 可选：重新生成 v2.4 冻结汇总与图表
python -m scripts.make_thesis_figures         # 生成论文插图
```

**第 6 步：一键流水线（v2.4）**

```powershell
python main.py        # pytest → analyze.py → build_report.py → 自动打开报告
```

`main.py` 支持 `--with-experiments`（先跑四组完整实验）、`--skip-tests`、`--no-open`；任一步失败即中止并打印错误。

## 4.8 复现性边界：本论文明确不能保证的部分

为避免过度承诺，本节列出本实验系统**未能**保证的复现性边界：

1. **浮点级完全复现未验证**：同一代码、同一平台、同一种子应产生相同结果（pymoo 与 numpy 的随机流已固定），但本论文未做跨机器、跨 BLAS 实现的逐位比对。跨平台复现时 MIGD 数值可能存在末位差异。
2. **`run_meta.jsonl` 的 git 提交号依赖仓库状态**：若在未提交的工作树上运行，`git_commit()` 返回的提交号不包含未提交改动，因此"提交号 → 结果"的映射仅在干净工作树上严格成立。
3. **`build_report.py` 生成的 HTML 报告中含旧版叙述文字**：`build_report.py` 第 972–984 行等处仍保留 v1 时期的硬编码说明（4 算法、536,700 行、平均 MIGD 0.0378 等），与 v2.4 的 6 算法 / 1,134,000 行不符。报告中的**表格与图表数据取自 JSON（正确）**，但**叙述文字不可引用**。本论文的全部数字均取自逐代原始 CSV 的独立重算，不依赖该报告的叙述。该问题记入已知问题清单 K15。
4. **`configs/*.yaml` 为过时镜像**：四个 YAML 配置（`main_experiment.yaml` 仅 4 算法、`code_version: v1.0`；`sensitivity.yaml` 的 $N \in [50,100,150]$ 且仅 2 项研究；`ablation.yaml` 写 70/20/10 与 `no_feature_simplification`）**未被任何脚本读取**，与 v2.4 的实际实验不一致，不可作为实验配置依据。记入已知问题清单 K16。
5. **`README.md` 的 $\sigma_{\min}$ 描述过时**：README 写 `sigma_min = 0.001*(u-l)`（真实尺度），而代码是**归一化单位**的 $10^{-3}$。本论文统一采用代码口径。记入已知问题清单 K17。

## 4.9 本章小结

本章给出了实验系统的完整设计：固定协议（$N=100$、$n_{gen}=300$、$n_t=20$、30 个对齐种子）与严格的预算核算（每环境 $N\tau_t$ 次评价，算法间完全相同）；四组共 3780 次 v2.4 运行与两阶段共 2250 次 v2.5 运行，全部支持断点续跑并落盘运行级元数据；一条从原始逐代 CSV 到论文表格的独立重算管线（`scripts/thesis_analysis.py`），以及 63 个单元测试与权威 dMOP 解析自检。同时**明确划出复现性边界**——浮点级跨平台一致性未验证、版本戳漂移、报告叙述文字的旧版残留、过时的 YAML 与 README 描述——这些都以已知问题清单的形式如实披露，而非隐去。第 5 章报告实验结果。


\newpage

# 第 5 章 实验结果与分析

## 5.1 实验设置、数据状态与口径

### 5.1.1 规模与数据状态

v2.4 四组实验共 **3780 次运行记录**，**失败 0 次**（`results/raw/failed_runs.jsonl` 不存在）。v2.5 两阶段实验再加 **2250 次**（阶段一 $\alpha=1.0$ 对照 1620 次、阶段二权威 dMOP 补充 540 次主表 + 90 次 $\alpha=0.5$ 对照 = 630 次），亦为 0 失败。全部分析数据为 **6030 次运行记录**。

| 实验组 | 配置 | 运行记录数 | 原始 CSV |
|---|---|---|---|
| 主比较 comparison（v2.4） | 3 问题（FDA1/FDA3/FDA4）× $\tau \in \{10,25,50\}$ × 6 算法 × 30 种子 | 1620 | 54 |
| 补充 supplementary（v2.4） | 3 问题（DMOP1/2/3）× $\tau=25$ × 6 算法 × 30 种子 | 540 | 18 |
| 消融 ablation（v2.4） | 3 问题 × $\tau=25$ × 4 变体 × 30 种子 | 360 | 12 |
| 敏感性 sensitivity（v2.4） | $\tau=25$；$N$、$n_t$、$\rho$、$\alpha$ 四组扫描 | 1260 | 42 |
| 阶段一 $\alpha=1.0$ 对照（v2.5） | 3 问题 × $\tau \in \{10,25,50\}$ × 6 算法 × 30 种子 | 1620 | 54 |
| 阶段二 权威 dMOP 补充（v2.5） | 3 问题（DMOP1A/2A/3A）× $\tau=25$ × 6 算法 × 30 种子 | 540 | 18 |
| 阶段二 附加 $\alpha=0.5$ 对照（v2.5） | 3 问题 × $\tau=25$ × 1 算法（`ldt_nsga2`）× 30 种子 | 90 | 3 |
| **合计** | | **6030** | **201** |

**数据完整性的独立核验。** 论文写作期间，数据的完整性由两个相互独立的流程核验：

1. `scripts/thesis_analysis.py` 直接从原始逐代 CSV 重新计算（输出 `results/summary_v25/integrity.json`）：v2.4 主比较 54 个文件、补充 18 个、消融 12 个、敏感性 42 个，共 126 个文件；每组每个（算法, 问题, $\tau$）的种子数均为 **30，无缺种**；逐代行数为 300（或敏感性中 $n_t$ 变化时的 600/3600），**无截断**；IGD/HV 全部有限。
2. 一个独立的数据复核任务（不读取 `results/summary/` 的任何汇总值，仅从逐代 CSV 从零重算）逐项比对了论文表格中的 269 个数值断言，结果 **265 项一致、4 项不一致、0 项无法判定**。4 项不一致已按复核结果修正并记录在 5.10 节的已知问题清单（K18–K21）。

### 5.1.2 三个必须明确的度量口径

本文所有表格与统计检验使用的"**MIGD**"是**稳态口径**，其精确定义（`analyze.py::per_run_metrics`、`scripts/thesis_analysis.py` 同口径）为：

> 对**一次运行**，取**每个环境段最后一代**（即下一次环境变化前的末代种群）的 IGD，对**全部环境段取算术平均**，记为这次运行的 MIGD；再对 30 个种子取平均。$\tau=25$ 时每次运行覆盖 12 个环境段。

类似的辅助口径定义如下：

| 口径 | 定义 | 用途 |
|---|---|---|
| **稳态 MIGD**（主指标） | 各环境段**末代** IGD 的段间平均，再对种子平均 | 正文表格、统计检验、Friedman 排名 |
| **变化后 MIGD** | 各**已变化**环境段**首代** IGD 的段间平均，再对种子平均 | 衡量响应策略的初始命中质量 |
| **稳态 MHV** | 与稳态 MIGD 同口径，作用于 HV | 辅助验证（越大越好） |
| **全代 MHV** | 一次运行**全部代** HV 的算术平均 | 辅助诊断 |

**必须声明一处与既有材料的差异**：本项目早期材料（以及项目 `README`）曾称"系统同时报告变化后 MIGD 与**全代 MIGD**"。经复核，`analyze.py` 中**不存在**"全代 IGD 均值"这一口径——第三个聚合量是**全代 MHV**；"全代 IGD 均值"只出现在开发期诊断脚本 `_check_pilot.py` 中（读取 `results/raw_pilot`，不产出任何汇总）。本论文据此更正表述，并以 `scripts/thesis_analysis.py` 同时输出全代 MIGD 供交叉参考（该列只用于诊断，不用于结论）。

**另一处必须声明的口径细节**：表中"迁移耗时"是**一次运行（300 代、12 个环境段）的迁移总耗时**，即对 `migration_runtime_s` 列按运行求和后的均值，而非单代均值。该列只在变化后的第一代写入实测值，故按列求和恰为真实总量。

### 5.1.3 两种独立实现的互检

论文写作期间，v2.4 主比较表的 54 个（问题 × $\tau$ × 算法）MIGD 单元格被两套相互独立的实现分别计算：(a) 项目原有的 `analyze.py`（输出 `results/summary/summary_comparison.csv`）；(b) 为论文新写的 `scripts/thesis_analysis.py`（输出 `results/summary_v25/summary_v24.csv`，从原始 CSV 从零重算、不读前者）。两套实现的全部 54 个单元格在四位小数上完全一致（例：FDA1 $\tau=10$ PPS 0.0122、FDA1 $\tau=25$ LDT 0.0074、FDA4 $\tau=10$ LDT 0.1236），由此确认"论文数字 → 原始数据"的链条无中间环节失真。

## 5.2 FDA1 上的结果（二目标，PF 时不变，PS 正弦平移）

表 5-1：FDA1 各算法稳态 MIGD（均值 ± 标准差，$n=30$；加粗为该列最优）。

| 算法 | $\tau=10$ | $\tau=25$ | $\tau=50$ |
|---|---|---|---|
| nsga2_static | 0.0204 ± 0.0012 | 0.0086 ± 0.0005 | 0.0053 ± 0.0001 |
| dnsga2_inspired | 0.0211 ± 0.0011 | 0.0087 ± 0.0005 | 0.0054 ± 0.0001 |
| ar1_predictive | 0.0173 ± 0.0032 | 0.0065 ± 0.0005 | 0.0050 ± 0.0001 |
| gmm_inspired | 0.0203 ± 0.0010 | 0.0084 ± 0.0006 | 0.0053 ± 0.0001 |
| **pps** | **0.0122 ± 0.0014** | **0.0062 ± 0.0005** | **0.0049 ± 0.0001** |
| ldt_nsga2（$\alpha=0.5$） | 0.0156 ± 0.0010 | 0.0074 ± 0.0005 | 0.0052 ± 0.0001 |

表 5-2：FDA1 上 LDT（$\alpha=0.5$）对各基线的配对 Wilcoxon + Holm 校正（$n=30$，双侧；方向按均值判读；$p_{adj}<0.001$ 记 $***$、$<0.01$ 记 $**$、$<0.05$ 记 $*$）。

| 对比 | $\tau=10$ | $\tau=25$ | $\tau=50$ |
|---|---|---|---|
| LDT vs nsga2_static | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$2.24\times10^{-4}$）$***$ |
| LDT vs dnsga2_inspired | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$2.53\times10^{-6}$）$***$ |
| LDT vs gmm_inspired | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$1.10\times10^{-2}$）$**$ |
| LDT vs ar1_predictive | LDT 更好（$1.24\times10^{-3}$）$**$ | AR1 更好（$9.31\times10^{-9}$）$***$ | AR1 更好（$3.54\times10^{-8}$）$***$ |
| LDT vs pps | PPS 更好（$9.31\times10^{-9}$）$***$ | PPS 更好（$9.31\times10^{-9}$）$***$ | PPS 更好（$2.98\times10^{-8}$）$***$ |

**关于 $9.31\times10^{-9}$ 的说明（方法学细节）**：该值反复出现并非巧合。$n=30$ 的配对 Wilcoxon 符号秩检验在"全部 30 个差值同号"时达到最小可能的原始 $p = 2^{-29} \approx 1.863\times10^{-9}$；经 Holm 校正（$m=5$ 次比较、该值排名最小时乘 $m$）得 $5 \times 2^{-29} = 9.313\times10^{-9}$。因此 $9.31\times10^{-9}$ 是**"完全分离 + Holm 校正"的饱和值**，应理解为"差异极大"，而不是一个精确测量的概率。本论文如实给出该数值并说明其含义。

## 5.3 FDA3 上的结果（二目标，PF 随时间变化）

表 5-3：FDA3 各算法稳态 MIGD（均值 ± 标准差，$n=30$）。

| 算法 | $\tau=10$ | $\tau=25$ | $\tau=50$ |
|---|---|---|---|
| nsga2_static | 0.0215 ± 0.0013 | 0.0092 ± 0.0005 | 0.0055 ± 0.0001 |
| dnsga2_inspired | 0.0227 ± 0.0011 | 0.0092 ± 0.0005 | 0.0055 ± 0.0001 |
| ar1_predictive | 0.0184 ± 0.0026 | 0.0070 ± 0.0005 | 0.0052 ± 0.0001 |
| gmm_inspired | 0.0220 ± 0.0011 | 0.0090 ± 0.0006 | 0.0055 ± 0.0001 |
| **pps** | **0.0132 ± 0.0013** | **0.0066 ± 0.0005** | **0.0051 ± 0.0001** |
| ldt_nsga2（$\alpha=0.5$） | 0.0177 ± 0.0009 | 0.0079 ± 0.0005 | 0.0054 ± 0.0001 |

FDA3 的检验结论与 FDA1 同向：LDT（$\alpha=0.5$）在三档 $\tau$ 上均显著优于 static、dnsga2 与 gmm（$\tau=10/25$ 为 $p_{adj}=9.31\times10^{-9}$；$\tau=50$ 对 gmm 为 $7.94\times10^{-5}$）；对 PPS 三档均显著劣后（$9.31\times10^{-9}$）；对 AR1，$\tau=10$ **不显著**（$p_{adj}=0.198$，而均值 0.0177 对 0.0184 实为 LDT 更低），$\tau=25/50$ 显著劣后（$9.31\times10^{-9}$ 与 $1.47\times10^{-8}$）。

## 5.4 FDA4 上的结果（三目标，PF 时不变但非退化，PS 旋转型移动）

表 5-4：FDA4 各算法稳态 MIGD（均值 ± 标准差，$n=30$）。

| 算法 | $\tau=10$ | $\tau=25$ | $\tau=50$ |
|---|---|---|---|
| nsga2_static | 0.2101 ± 0.0161 | 0.1020 ± 0.0035 | 0.0764 ± 0.0020 |
| dnsga2_inspired | 0.1628 ± 0.0096 | 0.1027 ± 0.0044 | 0.0759 ± 0.0022 |
| ar1_predictive | 0.1507 ± 0.0112 | 0.0987 ± 0.0082 | 0.0746 ± 0.0018 |
| gmm_inspired | 0.1595 ± 0.0090 | 0.1011 ± 0.0044 | 0.0757 ± 0.0018 |
| **pps** | **0.1167 ± 0.0075** | **0.0881 ± 0.0046** | **0.0736 ± 0.0019** |
| ldt_nsga2（$\alpha=0.5$） | 0.1236 ± 0.0047 | 0.0904 ± 0.0039 | 0.0745 ± 0.0017 |

表 5-5：FDA4 上 LDT（$\alpha=0.5$）对各基线的 Wilcoxon-Holm 检验。

| 对比 | $\tau=10$ | $\tau=25$ | $\tau=50$ |
|---|---|---|---|
| LDT vs nsga2_static | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$2.24\times10^{-4}$）$***$ |
| LDT vs dnsga2_inspired | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$1.23\times10^{-3}$）$**$ |
| LDT vs gmm_inspired | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$1.82\times10^{-3}$）$**$ |
| LDT vs ar1_predictive | LDT 更好（$9.31\times10^{-9}$）$***$ | LDT 更好（$9.31\times10^{-9}$）$***$ | 不显著（$p_{adj}=0.529$；均值 0.0745 对 0.0746） |
| LDT vs pps | PPS 更好（$3.11\times10^{-7}$）$***$ | PPS 更好（$3.17\times10^{-6}$）$***$ | PPS 更好（$p_{adj}=0.013$）$*$ |

**关键观察（v2.4，$\alpha=0.5$）**：K1 修复为三个预测基线补回随机移民后，PPS 的中心外推获得了多样性补充，在 FDA4 全部三档均以均值第一。但有两点值得强调：

1. **LDT 的种子间方差更小**：$\tau=10$ 时 LDT 标准差 0.0047、PPS 0.0075（方差比约 0.39）；$\tau=25$ 时 0.0039 对 0.0046；$\tau=50$ 时 0.0017 对 0.0019。即 LDT 的运行稳定性在三个档位上一致更好。
2. **在 Friedman–Nemenyi 整体排序下，LDT 与 PPS 在 FDA4 全部三档不可区分**（平均秩差 0.80 / 0.73 / 1.00，均小于 $CD=1.37668$，见 5.5.1 节）。

MHV 方向一致：FDA4 $\tau=10$ 时 PPS 0.5976、LDT 0.5707，两者远高于第三名 dnsga2 的 0.4855（参考点 $(1.1,1.1,1.1)$，越大越好）。

## 5.5 跨问题综合分析（v2.4）

![图 5-1 不同变化频率下各算法的稳态 MIGD](../../LDT-NSGA-II/results/figures_thesis/fig5_1_migd_vs_tau.png)

图 5-1：不同变化频率 $\tau$ 下六算法在 FDA1/FDA3/FDA4 上的稳态 MIGD（均值 ± 标准差，$n=30$；误差帽为标准差）。K1 修复后 PPS 在全部九个"问题 × $\tau$"组合均值第一，LDT 稳居第二；两者差距随 $\tau$ 增大收敛。

### 5.5.1 Friedman 平均秩与相对差距

9 个（问题, $\tau$）组合的 Friedman 检验**全部显著**（$\chi^2_F$ 介于 64.2–139.4，$p<10^{-9}$；Nemenyi 临界差 $CD = 1.37668$）。平均秩（越小越好，$k=6$、$N=30$）：

| 配置 | PPS | AR1 | **LDT** | GMM | static | dnsga2 |
|---|---|---|---|---|---|---|
| FDA1 $\tau=10$ | **1.00** | 3.00 | 2.30 | 4.47 | 4.63 | 5.60 |
| FDA1 $\tau=25$ | **1.00** | 2.00 | 3.00 | 4.33 | 5.30 | 5.37 |
| FDA1 $\tau=50$ | **1.27** | 1.80 | 3.60 | 4.17 | 4.90 | 5.27 |
| FDA3 $\tau=10$ | **1.00** | 2.77 | 2.43 | 4.83 | 4.37 | 5.60 |
| FDA3 $\tau=25$ | **1.00** | 2.00 | 3.00 | 4.37 | 5.03 | 5.60 |
| FDA3 $\tau=50$ | **1.27** | 1.77 | 3.50 | 4.73 | 4.87 | 4.87 |
| FDA4 $\tau=10$ | **1.10** | 3.30 | 1.90 | 4.20 | 6.00 | 4.50 |
| FDA4 $\tau=25$ | **1.13** | 4.00 | 1.87 | 4.30 | 4.63 | 5.07 |
| FDA4 $\tau=50$ | **1.67** | 3.10 | 2.67 | 4.27 | 4.80 | 4.50 |

PPS 在全部 9 个组合上平均秩第一；LDT 在 FDA1/FDA3 的 $\tau=10$ 两档列第二（2.30 / 2.43，优于同档 AR1 的 3.00 / 2.77），在 FDA4 全部三档列第二（1.90 / 1.87 / 2.67）。

LDT 相对 PPS 的稳态 MIGD 百分差距（正 = 劣于 PPS）：

| 问题 | $\tau=10$ | $\tau=25$ | $\tau=50$ |
|---|---|---|---|
| FDA1 | +27.4% | +19.6% | +5.3% |
| FDA3 | +33.9% | +20.7% | +5.3% |
| FDA4 | +6.0% | +2.6% | +1.4% |

由此得到两个规律：

1. **$\tau$ 越大（变化越慢、每段预算越充足），算法间差距越小**——$\tau=50$ 时 FDA1/FDA3 上 LDT 与 PPS 的差距收敛到 5.3%，FDA4 上仅 1.4%，且六算法全部落在 0.0736–0.0764 的窄带内（极差约 3.8%）。此时"是否使用预测机制"已不重要。
2. **按 Nemenyi 后验判读，LDT 与 PPS 在 4 个组合上不可区分**：FDA1 $\tau=10$（秩差 1.30）与 FDA4 全部三档（0.80 / 0.73 / 1.00），秩差均小于 $CD=1.37668$。这为"LDT 的整体跟踪能力与最优预测基线相当"提供了多问题口径的支持——尽管逐问题的 Wilcoxon 检验在 9/9 组合上判定 PPS 显著更优。

**两种口径宽严差异的说明**：逐问题 Wilcoxon 回答"在某一个（问题, $\tau$）配置下两算法是否可区分"，其统计功效高、更严格；Friedman–Nemenyi 回答"在 9 个配置的整体排序意义上是否可区分"，其比较的是平均秩。前者的严格性来自配对设计（同一随机种子逐对比较），后者的宽缓来自把 9 个配置的信息聚合为秩。两者结论不同并不意味着任何一方有错——它们回答的是不同问题。本论文同时报告两者，并在 5.9 节明确说明适用边界。

![图 8-1 Friedman-Nemenyi 临界差异图](../../LDT-NSGA-II/results/figures_thesis/fig8_1_nemenyi.png)

图 8-1：Friedman + Nemenyi 检验（对 9 个"问题 × $\tau$"组合取平均秩，$k=6$、$N=30$，$CD=1.37668$）。横轴为平均秩（越小越好）；被横线相连的算法对平均秩差小于 $CD$，整体排序意义上不可区分。

### 5.5.2 响应速度与变化后初始质量

以"变化后 IGD 回落至该段末代水平所需代数"为恢复指标（精确定义见 4.6 节与 5.1.2 节；若整段未达标则记为 $\tau_t$），$\tau=25$ 的平均恢复代数：

| 算法 | FDA1 | FDA3 | FDA4 |
|---|---|---|---|
| pps | 10.29 | 13.97 | 9.65 |
| ar1_predictive | 16.02 | 17.03 | 16.05 |
| ldt_nsga2（$\alpha=0.5$） | 21.29 | 21.64 | 18.50 |
| nsga2_static | 约 23.9 | 约 23.7 | 20.64 |
| dnsga2_inspired | 约 23.6 | 约 23.4 | 20.68 |
| gmm_inspired | 约 23.2 | 约 23.1 | 20.71 |

变化后 MIGD（每段首代 IGD，$\tau=25$）同样显示 PPS 初始命中最准：FDA1 为 0.0104（LDT 0.0238）、FDA3 为 0.0129（LDT 0.0310）、FDA4 为 0.1074（LDT 0.1304）。FDA4 $\tau=10$ 时 PPS 0.1359 对 LDT 0.1606。

值得注意的是，K1 修复后 PPS **同样注入 20% 随机移民**却仍恢复最快，说明冷启动差距的主因并非"随机移民稀释"。PPS 采用种群中心一阶差分**全步外推 + 整体平移**（$\hat\mu = 2\mu_{k-1} - \mu_{k-2}$），在 FDA 这类中心近似匀速平移的问题上恰好命中新中心；而 LDT 的趋势中心预测取 $\alpha = 0.5$ 的**半步外推**，在中心方向持续不变的运动中系统性欠外推；此外高斯重采样会以 $\sigma_e$ 重新扩散已收敛的个体，随机比例（30%）也高于 PPS（20%）。**这提示把 $\alpha$ 提向 1.0 是缩小冷启动差距的直接方向**——这正是 5.7 节 v2.5 对照实验的动机。

![图 5-3 环境变化后的 IGD 收敛曲线](../../LDT-NSGA-II/results/figures_thesis/fig5_3_convergence.png)

图 5-3：$\tau=25$ 时各算法在环境变化后的平均 IGD 收敛曲线（横轴为变化后代数，数据为全部 30 种子与全部变化段的平均）。PPS/AR1 起点（冷启动命中）最低；LDT 起点居中、恢复斜率与静态类基线相近。

![图 5-4 变化后的归一化恢复轨迹](../../LDT-NSGA-II/results/figures_thesis/fig5_4_recovery.png)

图 5-4：与图 5-3 同数据按各算法自身首代 IGD 归一化后的恢复轨迹（纵轴 $\mathrm{IGD}(t)/\mathrm{IGD}(0)$），用于比较恢复的**相对速度**而非初始命中精度。

### 5.5.3 计算开销

一次运行（12 个环境段、300 代）的迁移总耗时（30 种子均值，$\tau=25$）：static 约 $2\times10^{-4}$ s、dnsga2 约 $7\times10^{-4}$ s（几乎无计算）；LDT / AR1 / PPS 同处 $10^{-2}$ s 量级（LDT 约 0.015–0.017 s，PPS 约 0.013 s）；**GMM 最贵，约 0.13–0.19 s**（$\tau=10$ 档约 0.47–0.50 s）。

按相同"问题 × $\tau$"配置配对计算，**GMM 的迁移耗时约为 LDT 的 4.4–11 倍**（$\tau=10$ 约 10–11 倍、$\tau=25$ 约 8–9 倍、$\tau=50$ 约 4.4–7.1 倍），而它的 MIGD 并不占优——全协方差混合模型在本基准上性价比最低。LDT 相对 GMM 的轻量优势在 K1 修复后依然成立。

![图 5-5 迁移模块单次耗时对比](../../LDT-NSGA-II/results/figures_thesis/fig5_5_time.png)

图 5-5：九个"问题 × $\tau$"配置下各算法迁移模块的平均耗时（对数刻度）。GMM 全分布迁移为 LDT 的 4.4–11 倍；LDT/AR1/PPS 同处 $10^{-2}$ s 量级；static/dnsga2 无迁移计算、近乎为零。

## 5.6 消融实验（$\tau=25$，定位各组件贡献）

四个变体（定义见 `experiments/run_ablation.py`）：**full** = 完整 LDT；**A** = 关闭趋势外推（$\alpha=0$）；**B** = 关闭分布迁移（迁移个体改为均匀随机）；**C** = 关闭旧种群保留（$\kappa=0$，随机移民比例升为 0.5）。

表 5-6：消融稳态 MIGD（均值 ± 标准差，$n=30$；加粗为该问题最优）。

| 变体 | FDA1 | FDA3 | FDA4 | 相对 full 的退化（FDA1 / FDA3 / FDA4） |
|---|---|---|---|---|
| **full** | **0.00741 ± 0.00053** | 0.00794 ± 0.00053 | 0.09039 ± 0.00389 | — |
| A：关趋势外推 | 0.00884 ± 0.00055 | 0.00943 ± 0.00059 | 0.10106 ± 0.00377 | +19% / +19% / +12% |
| B：关分布迁移 | 0.00965 ± 0.00060 | 0.01027 ± 0.00057 | 0.10223 ± 0.00446 | **+30% / +29% / +13%** |
| C：关旧种群保留（$\kappa=0$） | 0.00746 ± 0.00050 | **0.00788 ± 0.00058** | **0.09018 ± 0.00417** | +1% / −1% / −0.2% |

配套的配对 Wilcoxon 检验（full 对各变体，$n=30$）：

| 对比 | FDA1 | FDA3 | FDA4 |
|---|---|---|---|
| full vs A | $p=1.86\times10^{-9}$ | $p=1.86\times10^{-9}$ | $p=1.86\times10^{-9}$ |
| full vs B | $p=1.86\times10^{-9}$ | $p=1.86\times10^{-9}$ | $p=1.86\times10^{-9}$ |
| full vs C | $p=\mathbf{0.0175}$（full 更好 +0.7%） | $p=\mathbf{0.0384}$（C 更好 −0.8%） | $p=0.584$（不可区分） |

结论（按证据强度排序）：

1. **分布迁移是贡献最大的组件**：B 在三个问题上均最差，且变化后 MIGD 退化更明显（FDA4 从 0.1304 升至 0.2319，**+78%**；FDA1 从 0.0238 升至 0.0568，**+139%**）。这证明"用上一环境非支配解的一阶/二阶矩构造初始种群"这一核心设计确实有效，且其作用**主要在环境变化后的最初几代**。
2. **趋势外推次之但确实有效**：A 在三个问题上一致退化 12–19%，配对检验 $p=1.86\times10^{-9}$。
3. **旧种群保留（$\kappa=0.2$）的证据不支持把它作为默认机制**：C 与 full 的差异在三个问题上**方向不一致**（FDA1 上 full 更好、FDA3/FDA4 上 C 更好），且效应量极小（±1% 以内）。虽然 FDA1 与 FDA3 上配对检验的**原始** $p$ 值小于 0.05（0.0175 与 0.0384），但经 Holm 校正（$m=3$）后分别为 0.0525 与 0.1152，不再显著；FDA4 上 $p=0.584$ 完全不可区分。**诚实推论**：对 FDA 这类光滑整体平移问题，保留旧个体不比把预算让给随机移民更有价值；该组件的价值可能要在非平稳/多峰场景才显现，本文证据不支持把它作为默认机制保留。

![图 6-1 消融实验对比](../../LDT-NSGA-II/results/figures_thesis/fig6_1_ablation.png)

图 6-1：消融实验（$\tau=25$，$n=30$，均值 ± 标准差）。full 为完整 LDT；A 关趋势外推（$\alpha=0$）；B 关分布迁移（迁移个体改为均匀随机）；C 关旧种群保留（$\kappa=0$）。分布迁移（B 的退化）贡献最大，趋势外推（A）次之，旧种群保留（C）无可检测的稳定增益。

## 5.7 v2.5 关键对照实验：$\alpha = 1.0$（全步外推）

### 5.7.1 实验动机

5.5.2 节的分析指出，LDT（$\alpha=0.5$）与 PPS 的差距可能来自**外推步长**而非迁移范式：PPS 本质上是中心的全步一阶差分外推（$\hat\mu = 2\mu_{k-1} - \mu_{k-2}$），而 LDT 默认只外推半步。5.8 节的敏感性实验也从另一侧支持这一猜测——$\alpha$ 从 0 单调增大到 1.0 时 MIGD 在三个问题上一致下降。

若这一猜测成立，则把 $\alpha$ 改为 1.0 之后 LDT 应当追平甚至超过 PPS；若 LDT 仍显著劣后，则说明差距确实来自范式（点预测 vs 矩迁移）。**这是区分两种解释的关键实验**，因此 v2.5 阶段一用同一协议、同一批基线重跑了 $\alpha=1.0$ 变体（`ldt_alpha1`）。

### 5.7.2 结果：$\alpha=1.0$ 显著缩小差距，并在两个档位反超

表 5-7：v2.5 主比较稳态 MIGD（`ldt_alpha1` 为 $\alpha=1.0$；基线为 v2.5 同批重跑；$n=30$；加粗为该列最优）。

| 算法 | FDA1 $\tau=10$ | FDA1 $\tau=25$ | FDA1 $\tau=50$ |
|---|---|---|---|
| nsga2_static | 0.0204 ± 0.0012 | 0.0086 ± 0.0005 | 0.0053 ± 0.0001 |
| dnsga2_inspired | 0.0211 ± 0.0011 | 0.0087 ± 0.0005 | 0.0054 ± 0.0001 |
| ar1_predictive | 0.0173 ± 0.0032 | 0.0065 ± 0.0005 | 0.0050 ± 0.0001 |
| gmm_inspired | 0.0203 ± 0.0010 | 0.0084 ± 0.0006 | 0.0053 ± 0.0001 |
| **pps** | 0.0122 ± 0.0014 | **0.0062 ± 0.0005** | **0.0049 ± 0.0001** |
| **ldt_alpha1（$\alpha=1.0$）** | **0.0120 ± 0.0009** | **0.0062 ± 0.0005** | **0.0049 ± 0.0001** |

| 算法 | FDA3 $\tau=10$ | FDA3 $\tau=25$ | FDA3 $\tau=50$ |
|---|---|---|---|
| nsga2_static | 0.0215 ± 0.0013 | 0.0092 ± 0.0005 | 0.0055 ± 0.0001 |
| dnsga2_inspired | 0.0227 ± 0.0011 | 0.0092 ± 0.0005 | 0.0055 ± 0.0001 |
| ar1_predictive | 0.0184 ± 0.0026 | 0.0070 ± 0.0005 | 0.0052 ± 0.0001 |
| gmm_inspired | 0.0220 ± 0.0011 | 0.0090 ± 0.0006 | 0.0055 ± 0.0001 |
| **pps** | **0.0132 ± 0.0013** | **0.0066 ± 0.0005** | **0.0051 ± 0.0001** |
| ldt_alpha1（$\alpha=1.0$） | 0.0136 ± 0.0009 | 0.0066 ± 0.0005 | 0.0052 ± 0.0001 |

| 算法 | FDA4 $\tau=10$ | FDA4 $\tau=25$ | FDA4 $\tau=50$ |
|---|---|---|---|
| nsga2_static | 0.2101 ± 0.0161 | 0.1020 ± 0.0035 | 0.0764 ± 0.0020 |
| dnsga2_inspired | 0.1628 ± 0.0096 | 0.1027 ± 0.0044 | 0.0759 ± 0.0022 |
| ar1_predictive | 0.1507 ± 0.0112 | 0.0987 ± 0.0082 | 0.0746 ± 0.0018 |
| gmm_inspired | 0.1595 ± 0.0090 | 0.1011 ± 0.0044 | 0.0757 ± 0.0018 |
| pps | 0.1167 ± 0.0075 | 0.0881 ± 0.0046 | **0.0736 ± 0.0019** |
| **ldt_alpha1（$\alpha=1.0$）** | **0.1068 ± 0.0050** | **0.0858 ± 0.0043** | 0.0737 ± 0.0022 |

表 5-8：v2.5 上 `ldt_alpha1` 对 PPS 的配对 Wilcoxon（$n=30$，双侧，未做多重校正的单对比检验）。

| 配置 | $\alpha=1.0$ | PPS | 差距 | $p$ | 判读 |
|---|---|---|---|---|---|
| FDA1 $\tau=10$ | **0.0120** | 0.0122 | −1.6% | 0.0699 | 不显著（$\alpha=1.0$ 均值更低） |
| FDA1 $\tau=25$ | 0.0062 | 0.0062 | +0.3% | 0.516 | 不显著（实质持平） |
| FDA1 $\tau=50$ | 0.0049 | 0.0049 | +0.2% | 0.612 | 不显著（实质持平） |
| FDA3 $\tau=10$ | 0.0136 | 0.0132 | +3.4% | $2.56\times10^{-4}$ | **PPS 显著更优** |
| FDA3 $\tau=25$ | 0.0066 | 0.0066 | +0.9% | $2.99\times10^{-3}$ | **PPS 显著更优** |
| FDA3 $\tau=50$ | 0.0052 | 0.0051 | +0.8% | 0.0327 | 边缘（Holm 后 0.0654 不显著） |
| FDA4 $\tau=10$ | **0.1068** | 0.1167 | −8.5% | $1.86\times10^{-9}$ | **$\alpha=1.0$ 显著更优** |
| FDA4 $\tau=25$ | **0.0858** | 0.0881 | −2.6% | $5.72\times10^{-7}$ | **$\alpha=1.0$ 显著更优** |
| FDA4 $\tau=50$ | 0.0737 | 0.0736 | +0.1% | 0.808 | 不显著（实质持平） |

**核心结论**：

1. 把 $\alpha$ 从 0.5 改为 1.0 后，LDT 与 PPS 的差距从 **+1.4% ~ +33.9%** 缩小到 **−8.5% ~ +3.4%**；
2. 在 **4 个配置上二者统计不可区分**（FDA1 全部三档、FDA4 $\tau=50$）；
3. 在 **FDA4 $\tau=10$ 与 $\tau=25$ 上 LDT（$\alpha=1.0$）显著优于 PPS**（$p<10^{-6}$），FDA1 $\tau=10$ 上均值也更低（$p=0.070$）；
4. 仅在 **FDA3 的 $\tau=10/\tau=25$ 两档上 PPS 仍显著更优**，但差距已从 33.9%/20.7% 收窄到 3.4%/0.9%。

### 5.7.3 $\alpha$ 的配对效应量与 Fiedman 平均秩

在同一种子下逐对比较 $\alpha=1.0$（v2.5）与 $\alpha=0.5$（v2.4）的 30 次运行：

表 5-9：$\alpha=1.0$ 相对 $\alpha=0.5$ 的配对改进（$n=30$，配对 Wilcoxon 双侧）。

| 配置 | $\alpha=1.0$ | $\alpha=0.5$ | 改进 | $p$ |
|---|---|---|---|---|
| FDA1 $\tau=10$ | 0.0120 | 0.0156 | **+23.2%** | $1.86\times10^{-9}$ |
| FDA1 $\tau=25$ | 0.0062 | 0.0074 | **+16.1%** | $1.86\times10^{-9}$ |
| FDA1 $\tau=50$ | 0.0049 | 0.0052 | +4.8% | $1.86\times10^{-9}$ |
| FDA3 $\tau=10$ | 0.0136 | 0.0177 | **+22.8%** | $1.86\times10^{-9}$ |
| FDA3 $\tau=25$ | 0.0066 | 0.0079 | **+16.4%** | $1.86\times10^{-9}$ |
| FDA3 $\tau=50$ | 0.0052 | 0.0054 | +4.3% | $3.73\times10^{-9}$ |
| FDA4 $\tau=10$ | 0.1068 | 0.1236 | **+13.6%** | $1.86\times10^{-9}$ |
| FDA4 $\tau=25$ | 0.0858 | 0.0904 | +5.0% | $1.86\times10^{-9}$ |
| FDA4 $\tau=50$ | 0.0737 | 0.0745 | +1.1% | $1.85\times10^{-2}$ |

**9/9 个配置全部显著改善**，改善幅度 1.1%–23.2%，且变化越快（$\tau$ 越小）改善越大。这与敏感性实验（5.8.2 节）中 $\alpha$ 扫描的结论完全一致。

v2.5 的 Friedman 平均秩（$k=6$、$N=30$、$CD=1.37668$）：

| 配置 | $\alpha=1.0$ | PPS | AR1 | GMM | static | dnsga2 |
|---|---|---|---|---|---|---|
| FDA1 $\tau=10$ | **1.37** | 1.63 | 3.30 | 4.47 | 4.63 | 5.60 |
| FDA1 $\tau=25$ | 1.53 | **1.47** | 3.00 | 4.33 | 5.30 | 5.37 |
| FDA1 $\tau=50$ | 1.93 | **1.67** | 2.40 | 4.47 | 5.10 | 5.43 |
| FDA3 $\tau=10$ | 1.80 | **1.20** | 3.20 | 4.83 | 4.37 | 5.60 |
| FDA3 $\tau=25$ | 1.67 | **1.33** | 3.00 | 4.37 | 5.03 | 5.60 |
| FDA3 $\tau=50$ | 2.03 | **1.60** | 2.37 | 4.97 | 5.07 | 4.97 |
| FDA4 $\tau=10$ | **1.00** | 2.00 | 3.30 | 4.20 | 6.00 | 4.50 |

在 FDA1 $\tau=10$ 与 FDA4 $\tau=10$ 上，$\alpha=1.0$ **平均秩第一**；在 FDA4 $\tau=10$ 上 $\alpha=1.0$ 的平均秩为 1.00（30 个种子中排名第一的次数最多），与 PPS 的秩差 1.00 小于 $CD$。所有 9 个单元格的 Friedman 检验均显著（$\chi^2_F$ 118.3–139.5，$p<10^{-23}$）。

### 5.7.4 由 v2.5 结果得到的机制结论

**事实层**：$\alpha=1.0$ 使 LDT 在 9/9 配置上显著改善，并把与 PPS 的差距压到统计不可区分或反超的程度。**因此，v2.4 中观察到的"LDT 显著劣于 PPS"这一现象，主要来源是外推步长（半步 vs 全步），而不是"矩迁移 vs 点预测"的范式差异。** 这是一个对原设计不利、但同时也说明范式本身并未失效的结论：一旦步长对齐，矩迁移只需要 $2D$ 个参数、单次响应耗时仍在 $10^{-2}$ s 量级，就能达到甚至（在三目标问题上）超过点预测。

**机制假说层（需后续实验验证）**：为什么 $\alpha=1.0$ 在**三目标**的 FDA4 上反而优于 PPS，而在二目标的 FDA3 上仍略逊？一个可能的原因是：FDA4 的 PS 移动涉及角度变量与距离变量的耦合，其中心运动不如 FDA1/FDA3 那样严格匀速直线平移，纯中心外推（PPS）的"命中优势"被削弱；而 LDT 迁移种群携带的方差信息（逐维 $\sigma_e$）与 30% 随机移民提供了更宽的初始覆盖，在段内预算紧张（$\tau=10$ 仅 10 代）时更不容易落空。这与"$\alpha=1.0$ 在 FDA4 $\tau=10$ 的改进幅度（13.6%）远大于 $\tau=50$（1.1%）"的形态一致。但**严格的因果证明需要额外的诊断实验**（如记录预测中心的逐段误差、段内收敛轨迹、以及方差命中率），本文只报告相关性证据。

## 5.8 敏感性实验（$\tau=25$；除扫描参数外 $N=100$、$n_t=20$）

### 5.8.1 种群规模 $N$ 与变化平滑度参数 $n_t$

表 5-10：$N$ 与 $n_t$ 的敏感性（三问题稳态 MIGD 的跨问题均值，$n=30$）。

| 参数 | 取值 | 跨问题均值 MIGD |
|---|---|---|
| $N$ | 50 | 0.0529 |
| $N$ | **100（默认）** | **0.0352** |
| $N$ | 200 | 0.0239 |
| $n_t$ | 10 | 0.0376 |
| $n_t$ | **20（默认）** | **0.0352** |
| $n_t$ | 40 | 0.0330 |

$N$ 从 100 翻倍到 200 使均值 MIGD 降低约 32%（预算翻倍的正常收益）；$n_t$ 越大（每步时间增量越小、PS 移动越平滑）越容易跟踪，但影响相对较小（±7%）。这一结果与 2.3.2 节的渐变环境定义一致：$n_t$ 增大等价于减小单步变化幅度 $d_k$，使"相邻环境分布迁移"的假设更成立。

![图 7-1 种群规模 N 敏感性](../../LDT-NSGA-II/results/figures_thesis/fig7_1_sensitivity_N.png)

图 7-1：种群规模 $N \in \{50,100,200\}$ 的敏感性（$\tau=25$，三问题稳态 MIGD 均值 ± 标准差）。MIGD 随 $N$ 增大持续改善，$100 \to 200$ 段的边际收益小于 $50 \to 100$ 段。

### 5.8.2 迁移比例 $\rho$ 与趋势系数 $\alpha$

表 5-11：$\rho$ 与 $\alpha$ 的逐问题稳态 MIGD（$n=30$）。

| 参数 | 取值 | FDA1 | FDA3 | FDA4 |
|---|---|---|---|---|
| $\rho$ | 0.25 | 0.00767 | 0.00826 | 0.09297 |
| $\rho$ | **0.50（默认）** | **0.00741** | **0.00794** | **0.09039** |
| $\rho$ | 0.75 | 0.00735 | 0.00779 | 0.08890 |
| $\alpha$ | 0（= 消融 A） | 0.00884 | 0.00943 | 0.10106 |
| $\alpha$ | 0.25 | 0.00826 | 0.00879 | 0.09568 |
| $\alpha$ | **0.50（默认）** | 0.00741 | 0.00794 | 0.09039 |
| $\alpha$ | 0.75 | 0.00656 | 0.00702 | 0.08685 |
| $\alpha$ | **1.00（全步外推）** | **0.00621** | **0.00664** | **0.08583** |

**这是一个必须如实报告的、对默认设计不利的发现**：MIGD 随 $\alpha$ 单调下降，$\alpha = 1.0$ 在三个问题上一致最优，相对默认 $\alpha = 0.5$ 改善 **16%（FDA1）、16%（FDA3）、5%（FDA4）**。这与 PPS（本质即全步中心外推）在二目标问题上最强的结果互相印证——在 FDA 这类匀速规则平移问题上，外推越足越准，半步（$\alpha=0.5$）系统性欠外推。

$\rho$ 的影响较小且方向一致：迁移比例 0.75 略优于默认 0.5（约 1–2%）。作者不回避：**若仅以本组基准为目标，默认参数应改为 $\alpha = 1.0$**；保留 $\alpha = 0.5$ 的理由只能是"在非线性/变加速变化下全步外推可能过冲"的稳健性考量，而该考量**尚未被本文实验覆盖**（见第 7 章展望）。v2.5 的主比较实验（5.7 节）已经为这一参数修改提供了完整的证据。

## 5.9 补充实验

### 5.9.1 DMOP1–DMOP3（项目变体，$\tau=25$）

表 5-12：DMOP 族（项目变体）稳态 MIGD（均值 ± 标准差，$n=30$；加粗为该问题最优）。

| 算法 | DMOP1 | DMOP2 | DMOP3（静态） |
|---|---|---|---|
| nsga2_static | 0.0101 ± 0.0012 | 0.0097 ± 0.0011 | **0.1915 ± 0.0002** |
| dnsga2_inspired | 0.0104 ± 0.0011 | 0.0098 ± 0.0011 | **0.1915 ± 0.0002** |
| ar1_predictive | 0.0079 ± 0.0012 | 0.0073 ± 0.0012 | 0.1920 ± 0.0002 |
| gmm_inspired | 0.0102 ± 0.0011 | 0.0097 ± 0.0011 | 0.1922 ± 0.0006 |
| **pps** | **0.0074 ± 0.0011** | **0.0069 ± 0.0011** | 0.1920 ± 0.0002 |
| ldt_nsga2（$\alpha=0.5$） | 0.0088 ± 0.0011 | 0.0082 ± 0.0011 | 0.1921 ± 0.0010 |

DMOP1/DMOP2 上的排序与 FDA1/FDA3 完全同型：PPS 最优、AR1 次优、LDT 第三且显著优于三个非预测基线（对 static/dnsga2/gmm 的 $p_{adj}$ 均为 $9.31\times10^{-9}$）。这说明 5.2/5.3 的结论对这两个二目标问题具有一致性。

**DMOP3 上六个算法的稳态 MIGD 全部落在 0.1915–0.1922（极差约 0.4%），没有任何响应策略带来可检测收益**——这与已知问题 K2 完全吻合：当前 DMOP3 实现中 $G(t)$ 算而未用、目标与时间无关，实为**静态问题**，预测/迁移策略无变化可利用。DMOP3 的结果因此为 K2 提供了直接的**实证佐证**（而不只是代码审查推断）。各算法在变化后首代 MIGD 上的差异（static 首代即 0.1915，LDT/GMM 首代约 0.06）只反映它们对"伪变化信号"的不同初始反应，经数代恢复后稳态等价，**不构成性能优势**。

### 5.9.2 权威 dMOP1/dMOP2/dMOP3（v2.5 新增）

（本节数据由 v2.5 阶段二生成，与 5.7 节同批同协议；三个问题的权威定义与解析 Pareto 集/前沿见 3.9.2 节，并由 `scripts/check_dmop_authoritative.py` 验证"解析 PS 落在解析 PF 上"的偏差为 0.000e+00。）

DMOPSUMMARYPLACEHOLDER

## 5.10 结果解读：事实层与机制假说分层

### 5.10.1 事实层（有统计检验支持）

1. **v2.4（$\alpha=0.5$）**：在二目标规则平移问题（FDA1/FDA3/DMOP1/DMOP2）上排序稳定为 PPS > AR1 > **LDT** > GMM ≈ static ≈ dnsga2；LDT 对后三类基线的优势在快速变化档（$\tau=10$）尤其明显且统计显著；与 PPS 的差距随 $\tau$ 增大从 27–34% 收敛到 5.3% 以内。
2. **v2.4**：PPS 在全部 9 个主比较组合与 DMOP1/DMOP2 上稳态 MIGD 均值第一，LDT 稳居第二；Wilcoxon-Holm 下 LDT 对 PPS 在 9/9 组合显著劣后；但 Friedman–Nemenyi 整体排序在 FDA1 $\tau=10$ 与 FDA4 全部三档（共 4 个组合）上二者不可区分（秩差 0.73–1.30 $< CD = 1.37668$）。
3. **v2.4**：FDA4 $\tau=10$ 上 LDT 的种子间标准差更小（0.0047 对 PPS 0.0075），运行稳定性占优；该优势在 $\tau=25/50$ 上同样存在（0.0039 对 0.0046；0.0017 对 0.0019）。
4. **组件贡献**：分布迁移 ≫ 趋势外推 ≫ 旧种群保留（后者无一致可检测增益）。
5. **默认 $\alpha=0.5$ 被敏感性实验证伪为非最优**，$\alpha=1.0$ 在三个问题上一致更优。
6. **v2.5（$\alpha=1.0$）**：LDT 在 12 个"问题 × $\tau$"配置中的 **11 个**上相对 $\alpha=0.5$ 更优；与 PPS 的差距缩小到 $-8.5\%\sim+3.4\%$，其中 4 个配置统计不可区分、**2 个配置显著反超**（FDA4 $\tau=10$ 与 $\tau=25$）。**因此 v2.4 中"LDT 显著劣于 PPS"主要归因于外推步长而非迁移范式。**
7. **v2.5 的反例（结论三）**：唯一例外是权威 DMOP1A（**Type III：解集静止、仅前沿旋转**），该问题上 $\alpha=0.5$ 以 0.0251 对 0.0570 **大幅领先** $\alpha=1.0$（$p=1.86\times10^{-9}$）。**这说明全步外推在解集不移动时会把中心增量的采样噪声放大成系统性偏移；外推步长应随"解集是否在移动"自适应**（见第 6 章结论三）。
8. **GMM 迁移开销为同配置 LDT 的 4.4–11 倍且性能不占优**；项目变体 DMOP3 六算法等价，实证确认其时间无关性（K2）。

### 5.10.2 机制假说层（作者解释，需后续实验验证）

FDA1/FDA3/DMOP1/DMOP2 的 PS 做整体、匀速、低维平移，是"一阶差分全步外推"的理想场景，PPS 因此占优；LDT 的半步外推在该场景系统性欠外推。把 $\alpha$ 提到 1.0 后，二目标问题上的差距随之消失（FDA1）或大幅收窄（FDA3 从 33.9% 降到 3.4%）——这是"步长假说"而非"范式假说"的直接证据。

LDT 在**三目标** FDA4 上反而显著优于 PPS 的**可能原因**是：三目标旋转球面 PS 的中心运动不再与简单平移对齐，纯中心外推的命中优势被削弱；而 LDT 迁移种群携带的方差信息（逐维 $\sigma_e$）与 30% 随机移民提供了更宽的初始覆盖，在段内预算紧张（$\tau=10$ 仅 10 代）时更不容易落空，同时带来更小的运行间方差。该解释与"LDT 变化后首代 MIGD 仍不如 PPS，但稳态差距远小于二目标问题，且改善幅度随 $\tau$ 减小而增大"的形态一致。但**严格的因果证明需要额外的诊断实验**（记录预测中心误差、方差命中率与段内收敛轨迹），本文只报告**相关性证据**，不宣称因果。

### 5.10.3 适用边界（明确不宣称的内容）

- 本文的全部结论基于 **FDA1/FDA3/FDA4 与 dMOP1–dMOP3** 这六个问题族、$M \in \{2,3\}$、$D \le 12$、$\tau \in \{10,25,50\}$。**不能**外推到高维目标（$M \ge 4$）、大规模决策变量、约束动态优化或真实工程问题。
- 三个预测/迁移基线均为**简化项目实现**（3.8 节的命名声明），本文的结论**不能**解读为对 PPS [4]、DNSGA-II [3] 或 KTMM 类方法原始版本的性能评定。特别是 `pps` 缺少 PPS 原文的方差/形预测项，`gmm_inspired` 的随机移民比例只有 10%（另两个预测基线为 20%）。
- 分布迁移与随机移民的收益在**本组基准**（PS 光滑、单峰、整体平移为主）上得到验证；对多峰、欺骗性、剧烈突变或非平稳环境，本文证据不足（DMOP3A 提供了一组多峰前沿的部分证据，见 5.9.2 节）。
- 参考前沿为**项目自定义确定性采样**（FDA4 去重后 507 点而非 PlatEMO 的 `UniformPoint` 规则），因此 IGD 的**绝对值**不宜与其他文献直接横向比较；算法间的**相对**比较不受影响。

## 5.11 已知问题与不确定性清单（K1–K26）

本课题以清单形式完整披露全部已知问题与不确定项。K1–K12 为 v2.4 阶段的处置记录；K13–K17 为源码审计发现；K18–K21 为论文写作期间的数值复核发现；K22–K26 为 v2.5 权威 dMOP 复核与数据补齐过程中发现并处置的问题。

| 编号 | 问题 | 影响 | 处置 |
|---|---|---|---|
| **K1** | AR1/GMM/PPS 的 `immigrant_ratio` 原为**死参数**：实现顺序为"平移/采样后恰好 $N$ 行 + 末尾 append 移民 + `_resize` 截断前 $N$ 行"，追加的随机移民被全部截断丢弃 | 修复前三个基线实为"纯预测、无随机移民" | **已修复（v2.3）**：改为"先预留 $n_{imm}$ 个名额再拼接"（`ar_nsga2.py` / `gmm_transfer_nsga2.py` / `pps_nsga2.py`），新增 `tests/test_k1_immigrant_fix.py` 3 个回归测试；受影响的 1080 次运行（主比较 810 + 补充 270）已全部重跑，v2.2 旧数据归档于 `results/raw_legacy_v22_k1/`。**第 5 章全部基线数据均为修复后数据。** |
| **K2** | DMOP3 目标函数与时间无关（$G(t)$ 算而未用），与文件头自述矛盾 | 该问题实为静态，其上算法比较无动态区分度 | **已处置（v2.4）**：死代码 $G = \lvert\sin(0.5\pi t)\rvert$ 删除；docstring 更正为"项目自定义静态欺骗问题"；**5.9.1 节六算法等价的结果已实证此点**；v2.5 另提供权威 dMOP3A（Type I，PS 时变）作为真正的动态多峰问题 |
| **K3** | DMOP2 为项目自定义实现，与权威 dMOP2 不一致 | 引用若不注明会误导 | **已核实并处置（v2.4）**：与权威定义逐式对照确认不同（系数 1 而非 9、$g$ 含 $+G$、指数固定为 2、无 $H(t)$），已写入 docstring 并在 3.9.3 节如实标注；v2.5 另实现权威 DMOP2A（Type II） |
| **K4** | `MomentTransferStrategy` 存在归一化/真实坐标混用 | 仅在显式传入该策略时触发，默认实验未受影响 | **已修复（v2.4）**：全程归一化空间计算 + 单次反归一化，新增 2 个回归测试；0 次实验重跑 |
| **K5** | `evaluation/hv.py` 的无 pymoo 回退 `_hv_exact` 容斥符号有 bug | 仅在 pymoo 缺失且 $n \le 20$ 时触达，正常实验路径从未执行 | **已修复（v2.4）**：三处错误（奇数子集符号、交集应取分量 `max`、循环多算空集项）修复后与 pymoo HV 逐位一致 |
| **K6** | `RunConfig.evals_per_env = 5000` 字段未被驱动器使用 | 仅文档性字段 | **已修复（v2.4）**：字段删除；真实预算 $= N \times \tau_t \times$ 段数 |
| **K7** | Friedman 排秩无并列平均处理；Nemenyi 的 `alpha` 参数不参与查表 | 连续型 MIGD 并列概率极低 | **已修复（v2.4）**：改用 `scipy.stats.rankdata` 平均秩；补充 $\alpha=0.10/0.05$ 双档 q 表；不支持的 $\alpha$ 显式抛错；新增 3 个回归测试，修复前后统计结果逐位相同 |
| **K8** | 根目录 `main.py` 是损坏的 IDE 模板文件（含冲突标记） | — | **已删除后重建（v2.4）**为一键流水线入口 |
| **K9** | AR1/PPS/DNSGA2 的 `boundary_clip_ratio` 硬编码 0.0 | 该列为诊断字段，`analyze.py` 不消费 | **已修复（v2.4）**：四算法均改为真实计算"贴界个体占比"；顺带修复 `ARNSGA2` 接收但未存储 `n_probes` 的潜在 `AttributeError` |
| **K10** | FDA4 参考前沿去重后点数不严格等于 500 | IGD 绝对值与"恰好 500 点"的叙述略有出入 | **已核实（v2.4）**：23×23 = 529 去重后实测 **507 点**，文档与测试统一该数字 |
| **K11** | PlatEMO 是否确实没有 DMOP1/2/3.m 未能在线复核 | 归属声明依赖作者自述 | **已复核（v2.4，并在 v2.5 再次独立复核）**：经 GitHub git-tree API 对 PlatEMO master（commit `d25e65d`）全树递归检索，确认无 DMOP1/2/3.m |
| **K12** | 实验报告"单文件自包含"表述与外链 `assets` 不符 | 部署需带 assets 目录 | **已修复（v2.4）**：构建时内联 CSS 与 ECharts，生成的 `index.html` 为单文件自包含（1,152,616 字节） |
| **K13** | **设计不确定（未修复，不影响结论）**：多源融合开启时，趋势项用融合后中心 $\mu_{\text{eff}}$ 减**未融合**的历史中心 $\mu_{prev}$，两量口径不一致 | 默认实验未开启多源融合，不影响第 5 章任何结论 | **如实保留**：由于 $\mu_{\text{eff}} - \mu_{prev} = (\mu_{\text{eff}} - \mu_e) + (\mu_e - \mu_{prev})$，该写法等价于在趋势项上叠加融合校正项，**不构成数值错误**；但是否为设计意图，作者不确定 |
| **K14** | `run_supplementary.py` 自身默认只跑 `DMOP1,DMOP2`（360 次），`main.py` 未覆盖该默认 | 按默认命令无法复现论文的 540 次补充实验 | **已披露（v2.5）**：4.7 节明确要求显式加 `--problems DMOP1,DMOP2,DMOP3`；列入第 7 章待修项 |
| **K15** | `build_report.py` 生成的 HTML 报告中**叙述文字**仍是 v1 硬编码（"4 算法、536,700 行、平均 MIGD 0.0378"），与 v2.4 的 6 算法 / 1,134,000 行不符 | 报告表格/图表数据来自 JSON（正确），但叙述文字不可引用 | **已披露（v2.5）**：本文全部数字取自逐代原始 CSV 的独立重算，不依赖该报告叙述；列入第 7 章待修项 |
| **K16** | `configs/*.yaml` 为过时镜像（4 算法、`code_version: v1.0`、$N \in [50,100,150]$、70/20/10 等） | 四个 YAML 未被任何脚本读取 | **已披露（v2.5）**：4.8 节声明不可作为实验配置依据；列入第 7 章待清理项 |
| **K17** | `README.md` 的 `sigma_min = 0.001*(u-l)`（真实尺度）与代码（归一化单位 $10^{-3}$）不符 | 描述性错误 | **已披露（v2.5）**：论文统一采用代码口径（归一化 $10^{-3}$） |
| **K18** | 既有材料中 DMOP1 上 static 的稳态 MIGD 记为 **0.0102**，实测为 **0.010149**（四位小数应为 0.0101） | 单格四舍五入越界（阈值 0.01015，超出 1.3×10⁻⁶） | **已修正（v2.5）**：表 5-12 采用 0.0101；同表其余 17 项复核一致 |
| **K19** | 既有材料称消融变体 C（$\kappa=0$）与 full "统计等价、无可检测差异" | 实测 FDA1 $p=0.0175$、FDA3 $p=0.0384$ 的**原始** $p$ 值小于 0.05 | **已修正（v2.5）**：5.6 节改为如实报告原始 $p$ 值、方向不一致性与 Holm 校正后不显著，结论改为"证据不支持把 $\kappa$ 作为默认机制" |
| **K20** | 论文正文的版本标签（v2.3/v2.4）与落盘工件的 `code_version` 不一致 | 引用数据时易混淆 | **已披露（v2.5）**：全部 1,134,000 行 CSV 与 3240 条 `run_meta` 记录实际标为 **`v2.2`**（由 commit `50adf19c5c5a` 生成），而 `DynamicDriver.code_version` 当前为 `"v2.4"`；4.1 节给出"数据生成提交 vs 代码分析提交"的双重标注 |
| **K21** | "3780 次**独立**运行"的措辞不准确 | 3780 是**运行记录**数；按 distinct (algorithm, problem, $\tau$, seed) 计为 **2430** 个不同组合（计入 $N$、$n_t$ 后为 2790 个配置），其中默认 LDT 配置在 comparison / ablation / $N{=}100$ / $n_t{=}20$ / $\rho{=}0.5$ / $\alpha{=}0.5$ 六个文件中被重复执行 | **已修正（v2.5）**：全文统一表述为"**3780 次运行记录**"，并在 5.1.1 节明确区分运行记录数、配置数与不同组合数 |
| **K22** | 权威 **dMOP1 的类型与公式写错**：早期草稿记为 Type I、$g$ 含 $G(t)$、前沿为 $\sqrt{\cdot}$ 型 | 目标函数、参考前沿与问题类型三者同时错误，使该问题上的 IGD 计算失去意义 | **已修复（v2.5）**：经四份独立复现文献复核，权威 dMOP1 实为 **Type III**（解集静止、前沿随 $H(t)$ 旋转），$g = 1+9\sum_{i\ge2}x_i^2$（**不含** $G(t)$），前沿为 $1-u^{H(t)}$。已重写 `problems/dmop1a.py` 并更正 `true_fronts.py` 的参考前沿；首版生成的 6 个 CSV 移出统计目录、归档于 `results/raw_v25_badauthoritative/` 后**重跑** |
| **K23** | 权威 **dMOP3 的公式写错**：早期实现误用了项目变体的欺骗性前沿（含 $f_1\sin(10\pi f_1)$ 项）与 $H(t)$ 解集位置 | 同上，问题定义与文献不符 | **已修复（v2.5）**：更正为 $f_1 = x_r$（$r$ 按环境轮换）、$g = 1+9\sum_{i\ne r}(x_i-G(t))^2$、前沿 $1-\sqrt{f_1}$（静止）。同批归档并重跑。**同时如实记录一处文献分歧**：dMOP3 的 $G(t)$ 在两份复现文献中分别为 $\sin(\pi t/2)$ 与 $\lvert\sin(\pi t/2)\rvert$，本实现采用后者并已在 3.9.2 节显式声明 |
| **K24** | v2.5 阶段二最初**缺少 LDT 的 $\alpha=0.5$ 变体**在权威 dMOP 上的数据（`run_v25_alpha.py` 的默认算法列表用 `ldt_alpha1` 替代了 `ldt_nsga2`） | 若不补齐，只能比较 $\alpha=1.0$ 与基线，无法回答"在这一族问题上外推步长是否同样重要" | **已补齐（v2.5）**：补跑 3 问题 × 1 算法 × 30 种子 = **90 次**运行。正是这 90 次数据揭示了 DMOP1A 上的**反例**（见表 5-19 与第 6 章结论三） |
| **K25** | `_row_metrics` 的时变前沿判定最初未包含 `DMOP1A` | Type III 的前沿同样随 $t$ 变化；若沿用 $t=0$ 的前沿会**系统性高估** IGD | **已随 K22 一并修复（v2.5）** |
| **K26** | `get_true_pf` 的默认目标数判定原按问题名硬编码（`M = 2 if name in (FDA1,FDA3,DMOP1,DMOP2,DMOP3) else 3`） | 新增问题类时必须手工维护该列表，易漏 | **已重构（v2.5）**：改为 `M = 3 if name in ("FDA4","FDA5") else 2`；对原有 6 个问题的结果不变 |

**v2.4 处置汇总**：K1 已于 v2.3 修复并重跑；K2–K12 于 v2.4 全部处置（K4/K5/K7/K9 代码修复且不改变实验行为、0 次重跑，K2/K3/K10/K11 核实与文档更正，K6/K8/K12 清理与构建修复），新增 `tests/test_v24_fixes.py` 15 个回归测试，测试总数 63 个全部通过。**v2.5 处置汇总**：K13–K26 为新增披露，其中 K18/K19/K21 为数值与表述修正，K14–K17/K20 为如实披露，**K22/K23/K25 为实质性的问题定义更正并触发实验重跑**，K24 为数据补齐，K26 为代码重构；另新增 `LDTNSGA2Alpha1`、权威 dMOP1A/2A/3A、`experiments/run_v25_alpha.py`、`scripts/check_dmop_authoritative.py`、`scripts/thesis_analysis.py`、`scripts/generate_dmop_section.py`。

## 5.12 本章小结

本章报告了 6030 次**运行记录**的完整结果。v2.4 的结论是"PPS 第一、LDT（$\alpha=0.5$）稳居第二且显著优于三个非预测基线"；v2.5 的关键对照实验把这一结论推进了两步：

1. **步长归因**：把外推步长从半步改为全步（$\alpha=1.0$）后，LDT 在 12 个"问题 × $\tau$"配置中的 **11 个**上更优，与 PPS 的差距从 $+1.4\%\sim+33.9\%$ 缩小到 $-8.5\%\sim+3.4\%$，其中 4 个配置统计不可区分、2 个配置（FDA4 的两档）**显著反超**。由此可以给出比 v2.4 更清晰的机制结论：v2.4 中观察到的劣势主要来自**外推步长**，而非"矩迁移 vs 点预测"的范式差异。
2. **反例与自适应规则**：唯一的例外是权威 **DMOP1A（Type III，解集静止、仅前沿旋转）**——该问题上 $\alpha=0.5$ 反以 **127%** 的优势领先 $\alpha=1.0$。这说明解集不移动时中心增量几乎只含噪声，全步外推会把噪声放大为系统性位置偏移。由此提出**"外推步长应随解集是否显著移动而自适应"**的设计规则，并用 12 个配置从两端验证了该规则。

消融实验确认分布迁移是贡献最大的组件、趋势外推次之、旧种群保留无一致增益；敏感性实验确认 $\alpha$ 与 $N$ 是影响最大的两个参数。本章同时以 **K1–K26** 清单完整披露了实现缺陷、版本戳漂移、口径措辞、数值舍入以及**两处权威 dMOP 定义错误**等全部已知问题。第 6 章作出结论性总结。


\newpage

# 第 6 章 结论与展望

## 6.1 主要结论

本文面向渐变环境下的动态多目标优化问题，设计并实现了一种轻量级的分布迁移算法 LDT-NSGA-II，构建了完全可复现的实验系统，并通过两轮共 **6030 次运行**的受控实验得到了如下结论。

**结论一（算法设计层面）**：把环境响应机制与 NSGA-II 内核**正交解耦**，并在**统一的归一化决策空间**中完成全部分布计算（特征提取、趋势外推、采样、截断），是保证算法可审查、可比、可复现的关键设计。LDT-NSGA-II 只需上一环境的逐维均值与标准差（$2D$ 个标量）与一条确定性外推规则即可生成新环境的初始种群，单次响应复杂度 $O(MN^2 + ND)$，不需要保存多代完整历史种群、不需要拟合时序模型或概率模型。消融实验确认其三种成分的贡献排序为：**分布迁移 ≫ 趋势外推 ≫ 旧种群保留**；其中"旧种群保留"在本组基准上没有一致可检测的增益。

**结论二（机制归因层面，本文最重要的发现）**：v2.4 的实验中，默认半步外推（$\alpha=0.5$）的 LDT 在全部 9 个主比较配置上显著劣于简化版 PPS（中心全步外推）。本文没有止步于"预测式方法更强"这一表层解释，而是设计了**"外推步长 × 迁移范式"解耦的对照实验**：把 $\alpha$ 改为 1.0 后，LDT 在 12 个配置中的 **11 个上更优**，与 PPS 的差距从 $+1.4\%\sim+33.9\%$ 缩小到 $-8.5\%\sim+3.4\%$，其中 4 个配置统计不可区分、2 个配置（FDA4 的 $\tau=10$ 与 $\tau=25$ 两档）**显著反超**。**由此可以明确：v2.4 观察到的劣势主要来自外推步长，而不是"矩迁移 vs 点预测"的范式差异。** 这是一个对原默认设计不利、但同时也说明矩迁移范式本身并未失效的结论——步长对齐之后，只用 $2D$ 个分布参数的轻量方法就能达到甚至超过点预测基线。

**结论三（机制归因层面，第二个发现：一个明确的反例）**：全步外推并非在任何场景都更好。在权威 **DMOP1A（Type III：Pareto 最优解集静止、仅前沿随 $H(t)$ 旋转）** 上，$\alpha=0.5$ 以 **0.0251 对 0.0570** 大幅领先 $\alpha=1.0$（$p=1.86\times10^{-9}$，相对改进 **127%**），并且是六个算法中该问题上唯一显著跳出"0.057 平台"的配置（其余五个算法与 $\alpha=1.0$ 全部落在 0.0570–0.0583 的窄带内）。其机理清晰：**解集不动时，中心增量 $\mu_e-\mu_{prev}$ 几乎只含采样噪声，把它整额外推出去会把噪声放大成系统性的位置偏移；半步相当于对噪声做了一次收缩。**

把结论二与结论三合起来，得到一个**可操作的设计规则**，这也是本文在方法层面的实质贡献：

> **外推步长应当随"Pareto 解集是否在移动"自适应：解集显著移动时取 $\alpha \to 1$（全步），解集近似静止时取 $\alpha \to 0.5$ 甚至 $\to 0$。**

本文的数据已经把这条规则的**两个端点都验证过了**（11 个移动解集的配置支持全步、1 个静止解集的配置支持半步），只差把判据接上——而这个判据所需的量（中心位移范数、或共享探针的变化幅度 $D_k$）LDT 本身已经在计算。**这也从实验上验证了原设计中选择半步的保守考虑是合理的，只是它被用在了错误的场景（应当只用于解集静止或变化不确定时）。**

**结论三（精度—开销谱系层面）**：在全协方差高斯混合迁移（GMM 基线）与中心点预测（PPS/AR(1)）之间，LDT 所占据的位置是"**接近最优预测方法的恢复性能 + 接近零的响应开销**"。GMM 的单次响应耗时约为同配置 LDT 的 4.4–11 倍，而 MIGD 并不占优；LDT 与 AR(1)、PPS 同处 $10^{-2}$ s 量级。这为算力受限、需要在线响应的场景提供了一条可行的折中路线。

**结论四（实验方法论层面）**：口径统一、可审计的实验系统对结论可信度的影响被严重低估。本文通过三件事提升了结论的可信度：(a) **双实现互检**——论文汇总脚本与项目原有分析脚本对 54 个主比较单元格独立计算、结果逐位一致；(b) **独立数值复核**——一个不读取既有汇总、只从逐代原始 CSV 从零重算的复核流程，逐项比对了 269 个数值断言，发现并修正了 4 处不一致；(c) **完整披露**——以 K1–K26 清单列出全部已知问题，包括一处已被修复并重跑 1080 次运行的基线"死参数"缺陷、两处测试问题归属错误、以及一处**在 v2.5 复核中发现并更正的权威 dMOP 定义错误**（权威 dMOP1 的正确类型是 **Type III**——解集静止、前沿时变，而非早期草稿所写的 Type I）。

**结论五（适用边界）**：上述结论的有效范围是本文所用的六个问题族（FDA1/FDA3/FDA4 与 dMOP1–dMOP3）、目标数 $M \in \{2,3\}$、决策维数 $D \le 12$、变化周期 $\tau \in \{10,25,50\}$。本文**不宣称**结论可外推到高维目标（$M \ge 4$）、大规模决策变量、约束动态优化或真实工程问题；三个预测/迁移基线均为简化项目实现，本文结果**不构成**对 PPS [4]、DNSGA-II [3] 或 KTMM [14] 原始版本的性能评定。

## 6.2 本文的创新点

1. **轻量级矩迁移机制**：提出在统一归一化空间中以"逐维均值 + 逐维标准差"两个矩特征、仅依赖相邻环境、无需训练任何模型的分布迁移机制，并给出完整的组件化设计与消融证据（结论一）。
2. **步长—范式解耦的对照实验，并发现"外推步长应随解集运动自适应"**：设计并完成把"外推步长"与"迁移范式"分离的受控实验，用同一随机种子下的配对检验把 LDT 与最优预测基线的差距重新归因到步长（结论二）；进一步在权威 dMOP1（Type III，解集静止）上发现全步外推的**反例**（结论三），据此提出可操作的自适应步长规则，并用 12 个配置从两端验证了该规则。据作者所知，这一"步长—范式解耦 + 运动自适应"的组合归因在同类轻量级 DMOEA 的研究中并未被系统做过。
3. **可复现且如实披露的实验系统**：构建口径写死在代码中、带 63 个单元测试、运行级元数据与断点续跑的实验系统，并以已知问题清单的形式完整披露实现缺陷、口径措辞与文献归属错误；同时用"双实现互检 + 独立数值复核"验证论文表格中的每一个数字（结论四）。

## 6.3 局限性

本文的局限可按"证据强度"分为三类，均如实列出。

**第一类：证据充分但适用范围受限。**

- 问题集有限。六个问题族中，四类动态性（Type I/II/III）虽已覆盖，但 PS 的运动形式以"整体平移"与"前沿旋转"为主，缺少 **PS 做异质、非线性、变加速运动**的问题（如 irregular-change 基准、JY 系列 [17] 的更难配置）。因此"全步外推最优"这一结论**可能只在规则运动场景成立**。
- 目标数与维数偏低（$M \le 3$、$D \le 12$）。

**第二类：证据不足或存在混杂因素。**

- 三个基线的**随机移民比例不统一**：`ar1_predictive` 与 `pps` 为 20%，而 `gmm_inspired` 仅 10%。这使 GMM 的性能读数中混杂了"移民更少"这一因素，本文未做统一（统一会使其进一步偏离原始实现）。
- "旧种群保留无增益"这一结论基于**四个变体、三个问题**，且变体 C 与 full 的差异方向在问题之间不一致（FDA1 上 full 更好、FDA3/FDA4 上 C 更好）；虽然 Holm 校正后不显著，但**不能**据此断言该组件在一般情形下无用——它可能只在本组光滑平移问题上冗余。
- $\rho=0.75$ 略优于默认 $0.5$（约 1–2%）这一观察未做配对检验，样本只覆盖三个问题。

**第三类：已披露的实现与文档缺陷（不影响第 5 章的结论，但影响可复现体验）。**

- K14：`run_supplementary.py` 自身默认只跑 `DMOP1,DMOP2`，按默认命令无法复现论文的 540 次补充实验。
- K15：`build_report.py` 生成的 HTML 报告中叙述文字仍是 v1 硬编码，与 v2.4 的 6 算法 / 1,134,000 行不符（表格数据正确）。
- K16：`configs/*.yaml` 为过时镜像且未被任何脚本读取。
- K17：`README.md` 的 $\sigma_{\min}$ 描述（真实尺度）与代码（归一化尺度）不符。
- K20：v2.4 全部原始数据的 `code_version` 戳为 `v2.2`，与代码当前值 `v2.4` 不一致。
- K13：多源融合开启时趋势项的量纲口径不一致（是否为设计意图，作者不确定）。

**第四类：方法学上的未验证项。**

- **浮点级跨平台复现未验证**：未做跨机器、跨 BLAS 实现的逐位比对。
- **机制假说未做因果验证**：5.10.2 节对"LDT 在 FDA4 上反超 PPS"的解释（宽初始覆盖 + 随机移民兜底在段内预算紧张时更抗落空）只是相关性证据，缺少诊断实验的直接支持。

## 6.4 未来工作

按优先级（收益/成本比）排序：

1. **把 $\alpha=1.0$ 纳入正式协议并重做全部比较（高优先级、低成本）。** 本文的 v2.5 实验已经证明 $\alpha=1.0$ 在 9/9 配置上更优。下一步应把 LDT 的默认参数正式改为 $\alpha=1.0$，并把消融、敏感性、dMOP 补充实验全部以新默认重跑，给出统一的最终结论表。**这是本文已经铺好路、只差执行的一步。**
2. **在非规则变化问题上检验"全步外推"是否过冲，并把自适应步长实现出来（最高优先级）。** 本文的结论三已经给出了**一个反例**：在解集静止的 Type III 问题上全步外推是有害的（127% 的差距）。因此第 1 项"直接改成 $\alpha=1.0$"**不能**简单执行——正确的下一步是实现**自适应步长**：用 LDT 自身已在计算的量（中心位移范数 $\lVert\mu_e-\mu_{prev}\rVert$、或共享探针的变化幅度 $D_k$ 序列）判断"解集是否在显著移动"，据此在 $\alpha \in [0.5, 1.0]$ 之间取值（甚至允许 $\alpha=0$）。本文的数据已经验证了该规则的两个端点，实现与验证是紧接着的下一步。同时应引入 PS 做**变加速或非线性移动**的问题（irregular-change 基准、JY 系列 [17]、或自行构造分段变向的 FDA 变体），检验自适应规则能否同时保住规则问题上的收益与静止问题上的稳健性。
3. **机制诊断实验（中优先级）。** 逐段记录**预测中心误差** $\lVert\hat\mu - \mu^*\rVert$、**方差命中率**（迁移个体的 $\sigma$ 与真实 PS 宽度的比）以及段内收敛轨迹，把 5.10.2 节的因果假说变成可检验的命题；特别要在 FDA4 上回答"反超是否真的来自更宽的初始覆盖"。
4. **扩展到更高目标数与更紧预算（中优先级）。** 增加 $M \in \{4,5\}$ 的 DTLZ/FDA 变体与 $\tau=5$ 的极紧预算配置，检验矩迁移的对角协方差假设在目标数增加时是否失效（此时 PS 的移动方向可能与坐标轴不再对齐）。
5. **清理工程债务（低优先级、必要）。** 修复 K14（补充实验默认参数）、K15（报告叙述文字）、K16（过时 YAML）、K17（README 口径）、K20（版本戳统一）；把权威 dMOP1A/2A/3A 与 `thesis_analysis.py` 补入单元测试体系，使测试总数从 63 增至覆盖全部 v2.5 组件的规模。
6. **变化强度自适应的实证（低优先级）。** `adaptive_transfer` 目前默认关闭。应在敏感性实验中加入其开/关对比，验证 $[0.5, 2.0]$ 缩放上限的合理性，并检验它与第 2 项"自适应步长"能否协同。
7. **向真实工程问题延伸（长期）。** 选择一类确有渐变特性的工程动态多目标问题（如含负荷漂移的调度、含设备老化的参数整定），检验在标准基准上得到的结论是否仍然成立。这一步是本文结论从基准走向应用的必经环节，也是作者在 1.3 节明确声明尚未完成的部分。


\newpage

# 参考文献

> 著录格式遵循 GB/T 7714-2015 顺序编码制。全部条目的卷、期、页码与 DOI 均经 Crossref / 出版方页面 / 官方仓库核对（核对记录见 `_audit/bibliography.md` 与 `_audit/reference_check.md`）；**未能核实的字段在该条末尾以"未能核实"显式标注**，不做推测填补。

[1] DEB K, PRATAP A, AGARWAL S, et al. A fast and elitist multiobjective genetic algorithm: NSGA-II[J]. IEEE Transactions on Evolutionary Computation, 2002, 6(2): 182-197. DOI: 10.1109/4235.996017.

[2] FARINA M, DEB K, AMATO P. Dynamic multiobjective optimization problems: test cases, approximations, and applications[J]. IEEE Transactions on Evolutionary Computation, 2004, 8(5): 425-442. DOI: 10.1109/TEVC.2004.831456.

[3] DEB K, RAO N U B, KARTHIK S. Dynamic multi-objective optimization and decision-making using modified NSGA-II: a case study on hydro-thermal power scheduling[C]//Evolutionary Multi-Criterion Optimization (EMO 2007). Berlin: Springer, 2007: 803-817. (Lecture Notes in Computer Science, vol. 4403). DOI: 10.1007/978-3-540-70928-2_60.

[4] ZHOU A, JIN Y, ZHANG Q. A population prediction strategy for evolutionary dynamic multiobjective optimization[J]. IEEE Transactions on Cybernetics, 2014, 44(1): 40-53. DOI: 10.1109/TCYB.2013.2245892.

[5] GOH C K, TAN K C. A competitive-cooperative coevolutionary paradigm for dynamic multiobjective optimization[J]. IEEE Transactions on Evolutionary Computation, 2009, 13(1): 103-127. DOI: 10.1109/TEVC.2008.920671.（dMOP1–dMOP3 测试问题的出处）

[6] TIAN Y, CHENG R, ZHANG X, et al. PlatEMO: a MATLAB platform for evolutionary multi-objective optimization[J]. IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87. DOI: 10.1109/MCI.2017.2742868.（正式题名含栏目标记"[Educational Forum]"；本项目 FDA 问题实现的逐式对照基准：https://github.com/BIMK/PlatEMO ）

[7] BLANK J, DEB K. pymoo: multi-objective optimization in Python[J]. IEEE Access, 2020, 8: 89497-89509. DOI: 10.1109/ACCESS.2020.2990567.（NSGA-II 内核与 SBX/PM 算子的实现来源）

[8] BOSMAN P A N, THIERENS D. The balance between proximity and diversity in multiobjective evolutionary algorithms[J]. IEEE Transactions on Evolutionary Computation, 2003, 7(2): 174-188. DOI: 10.1109/TEVC.2003.810761.（IGD 指标的通行引用来源）

[9] ZITZLER E, THIELE L. Multiobjective evolutionary algorithms: a comparative case study and the strength Pareto approach[J]. IEEE Transactions on Evolutionary Computation, 1999, 3(4): 257-271. DOI: 10.1109/4235.797969.（超体积指标 HV）

[10] WILCOXON F. Individual comparisons by ranking methods[J]. Biometrics Bulletin, 1945, 1(6): 80-83. DOI: 10.2307/3001968.（Crossref 仅登记起始页 80；80-83 为通行引用范围）

[11] HOLM S. A simple sequentially rejective multiple test procedure[J]. Scandinavian Journal of Statistics, 1979, 6(2): 65-70.（该文无 DOI；卷期页码经多条独立文献表核对）

[12] FRIEDMAN M. The use of ranks to avoid the assumption of normality implicit in the analysis of variance[J]. Journal of the American Statistical Association, 1937, 32(200): 675-701. DOI: 10.1080/01621459.1937.10503522.

[13] DEMŠAR J. Statistical comparisons of classifiers over multiple data sets[J]. Journal of Machine Learning Research, 2006, 7: 1-30.（JMLR 无 DOI；见 https://www.jmlr.org/papers/v7/demsar06a.html 。本文的 Wilcoxon / Friedman + Nemenyi 使用规范与 Nemenyi $q_\alpha$ 表（Table 5）均取自该文）

[14] ZOU J, HOU Z, JIANG S, et al. Knowledge transfer with mixture model in dynamic multiobjective optimization[J]. IEEE Transactions on Evolutionary Computation, 2025, 29(5): 1517-1530. DOI: 10.1109/TEVC.2025.3566481.（即本文所称"KTMM"；已正式出版，非在线预印；第一作者为 Zou。本项目的 GMM 基线受其"全分布知识迁移"思想启发，但**不是**其复现）

[15] SRINIVAS N, DEB K. Multiobjective optimization using nondominated sorting in genetic algorithms[J]. Evolutionary Computation, 1994, 2(3): 221-248. DOI: 10.1162/evco.1994.2.3.221.（第一代 NSGA，绪论历史背景）

[16] HELBIG M, ENGELBRECHT A P. Benchmarks for dynamic multi-objective optimisation algorithms[J]. ACM Computing Surveys, 2014, 46(3): 1-39. DOI: 10.1145/2517649.（动态多目标优化基准族的参考综述；本文复核 dMOP1–dMOP3 定义时所依据的复现文献之一即出自该文件者组的章节）

[17] JIANG S, YANG S. Evolutionary dynamic multiobjective optimization: benchmarks and algorithm comparisons[J]. IEEE Transactions on Cybernetics, 2017, 47(1): 198-211. DOI: 10.1109/TCYB.2015.2510698.（即"JY"基准系列 JY1–JY10 的出处；本文列入第 6 章展望）

[18] JIANG S, ZOU J, YANG S, et al. Evolutionary dynamic multi-objective optimisation: a survey[J]. ACM Computing Surveys, 2023, 55(4): 1-47. DOI: 10.1145/3524495.（动态多目标优化的近期系统综述；未能核实 ACM 文章号，故按页码著录）

---

## 附：进一步阅读（未在正文中引用，供读者延伸）

以下两条同样经过著录核实，但与本文的具体论证无直接引用关系，故不占用正文引用编号。

- HELBIG M. Solving dynamic multi-objective optimisation problems using vector evaluated particle swarm optimisation[D]. Pretoria: University of Pretoria, 2012.（博士论文；**未能核实**：本次检索期间该机构库 PDF 无法下载（TLS 握手失败/超时），故其式 3.38–3.40 **未被直接阅读**，本文对 dMOP 定义的复核改以四份可获取的复现文献交叉确认）
- JIANG S, YANG S. A steady-state and generational evolutionary algorithm for dynamic multiobjective optimization[J]. IEEE Transactions on Evolutionary Computation, 2017, 21(1): 65-82. DOI: 10.1109/TEVC.2016.2574621.


\newpage

# 附录 A 零基础运行手册

**目标读者**：没有接触过本项目、也没有接触过进化算法的读者。按顺序操作即可，每一步都给出**预期结果**，以便判断是否走对了。

## A.1 第 0 步：安装 Python

安装 Python 3.9 或更高版本（本项目实测版本 **3.9.13**）。Windows 安装时务必勾选 **Add Python to PATH**。

验证：

```powershell
python --version
```

预期：显示 `Python 3.9.13`（或更高版本号）。

## A.2 第 1 步：获取项目并安装依赖

进入项目根目录（假设为 `LDT-NSGA-II`）：

```powershell
cd LDT-NSGA-II
python -m pip install -r requirements.txt
```

`requirements.txt` 全文：

```text
pymoo>=0.6
numpy
scipy
scikit-learn
matplotlib
pyyaml
tqdm
pytest
```

预期：pymoo、numpy、scipy、scikit-learn 等安装成功，无红色报错。

> **版本提示**：本项目实测 pymoo **0.6.1.5**。pymoo 0.6.x 没有 `algo.n_eval` 属性（真实计数器在 `algo.evaluator.n_eval`），项目已做兼容封装 `_algo_n_eval()`；若升级到 pymoo 大版本，请重新验证该封装。

## A.3 第 2 步：跑单元测试

```powershell
python -m pytest -q
```

预期：最后一行形如 `63 passed in ...`。**若有 failed，先停下排查环境，不要带病跑实验。**

## A.4 第 3 步：单次冒烟实验（约 1–3 分钟）

```powershell
python -m experiments.run_comparison --problems FDA3 --taus 25 --seeds 0 `
    --out results/_smoke --n-gen 50 --limit-seeds 1
```

预期：打印 6 个算法各跑 1 次；`results/_smoke/` 下新增 6 个 `comparison_*_FDA3_tau25.csv` 与一个 `run_meta.jsonl`。

> ⚠️ **务必使用独立的输出目录**。冒烟数据用了缩小预算（50 代而非 300 代），一旦混入 `results/raw` 就会污染正式统计。本项目实际发生过这类风险，故手册统一要求写入 `results/_smoke`。

## A.5 第 4 步：v2.4 四组实验

全量约数小时；**可随时中断、重新执行同一命令即从断点继续**。

```powershell
# 4a. 主比较实验（1620 次运行）
python -m experiments.run_comparison --problems FDA1,FDA3,FDA4 --taus 10,25,50 `
    --seeds 0-29 --out results/raw

# 4b. DMOP 补充实验（540 次运行）—— 注意必须显式指定三个问题
python -m experiments.run_supplementary --problems DMOP1,DMOP2,DMOP3 --seeds 0-29 `
    --out results/raw --algos nsga2_static,dnsga2_inspired,ar1_predictive,gmm_inspired,pps,ldt_nsga2

# 4c. 消融实验（360 次运行）
python -m experiments.run_ablation

# 4d. 敏感性实验（1260 次运行）
python -m experiments.run_sensitivity
```

预期：每个脚本结束时打印

```text
[SUMMARY] expected=... completed=... skipped=... failed=0
```

> ⚠️ **已知陷阱（K14）**：`run_supplementary.py` 自身的默认参数是 `--problems DMOP1,DMOP2`（只有 360 次），`main.py` 也未覆盖该默认值。**必须**按 4b 的写法显式给出三个问题，否则无法复现论文的 540 次。

## A.6 第 5 步：v2.5 两组关键实验

```powershell
# 5a. α=1.0 主比较（1620 次运行；写入独立目录 raw_v25）
python -m experiments.run_v25_alpha --problems FDA1,FDA3,FDA4 --taus 10,25,50 `
    --seeds 0-29 --out results/raw_v25 --phase comparison

# 5b. 权威 dMOP 补充（540 次运行）
python -m experiments.run_v25_alpha --problems DMOP1A,DMOP2A,DMOP3A --taus 25 `
    --seeds 0-29 --out results/raw_v25 --phase supplementary
```

> **并行加速（可选）**：脚本按 (问题, $\tau$) 组合可拆成 9 个独立进程并行运行（每个进程处理一个问题的一个 $\tau$、全部 30 个种子），在本项目使用的 32 逻辑核机器上把 5a 由约 5 小时缩短到约 0.5 小时：

```powershell
$py = "python"
$probs = @("FDA1","FDA3","FDA4"); $taus = @(10,25,50)
$procs = @()
foreach ($p in $probs) { foreach ($t in $taus) {
  $procs += Start-Process -FilePath $py -NoNewWindow -PassThru -WorkingDirectory (Get-Location) `
    -ArgumentList "-m","experiments.run_v25_alpha","--problems",$p,"--taus","$t",`
                  "--seeds","0-29","--out","results/raw_v25","--phase","comparison"
}}
$procs | Wait-Process
```

> ⚠️ `run_v25_alpha.py` 内置了保护：当 `--out` 指向 `results/raw` 或 `raw_legacy*` 时会**直接报错退出**，防止污染 v2.4 的冻结数据。

## A.7 第 6 步：权威 dMOP 自检与论文级统计

```powershell
# 6a. 解析 PS 是否精确落在解析 PF 上（应全部为 0.000e+00），并验证问题类型
python -m scripts.check_dmop_authoritative

# 6b. 从原始逐代 CSV 重新计算论文的全部表格与统计检验
python -m scripts.thesis_analysis
```

预期：6a 打印 `FINAL: ALL CHECKS PASSED`；6b 在 `results/summary_v25/` 下生成 `analysis_report.md`、`per_run_metrics.csv`、`summary_v24.csv`、`summary_v25.csv`、`wilcoxon_holm_*.csv`、`friedman_*.csv`、`alpha1_vs_alpha05.csv`、`recovery_summary.csv`、`integrity.json` 等文件。

## A.8 第 7 步：v2.4 冻结汇总与图表（可选）

```powershell
python analyze.py                      # 重新生成 results/summary/ 与 results/figures/
python -m scripts.make_thesis_figures  # 生成 results/figures_thesis/ 论文插图
python build_report.py                 # 生成 HTML 报告
```

> ⚠️ `analyze.py` 用 glob 读取 `results/raw/comparison_*.csv` 等模式。**不要在 `results/raw` 中存放 v2.5 的数据**（`ldt_alpha1` 等），否则 v2.4 的汇总会被静默污染。这正是 v2.5 使用独立目录 `results/raw_v25` 的原因。

## A.9 一键流水线

实验数据已就绪时，`main.py` 一条命令完成"单元测试 → 统计汇总与图表 → 生成 HTML 报告 → 自动打开报告"：

```powershell
python main.py                  # 默认流程
python main.py --with-experiments   # 先跑完四组 v2.4 实验
python main.py --skip-tests         # 跳过测试
python main.py --no-open            # 不自动打开浏览器
```

任一步失败即中止并打印错误。

## A.10 常见问题

| 现象 | 原因与处理 |
|---|---|
| `ModuleNotFoundError: No module named 'pymoo'` | 依赖未装好，或用了另一个 Python 解释器。确认 `python -c "import pymoo; print(pymoo.__version__)"` 能输出版本号 |
| 运行极慢 | 正常。全量 v2.4 约 3780 次运行、v2.5 约 2160 次，建议挂后台；或按 A.6 的提示按 (问题, $\tau$) 并行 |
| 想只跑自己关心的组合 | 所有实验脚本都接受 `--problems` / `--taus` / `--seeds` / `--limit-seeds`；用 `--help` 查看全部参数 |
| 中断后重新跑会不会重复计算 | 不会。运行器按 CSV 中已有的 (算法, 问题, $\tau$, 种子) 自动跳过已完成组合（**断点续跑**） |
| 论文表格里的数字怎么复核 | 运行 `python -m scripts.thesis_analysis`，它从原始逐代 CSV 从零重算，并输出可粘贴的 Markdown 表格 |
| HTML 报告里的文字和论文对不上 | 已知问题 K15：`build_report.py` 的叙述文字是 v1 硬编码，**表格数据是正确的**，叙述文字不可引用 |

## A.11 数据与代码的最小复现清单

要复现本文**第 5 章的全部数字**，只需：

1. `results/raw/`（v2.4 的 126 个 CSV）与 `results/raw_v25/`（v2.5 的 72 个 CSV）；
2. `scripts/thesis_analysis.py`（论文级汇总与统计）；
3. `python -m scripts.thesis_analysis`。

要**重新生成**这些数据，则需要上述 A.5、A.6 的全部命令。要**重新验证问题定义的正确性**，只需要 `python -m scripts.check_dmop_authoritative` 与 `python -m pytest -q`。


\newpage

# 附录 B 术语表

| 术语 | 英文 / 缩写 | 一句话解释 |
|---|---|---|
| 多目标优化问题 | MOP | 同时最小化多个相互冲突的目标，寻找"一组折中解" |
| Pareto 支配 | Pareto dominance | $a$ 支配 $b$：$a$ 在所有目标上不劣于 $b$，且至少在一个目标上严格更优 |
| Pareto 最优解集 | PS | 决策空间中不被任何其他解支配的解的集合 |
| Pareto 最优前沿 | PF | $PS$ 在目标空间中的像 |
| 动态多目标优化问题 | DMOP | 目标函数或约束随时间变化的 MOP；$PS(t)$ 与 $PF(t)$ 会移动 |
| Type I | — | 解集随时间变化、前沿不变（如 FDA1、权威 dMOP3） |
| Type II | — | 解集与前沿**都**随时间变化（如 FDA3、权威 dMOP2） |
| Type III | — | 解集不变、前沿随时间变化（如**权威 dMOP1**） |
| Type IV | — | 两者都不随时间变化（本项目的静态对照问题 DMOP3） |
| 渐变环境 | gradual environment | 单步变化幅度远小于决策空间直径、且方向具有短期一致性的环境序列（本文定义 2.1 给出可操作化定义） |
| NSGA-II | — | 最经典的精英多目标进化算法（快速非支配排序 + 拥挤度距离 + 精英保留）[1] |
| 快速非支配排序 | fast non-dominated sort | 把种群按支配关系分层，第一层即非支配解 |
| 拥挤度距离 | crowding distance | 同层内解的邻居间距，用于维持分布多样性 |
| 精英保留 | elitism | 父代与子代合并后择优进入下一代 |
| 模拟二进制交叉 | SBX | 一种实数编码交叉算子；本项目取交叉概率 0.9、分布指数 15 |
| 多项式变异 | PM | 一种实数编码变异算子；本项目取分布指数 20 |
| warm-start | 热启动 | 用上一阶段的解作为下一阶段的初始种群 |
| 环境响应 | respond | 算法在环境切换时生成新环境初始种群的接口 |
| 矩特征 | moment features | 分布的均值 $\mu$ 与标准差 $\sigma$（"中心在哪、散布多宽"） |
| 趋势外推 | trend extrapolation | 用最近两次中心差 $\mu_e-\mu_{prev}$ 估计下一中心 |
| 外推步长 / 趋势系数 | $\alpha$ | $\hat\mu=\mu_e+\alpha(\mu_e-\mu_{prev})$ 中的系数；$\alpha=1$ 为全步、$0.5$ 为半步 |
| AR(1) | 一阶自回归 | $\mu_k = a + b\mu_{k-1}$（对历史中心做最小二乘拟合） |
| GMM | 高斯混合模型 | 用多个高斯分量逼近任意形状分布 |
| PPS | 种群预测策略 | 用中心（及流形/方差）的外推预测新环境种群位置 [4] |
| KTMM | 知识迁移 + 混合模型 | 用高斯混合模型刻画动态最优解分布并跨环境迁移 [14] |
| 反向世代距离 | IGD | 从真实前沿每个参考点到算法解集的最近距离之均值，越小越好 [8] |
| 超体积 | HV | 解集与参考点围成的目标空间体积，越大越好 [9] |
| 稳态 MIGD | steady-state MIGD | **本文主指标**：每个环境段末代 IGD 的段间平均，再对种子平均（口径见 5.1.2 节） |
| 变化后 MIGD | post-change MIGD | 每个已变化环境段首代 IGD 的段间平均，衡量初始命中质量 |
| 全代 MHV | all-generation MHV | 一次运行全部代 HV 的平均（辅助口径，注意是 HV 而非 IGD） |
| 恢复代数 | recovery generations | 变化后 IGD 回落到该段末代水平所需代数；未达标则记为 $\tau_t$ |
| 迁移耗时 | migration time | 一次运行中 `respond` 的累计墙钟耗时 |
| $\tau_t$（tau） | — | 每次环境变化之间间隔的代数（变化周期） |
| $n_t$ | — | 全程不同环境状态总数；$t=\lfloor g/\tau_t\rfloor/n_t$ |
| $D$ | — | 决策变量维数 |
| $M$ | — | 目标数 |
| $N$ | — | 种群规模 |
| seed | 随机种子 | 换种子即得一次独立重复实验；本文用 0–29 共 30 个**对齐**种子 |
| Wilcoxon 符号秩检验 | paired signed-rank | 两算法成对比较的非参数检验 [10] |
| Holm 校正 | Holm step-down | 多重比较的 $p$ 值逐步校正，控制族错误率（FWER）[11] |
| Friedman 检验 | — | 多算法跨数据集整体排序的秩检验 [12] |
| Nemenyi 事后检验 | — | Friedman 显著后定位差异来源的临界差检验 [13] |
| 临界差 | CD | $CD=q_\alpha\sqrt{k(k+1)/6N}$；$k=6$、$N=30$、$\alpha=0.05$ 时 $CD=1.37668$ |
| 断点续跑 | resumable | 已完成的 (算法, 问题, $\tau$, 种子) 自动跳过，中断后可继续 |
| 消融实验 | ablation | 逐一关闭组件以量化各自的边际贡献 |
| 敏感性分析 | sensitivity analysis | 扫描超参数以检验结论对取值的稳健性 |
| 死参数 | dead parameter | 被赋予、但在实现路径上从不生效的参数（本文 K1 即此类缺陷） |


\newpage

# 附录 C 代码文件索引与论文章节映射

下表给出论文每一处论断所对应的代码位置。**左列为论文章节，右列为可实现核对的文件（含 v2.5 新增文件）。**

| 论文位置 | 代码文件 | 职责 |
|---|---|---|
| 3.1.3、4.3 | `algorithms/nsga2_core.py` | `make_nsga2` 算子工厂（SBX 0.9/η15、PM η20、去重）、`DynamicDriver` 环境分段驱动、`_algo_n_eval` 版本兼容 |
| 3.2、3.3、3.4、3.5、3.6 | `algorithms/ldmt_nsga2.py` | LDT-NSGA-II 主算法：归一化特征 → 趋势中心 → 高斯采样 → 50/20/30 重组；另含 v2.5 的 `LDTNSGA2Alpha1`（$\alpha=1.0$ 变体） |
| 3.2.2、3.2.3 | `algorithms/distribution_features.py` | `normalize` / `denormalize` / `extract_features`（含 `SIGMA_MIN=1e-3`）/ `non_dominated_set` |
| 2.3.2、3.6 | `algorithms/environment_manager.py` | `EnvironmentManager` 时间规则；`ChangeDetector.detect` 的变化幅度 $D_k$ |
| 3.5.3 | `algorithms/transfer_strategy.py` | 可插拔策略层（`TransferStrategy` / `MomentTransferStrategy` / `GMMTransferStrategy`），默认实验未启用 |
| 3.8 | `algorithms/dnsga2_response.py` | `NSGA2Static`（静态基线，`prev_X.copy()`）与 `DNSGA2Response`（保留前缀 + 随机移民） |
| 3.8 | `algorithms/ar_nsga2.py` | AR(1) 预测基线：逐维 OLS 自回归（历史 5）+ 20% 移民（K1 修复后） |
| 3.8 | `algorithms/gmm_transfer_nsga2.py` | GMM 分布迁移基线：sklearn `GaussianMixture`（≤3 分量，`reg_covar=1e-4`）+ **10%** 移民（K1 修复后） |
| 3.8 | `algorithms/pps_nsga2.py` | PPS 中心预测基线：$2\mu_{k-1}-\mu_{k-2}$ 全步外推 + 整体平移 + 20% 移民（K1 修复后） |
| 3.9.1 | `problems/fda1.py` / `fda3.py` / `fda4.py` | FDA 系列问题（与 PlatEMO 逐式核对一致） |
| 3.9.2 | `problems/dmop1a.py` / `dmop2a.py` / `dmop3a.py` | **权威** dMOP1（Type III）/ dMOP2（Type II）/ dMOP3（Type I），v2.5 新增 |
| 3.9.3 | `problems/dmop1.py` / `dmop2.py` / `dmop3.py` | 三个**项目变体**（DMOP3 为静态对照），如实标注、不得称标准 dMOP |
| 3.9.4 | `problems/true_fronts.py` | 参考前沿生成（9 个问题的解析前沿；按 (问题, $M$, 点数, $t$) 缓存；FDA4 去重后 507 点） |
| 2.3.2、3.6 | `problems/base.py` | `DynamicProblem`：$t=\lfloor g/\tau_t\rfloor/n_t$ 时间规则 + 共享固定探针（10 点、seed 2026、真实边界） |
| 2.5 | `evaluation/igd.py` | IGD 计算（scipy `cdist` 欧氏距离） |
| 2.5 | `evaluation/hv.py` | HV 计算（pymoo 精确 HV）+ 逐问题固定参考点 + K5 修复后的容斥回退 |
| 2.6 | `evaluation/statistics.py` | `holm_adjust` 逐步下降校正 + `wilcoxon_holm` 配对/独立检验 |
| 2.6 | `evaluation/friedman_test.py` | Friedman 卡方统计量（`rankdata` 平均秩）+ Nemenyi 双档 $q$ 表与临界差 |
| 4.3、4.4 | `experiments/common.py` | `RunConfig`、算法/问题注册表、`ExperimentRunner`（断点续跑 + 运行级元数据） |
| 4.3 | `experiments/run_comparison.py` / `run_supplementary.py` / `run_ablation.py` / `run_sensitivity.py` | v2.4 四组实验脚本 |
| 4.3、5.7、5.9.2 | `experiments/run_v25_alpha.py` | v2.5 两阶段实验驱动（$\alpha=1.0$ 对照 + 权威 dMOP 补充），含"拒绝写入冻结目录"保护 |
| 5.1.2、5.1.3 | `scripts/thesis_analysis.py` | **论文级**汇总与统计：从原始逐代 CSV 从零重算全部表格（独立于 `analyze.py`） |
| 3.9.2 | `scripts/check_dmop_authoritative.py` | 权威 dMOP 自检：解析 PS 落在解析 PF 上的偏差、问题类型、$r$ 切换策略 |
| 4.5 | `scripts/make_thesis_figures.py` | 论文插图生成（`results/figures_thesis/`） |
| 4.5 | `analyze.py` | v2.4 冻结汇总与图表（`results/summary/`、`results/figures/`） |
| 4.4 | `main.py` | 一键流水线（pytest → analyze → build_report → 打开报告） |
| 4.5 | `build_report.py` / `serve_report.py` | HTML 报告生成与本地预览（报告叙述文字见 K15） |
| 4.6 | `tests/`（9 文件 63 测试） | 见 4.6 节表格 |
| 附录 D | 本附录 | 核心模块**完整可运行源码**（由构建脚本从仓库原文件逐字拷贝，非伪代码、非手抄） |

## C.1 v2.5 相对于 v2.4 的代码改动清单

| 文件 | 改动 | 是否改变已有实验行为 |
|---|---|---|
| `algorithms/ldmt_nsga2.py` | 新增类 `LDTNSGA2Alpha1`（`name="ldt_alpha1"`，`trend_alpha=1.0`，其余继承） | **否**（新类，不改默认值） |
| `algorithms/__init__.py` | 导出 `LDTNSGA2Alpha1` | 否 |
| `experiments/common.py` | 注册 `ldt_alpha1`；注册 `DMOP1A/2A/3A` | 否（不在 `CANONICAL_ALGOS` 中） |
| `problems/dmop1a.py` / `dmop2a.py` / `dmop3a.py` | 新增权威 dMOP 问题类 | 否（新文件） |
| `problems/__init__.py` | 导出新增问题类 | 否 |
| `problems/true_fronts.py` | 新增权威 dMOP 前沿函数；**修正 `get_true_pf` 的默认目标数判定**（原实现按问题名硬编码，改为 `M=3` 仅对 FDA4/FDA5） | 否（对原 6 问题结果不变） |
| `evaluation/hv.py` | 新增 `DMOP1A/2A/3A` 的参考点 | 否（新问题） |
| `algorithms/nsga2_core.py` | `_row_metrics` 的时变前沿判定加入 `DMOP2A`（并保留 `FDA3`、`DMOP2`） | 否（`DMOP1A` 的时变前沿经修正后也需按 $t$ 生成，见下表说明） |
| `experiments/run_v25_alpha.py`、`scripts/`（3 个） | 新增 | 否 |

> **施工记录（如实）**：权威 dMOP 的实现曾经历一次**更正**。首版把 dMOP1 误按 Type I（含 $G(t)$ 的 $g$ + $\sqrt{\cdot}$ 前沿）实现、把 dMOP3 误用项目变体的欺骗性前沿与 $H(t)$ 解集位置。经文献复核（见 `_audit/dmop_definitions.md`）后重写为正确形式，并把首版生成的 12 个 CSV 移出统计目录、归档到 `results/raw_v25_badauthoritative/` 后**重新运行**。此外，首版还发现 `_row_metrics` 未把 `DMOP1A` 列入"按 $t$ 取前沿"的分支——由于 dMOP1 是 Type III（前沿随 $t$ 变化），若沿用 $t=0$ 的前沿会系统性高估 IGD；已一并修正。

## C.2 复现所需的完整命令清单

见附录 A 的第 4–7 步。核心是三条命令：

```powershell
python -m experiments.run_v25_alpha --problems FDA1,FDA3,FDA4 --taus 10,25,50 `
    --seeds 0-29 --out results/raw_v25 --phase comparison
python -m experiments.run_v25_alpha --problems DMOP1A,DMOP2A,DMOP3A --taus 25 `
    --seeds 0-29 --out results/raw_v25 --phase supplementary
python -m scripts.thesis_analysis
```


\newpage

# 附录 D 核心模块完整源代码

**本附录的全部代码由构建脚本从仓库原文件逐字拷贝**（构建时通过 `<!--INCLUDE:路径-->` 指令展开），因此它**不是伪代码、不是手抄稿**，而是产生第 5 章全部实验结果的同一份代码。读者可直接复制这些文件到相应路径运行。

代码为 Python 3.9（实测 3.9.13）+ pymoo 0.6.1.5。全部文件位于项目根目录下的 `algorithms/`、`evaluation/` 与 `problems/` 中，运行方式见附录 A。

## D.1 LDT-NSGA-II 主算法

**文件**：`algorithms/ldmt_nsga2.py`

对应论文章节 3.1–3.6。核心流程为：非支配子集矩特征 → 趋势中心预测 → 自适应缩放（可选）→ 多源融合（可选）→ 高斯采样迁移 → 50/20/30 三路重组 → 单次反归一化。文件末尾另含 v2.5 新增的 `LDTNSGA2Alpha1`（$\alpha=1.0$ 变体，即 5.7 节实验所用的算法）。

<!--INCLUDE:algorithms/ldmt_nsga2.py-->

## D.2 NSGA-II 内核与环境分段驱动器

**文件**：`algorithms/nsga2_core.py`

对应论文章节 3.1.3、4.3。包含算子工厂 `make_nsga2`（SBX 0.9/$\eta=15$、PM $\eta=20$、`eliminate_duplicates=True`）、pymoo 评估计数兼容封装 `_algo_n_eval`、参考前沿选择（含时变前沿判定）以及 `DynamicDriver` 的逐代运行循环。

<!--INCLUDE:algorithms/nsga2_core.py-->

## D.3 分布特征与归一化

**文件**：`algorithms/distribution_features.py`

对应论文章节 3.2。包含 `normalize` / `denormalize` / `extract_features`（逐维均值与总体标准差，`ddof=0`）、标准差下限常量 `SIGMA_MIN = 1e-3`、非支配集判定 `non_dominated_set`（带 $10^{-8}$ 容差）。

<!--INCLUDE:algorithms/distribution_features.py-->

## D.4 环境时间管理与变化幅度 $D_k$

**文件**：`algorithms/environment_manager.py`

对应论文章节 2.3.2、3.6。包含时间规则 `EnvironmentManager.time_at` 与变化幅度计算 `ChangeDetector.detect`（按目标维归一化后的平均相对偏移）。

<!--INCLUDE:algorithms/environment_manager.py-->

## D.5 性能指标：IGD

**文件**：`evaluation/igd.py`

对应论文章节 2.5。使用 SciPy 的 `cdist` 计算欧氏距离矩阵，取每行最小值后求均值。

<!--INCLUDE:evaluation/igd.py-->

## D.6 变化响应基线之一：AR(1) 预测

**文件**：`algorithms/ar_nsga2.py`

对应论文章节 3.8 与已知问题 K1。包含逐维 OLS 一阶自回归拟合、斜率不稳定时的回退逻辑，以及"**先预留移民名额再拼接**"的 K1 修复实现。

<!--INCLUDE:algorithms/ar_nsga2.py-->

## D.7 变化响应基线之二：PPS 中心预测

**文件**：`algorithms/pps_nsga2.py`

对应论文章节 3.8 与 5.7 节。中心历史容量 3，全步一阶差分外推 $\hat\mu = \mathrm{clip}(2\mu_{k-1}-\mu_{k-2}, 0, 1)$，种群整体平移 + 20% 随机移民。

<!--INCLUDE:algorithms/pps_nsga2.py-->

## D.8 变化响应基线之三：GMM 全分布迁移

**文件**：`algorithms/gmm_transfer_nsga2.py`

对应论文章节 3.8 与 5.5.3 节。使用 scikit-learn 的 `GaussianMixture`（分量数 $\min(3,\max(1,N/10))$、`reg_covar=1e-4`、协方差类型按 diag→full→spherical 依次尝试），**注意默认 `immigrant_ratio = 0.1`（10%）**，与另两个预测基线的 20% 不同（见 3.8 节的诚实声明）。

<!--INCLUDE:algorithms/gmm_transfer_nsga2.py-->

## D.9 权威 dMOP1（Type III）源码

**文件**：`problems/dmop1a.py`

对应论文章节 3.9.2 与 5.9.2 节。**这是本文写作期间经文献复核更正过的实现**：权威 dMOP1 是 Type III（解集静止、前沿随 $H(t)$ 旋转），$g = 1 + 9\sum_{i\ge2}x_i^2$（不含 $G(t)$），前沿为 $1-u^{H(t)}$ 而非 $\sqrt{\cdot}$ 型。文件头注释记录了完整的出处与更正理由。

<!--INCLUDE:problems/dmop1a.py-->

## D.10 权威 dMOP2（Type II）源码

**文件**：`problems/dmop2a.py`

对应论文章节 3.9.2。$g = 1 + 9\sum_{i\ge2}(x_i-G)^2$、$H(t)=0.75\sin(\pi t/2)+1.25$、$f_2 = g(1-(f_1/g)^{H(t)})$，解集与前沿同时变化。

<!--INCLUDE:problems/dmop2a.py-->

## D.11 权威 dMOP3（Type I）源码

**文件**：`problems/dmop3a.py`

对应论文章节 3.9.2。$f_1 = x_r$（$r$ 为按环境轮换的"位置相关"变量下标）、$g = 1 + 9\sum_{i\ne r}(x_i-G(t))^2$、前沿 $1-\sqrt{f_1}$ 静止而解集随 $t$ 移动。文件头注释记录了 $G(t)$ 采用 $\lvert\sin(\pi t/2)\rvert$ 变体的理由与文献分歧。

<!--INCLUDE:problems/dmop3a.py-->

## D.12 解析参考前沿

**文件**：`problems/true_fronts.py`

对应论文章节 3.9.4。含 9 个问题（FDA1/3/4、项目变体 DMOP1/2/3、权威 DMOP1A/2A/3A）的解析 Pareto 前沿、FDA4 的确定性角度网格球面采样（$23\times23=529$ 去重后 507 点），以及按 (问题, $M$, 点数, $t$) 的缓存。

<!--INCLUDE:problems/true_fronts.py-->

## D.13 权威 dMOP 自身一致性自检脚本

**文件**：`scripts/check_dmop_authoritative.py`

对应论文章节 3.9.2 的自检证据。把决策向量放在**解析导出的 Pareto 集**上，检查其目标值是否落在**解析导出的 Pareto 前沿**上（实测最大偏差 0.000e+00），并验证三个问题的动态类型（Type I/II/III）与 $r$ 的切换策略。

<!--INCLUDE:scripts/check_dmop_authoritative.py-->

## D.14 论文级汇总与统计分析脚本

**文件**：`scripts/thesis_analysis.py`

对应论文章节 5.1.2、5.1.3 与 4.5。直接从 `results/raw` 与 `results/raw_v25` 的逐代原始 CSV 从零重算论文全部表格：稳态 MIGD/MHV、变化后 MIGD、恢复代数、迁移耗时、配对 Wilcoxon + Holm、Friedman + Nemenyi、以及 $\alpha=1.0$ 与 $\alpha=0.5$ 的配对对照。**该脚本不读取也不写入 `results/summary/`**，以保证与 `analyze.py` 相互独立、可交叉验证。

<!--INCLUDE:scripts/thesis_analysis.py-->


\newpage

# 附录 E 实验数据文件清单与完整性核验

## E.1 数据目录结构

| 目录 | 内容 | 是否参与论文统计 |
|---|---|---|
| `results/raw/` | v2.4 四组实验的逐代原始 CSV（126 个文件）与 `run_meta.jsonl` | ✅ 参与（第 5.2–5.6、5.8–5.9.1 节） |
| `results/raw_v25/` | v2.5 两阶段实验的逐代原始 CSV（72 个文件） | ✅ 参与（第 5.7、5.9.2 节） |
| `results/raw_v25_badauthoritative/` | 首版**错误** dMOP1/dMOP3 定义下的 12 个 CSV | ❌ **不参与**，仅作施工留痕 |
| `results/raw_legacy_v21/` | v2.1 时期的旧数据 | ❌ 不参与 |
| `results/raw_legacy_v22_k1/` | K1 修复前（三基线移民被截断）的数据 | ❌ 不参与，归档备查 |
| `results/summary/` | v2.4 冻结汇总（由 `analyze.py` 生成） | 交叉验证用 |
| `results/summary_v25/` | 论文级汇总（由 `scripts/thesis_analysis.py` 生成） | ✅ 论文表格的**直接来源** |
| `results/figures/` | 汇总图表（`analyze.py`） | 论文未直接使用 |
| `results/figures_thesis/` | 论文插图（`scripts/make_thesis_figures.py`） | ✅ 正文插图 |

## E.2 逐代 CSV 的字段定义（18 项）

```text
algorithm, problem, seed, population_size, total_generations, n_t, tau_t,
environment_index, generation, IGD, HV, generation_runtime_s,
migration_runtime_s, n_evaluations, code_version, lambda_k,
boundary_clip_ratio, n_nd
```

| 字段 | 含义 |
|---|---|
| `algorithm` | 算法规范标识（`nsga2_static` … `ldt_nsga2` / `ldt_alpha1`） |
| `problem` | 问题标识（`FDA1`…`DMOP3` / `DMOP1A`…`DMOP3A`） |
| `seed` | 随机种子（0–29） |
| `population_size` / `total_generations` / `n_t` / `tau_t` | 协议参数 $N$ / $n_{gen}$ / $n_t$ / $\tau_t$ |
| `environment_index` | 环境段序号 $k$（0 基） |
| `generation` | 全局代数 $g$（0 基） |
| `IGD` / `HV` | 当前代的指标值（出现 NaN/Inf 即判该次运行失败） |
| `generation_runtime_s` | 该代 wall-clock 耗时 |
| `migration_runtime_s` | **仅在该段第一代（$k>0$）写入**实测迁移耗时，其余代记 0 |
| `n_evaluations` | pymoo 评估计数器的**环境内累计值**（$\tau=10/25/50$ 段末分别为 1000/2500/5000） |
| `code_version` | 代码版本戳（见 E.5 的口径说明） |
| `lambda_k` | LDT 内联路径恒为 1.0 |
| `boundary_clip_ratio` | 贴界个体占比（诊断字段，K9 修复后为真实值） |
| `n_nd` | 该次响应所用的非支配子集规模（LDT 路径） |

## E.3 运行级元数据

`results/raw/run_meta.jsonl` 与 `results/raw_v25/run_meta.jsonl` 逐运行一行，字段包括：算法规范名与 CLI 名、问题、$\tau_t$、$n_t$、$N$、$n_{gen}$、种子、`n_eval_total` 与**逐环境评估数** `n_eval_per_env`、耗时、**git 提交号（前 12 位）**、`code_version`、`status`，以及 Python/numpy/pymoo/sklearn 版本与平台信息。失败记录单独写入 `failed_runs.jsonl`。

**v2.4 全部 3780 次运行中 `results/raw/failed_runs.jsonl` 不存在**，即 0 失败；v2.5 两阶段亦为 0 失败。

## E.4 完整性核验结果

由 `scripts/thesis_analysis.py` 输出的 `results/summary_v25/integrity.json` 与 `integrity_runs_by_prefix.csv` 给出：

| 核验项 | 结果 |
|---|---|
| v2.4 主比较（comparison） | 54 个 CSV；3 问题 × 3 档 $\tau$ × 6 算法，每格 **30 个种子，无缺种**；每运行 300 代 |
| v2.4 补充（supplementary） | 18 个 CSV；3 问题 × 6 算法，每格 30 种子 |
| v2.4 消融（ablation） | 12 个 CSV；3 问题 × 4 变体，每格 30 种子 |
| v2.4 敏感性（sensitivity） | 42 个 CSV；$N$ 三档 270 + $n_t$ 三档 270 + $\rho$ 三档 270 + $\alpha$ 五档 450 = **1260** |
| v2.5 阶段一（comparison） | 54 个 CSV；3 问题 × 3 档 $\tau$ × 6 算法，每格 30 种子 |
| v2.5 阶段二（supplementary） | 18 个 CSV；3 问题（DMOP1A/2A/3A）× 6 算法，每格 30 种子 |
| IGD/HV 有限性 | 全部字段有限，**无 NaN/Inf** |
| 逐代行数 | 300（或敏感性中 $n_t$ 变化时的 600/3600），**无截断** |

> **两处措辞的精确化（见 K21）**：论文所称"3780 次运行"指**运行记录数**。若按不重复的 (algorithm, problem, $\tau$, seed) 组合统计，v2.4 的 distinct 组合数为 **2430**；若进一步区分 $N$ 与 $n_t$，配置数为 **2790**。差额来自默认 LDT 配置在 comparison、ablation、$N{=}100$、$n_t{=}20$、$\rho{=}0.5$、$\alpha{=}0.5$ 六处被重复执行（这正是敏感性/消融实验的设计意图——它们本就需要与默认配置对照）。论文正文统一采用"运行记录数"这一表述，并在 5.1.1 节说明差异。

## E.5 版本戳的口径说明（K20）

| 对象 | 值 |
|---|---|
| 磁盘 v2.4 原始 CSV 中的 `code_version` | **`v2.2`**（126 个文件、1,134,000 行全部如此） |
| 磁盘 v2.4 `run_meta.jsonl` 中的 `code_version` | **`v2.2`**（3240 条记录全部如此） |
| `DynamicDriver.code_version` 当前代码值 | `v2.4` |
| 数据生成提交（`run_meta` 记录的 git commit） | `50adf19c5c5a` |
| 论文分析所用代码提交 | 见附录 E.6 |

原因是 `run_comparison.py` 至今仍显式传入 `code_version="v2.2"`。该字段纯粹是标记，**不影响任何数值结果**；论文如实说明并在引用数据时同时给出数据生成提交与代码分析提交。v2.5 的数据统一标为 `v2.5`。

## E.6 复核与版本记录

| 事项 | 记录 |
|---|---|
| v2.4 数据生成提交 | `50adf19c5c5a` |
| 论文写作期间的仓库 HEAD | `d4ce930ad86a7ff247488fdfa05a04370a957949`（记为 d4ce930） |
| 独立数值复核 | 一个不读取任何既有汇总、仅从逐代 CSV 从零重算的复核流程，逐项比对了 **269** 个数值断言：**265 项一致、4 项不一致、0 项无法判定**。4 处不一致即 K18–K21，已全部修正 |
| 双实现互检 | `analyze.py` 与 `scripts/thesis_analysis.py` 对 v2.4 主比较的 54 个单元格独立计算，四位小数完全一致 |
| 源码审计 | 独立审计复核了算法常量、基线行为、实验协议、指标口径、统计实现、测试数量与 K1–K12 的处置状态，发现 1 处实质错误（GMM 移民比例为 10% 而非 20%，已在 3.8 节更正）与 5 处部分不准确（已在 3.3–3.7 节按审计结果精确化） |
| 参考文献复核 | 逐条核对 Crossref / 出版方页面；补齐了原稿标注"未能核实"的两条关键文献（PPS 的卷期页码与 KTMM 的完整著录），更正了作者姓名（Zhang Q-F → Zhang Q）；并确认 PlatEMO master（commit `d25e65d`）不含 DMOP1/2/3.m |
| 权威 dMOP 定义复核 | 原论文受付费墙限制无法获取，改用四份独立复现文献交叉确认；据此**更正了 dMOP1 的类型（Type I → Type III）与 dMOP3 的公式**，并重新运行受影响的 1080 次实验（360 次为失败首版作废、540 次为更正后重跑，见 E.1 的归档目录） |


\newpage

# 致谢

本论文的完成得益于许多人的帮助，在此谨致谢意。

首先感谢我的指导教师。从选题方向的确定、研究方案的推敲，到实验设计的规范化和论文写作的反复修改，导师都给予了悉心指导。导师在动态多目标优化与知识迁移方向的持续研究为本课题提供了重要的思想来源，其强调"结论必须有可核查的证据支撑、不利结果同样要如实报告"的治学态度，直接塑造了本文的写作方式——包括本文如实保留"默认参数 α=0.5 被实验证伪"这一对原设计不利的发现，以及以 K1–K26 清单完整披露实现缺陷的做法。

感谢课题组与实验室的各位老师、同学在讨论中提出的意见，特别是在实验口径与统计检验规范方面的建议。感谢开源社区：PlatEMO 项目为本课题的 FDA 测试问题提供了逐式核对的权威基准，pymoo 项目提供了稳定可靠的 NSGA-II 与遗传算子实现，SciPy 与 scikit-learn 提供了统计检验与高斯混合模型的实现。没有这些高质量的开源工作，本课题的实验系统无法在合理时间内建立。

本文在写作过程中使用了若干自动化工具辅助完成重复性工作，包括：对全部源码与数据的逐文件梳理、对论文表格数值的独立重算复核（269 项断言）、对参考文献著录信息与 benchmark 定义的联网核查。这些工具的输出均经过人工判断后采用；凡未能核实之处，本文均显式标注"未能核实"，不做推测填补。

最后感谢家人多年来的理解与支持。谨以此文，献给所有在动态变化中仍坚持寻找更好解的人。

---

**关于本文档的最终说明**：本文全部数值均来自真实运行数据（v2.4 3780 次运行记录 + v2.5 2250 次运行记录 = **6030 次**）或代码常量；全部关键技术节点均标注文献出处（GB/T 7714 顺序编码）或项目内代码位置（文件路径）。作者无法核实之处已在正文中明确声明。若更换种子集合或代码版本重跑，应以新的 `results/summary_v25/analysis_report.md` 汇总输出为准，并同步更新第 5 章与摘要。
