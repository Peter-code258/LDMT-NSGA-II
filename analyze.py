"""Analyze all experiment CSVs: summary tables, Wilcoxon+Holm tests, figures,
and a narrative markdown report. Run after experiments complete.

    python analyze.py
"""
import os
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

warnings.filterwarnings("ignore")
RAW = "results/raw"
SUMM = "results/summary"
FIG = "results/figures"
for d in (SUMM, FIG):
    os.makedirs(d, exist_ok=True)

ALGO_ORDER = ["DNSGA2", "AR_NSGA2", "GMM_Transfer", "LDMT_NSGA2"]
ALGO_LABEL = {
    "DNSGA2": "DNSGA-II",
    "AR_NSGA2": "AR-NSGA-II",
    "GMM_Transfer": "GMM-Transfer",
    "LDMT_NSGA2": "LDMT-NSGA-II",
}
COLORS = {"DNSGA2": "#7f7f7f", "AR_NSGA2": "#1f77b4",
          "GMM_Transfer": "#ff7f0e", "LDMT_NSGA2": "#d62728"}


def load(pattern):
    files = sorted(glob.glob(os.path.join(RAW, pattern)))
    if not files:
        return None
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def per_run_metrics(df):
    """One row per run: MIGD, MHV, post-change IGD, times."""
    g = df.copy()
    g["gen_in_env"] = g["generation"] % g["tau_t"]
    run_keys = ["algorithm", "problem", "tau_t", "population_size", "n_t", "seed"]
    # last generation of each env (steady-state) and first gen after change,
    # grouped per RUN (algorithm, problem, tau_t, seed) so different runs
    # are not mixed together.
    last = g.groupby(run_keys + ["environment_index"]).apply(
        lambda x: x.loc[x["generation"].idxmax()])
    first = g.groupby(run_keys + ["environment_index"]).apply(
        lambda x: x.loc[x["generation"].idxmin()])
    first = first.reset_index(drop=True)
    last = last.reset_index(drop=True)
    rows = []
    for keys, sub in last.groupby(run_keys):
        algo, prob, tau, N, n_t, seed = keys
        first_m = first[(first["algorithm"] == algo)
                        & (first["problem"] == prob)
                        & (first["tau_t"] == tau)
                        & (first["population_size"] == N)
                        & (first["n_t"] == n_t)
                        & (first["seed"] == seed)]
        first_changed = first_m[first_m["environment_index"] > 0]
        run_full = g[(g["algorithm"] == algo) & (g["problem"] == prob)
                     & (g["tau_t"] == tau) & (g["population_size"] == N)
                     & (g["n_t"] == n_t) & (g["seed"] == seed)]
        rows.append({
            "algorithm": algo,
            "problem": prob,
            "tau_t": int(tau),
            "n_t": int(sub["n_t"].iloc[0]),
            "N": int(sub["population_size"].iloc[0]),
            "seed": int(seed),
            "MIGD": sub["IGD"].mean(),
            "MHV": sub["HV"].mean(),
            "post_change_MIGD": first_changed["IGD"].mean()
                if len(first_changed) else np.nan,
            "mean_gen_time_s": run_full["generation_runtime_s"].mean(),
            "total_mig_time_s": run_full["migration_runtime_s"].sum(),
        })
    return pd.DataFrame(rows)


def summarize(runs):
    """mean/std/median/IQR per group."""
    out = runs.groupby(["problem", "tau_t", "algorithm"]).agg(
        MIGD_mean=("MIGD", "mean"), MIGD_std=("MIGD", "std"),
        MIGD_med=("MIGD", "median"),
        MHV_mean=("MHV", "mean"), MHV_std=("MHV", "std"),
        post_MIGD_mean=("post_change_MIGD", "mean"),
        gen_time_mean=("mean_gen_time_s", "mean"),
        mig_time_mean=("total_mig_time_s", "mean"),
        n_runs=("seed", "nunique"),
    ).reset_index()
    return out


def wilcoxon_holm(runs, ref="LDMT_NSGA2"):
    """Pairwise Wilcoxon (LDMT vs each other) on MIGD per (problem,tau)."""
    rows = []
    for (prob, tau), grp in runs.groupby(["problem", "tau_t"]):
        ref_vals = grp[grp["algorithm"] == ref].set_index("seed")["MIGD"]
        comps = [a for a in ALGO_ORDER if a != ref]
        pvals = []
        for a in comps:
            o = grp[grp["algorithm"] == a].set_index("seed")["MIGD"]
            common = ref_vals.index.intersection(o.index)
            if len(common) < 5:
                pvals.append((a, np.nan))
                continue
            try:
                p = stats.wilcoxon(ref_vals.loc[common], o.loc[common]).pvalue
            except Exception:
                p = np.nan
            pvals.append((a, p))
        # Holm correction
        valid = [(a, p) for a, p in pvals if not np.isnan(p)]
        ps = sorted(range(len(valid)), key=lambda i: valid[i][1])
        m = len(valid)
        for rank, i in enumerate(ps):
            a, p = valid[i]
            adj = min(1.0, p * (m - rank))
            rows.append({"problem": prob, "tau_t": tau, "comparison": f"{ref} vs {a}",
                         "p_value": p, "p_holm": adj,
                         "sig": "***" if adj < 0.001 else ("**" if adj < 0.01 else ("*" if adj < 0.05 else "ns"))})
    return pd.DataFrame(rows)


def fmt(x, dec=4):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "-"
    return f"{x:.{dec}f}"


# ---------- load ----------
comp = load("comparison_*.csv")
abl = load("ablation_*.csv")
sens = load("sensitivity_*.csv")
supp = load("supplementary_*.csv")

report_lines = ["# LDMT-NSGA-II 实验分析报告\n",
                f"生成于实验完成后。原始数据见 results/raw，汇总见 results/summary，图见 results/figures。\n"]

# ---------- main comparison ----------
if comp is not None:
    runs = per_run_metrics(comp)
    runs.to_csv(os.path.join(SUMM, "comparison_runs.csv"), index=False)
    summ = summarize(runs)
    summ.to_csv(os.path.join(SUMM, "summary_comparison.csv"), index=False)
    tests = wilcoxon_holm(runs)
    tests.to_csv(os.path.join(SUMM, "statistical_tests.csv"), index=False)

    report_lines.append("## 1. 主对比实验 (FDA1/3/4 × τ_t∈{10,25,50} × 4算法 × 30次)\n")
    report_lines.append("MIGD(越小越好)、MHV(越大越好)、变化后MIGD(响应速度)、单代时间、迁移耗时。\n\n")
    report_lines.append("| problem | τ_t | algorithm | MIGD | MHV | postΔMIGD | gen_t(s) | mig_t(s) |")
    report_lines.append("\n|---|---|---|---|---|---|---|---|")
    for _, r in summ.iterrows():
        report_lines.append(
            f"\n| {r['problem']} | {int(r['tau_t'])} | {ALGO_LABEL.get(r['algorithm'], r['algorithm'])} | "
            f"{r['MIGD_mean']:.4f}±{r['MIGD_std']:.4f} | {r['MHV_mean']:.4f}±{r['MHV_std']:.4f} | "
            f"{fmt(r['post_MIGD_mean'])} | {fmt(r['gen_time_mean'],5)} | {fmt(r['mig_time_mean'],5)} |")

    # best baseline gap (relative MIGD vs best non-LDMT)
    report_lines.append("\n\n### 1.1 LDMT 相对最佳基线的 MIGD 差距\n")
    report_lines.append("| problem | τ_t | best_baseline | LDMT | rel_gap% |")
    report_lines.append("\n|---|---|---|---|---|")
    for (prob, tau), grp in summ.groupby(["problem", "tau_t"]):
        nonl = grp[grp["algorithm"] != "LDMT_NSGA2"]
        if nonl.empty:
            continue
        bb = nonl.loc[nonl["MIGD_mean"].idxmin()]
        ld = grp[grp["algorithm"] == "LDMT_NSGA2"]
        if ld.empty:
            continue
        ld = ld.iloc[0]
        gap = (ld["MIGD_mean"] - bb["MIGD_mean"]) / bb["MIGD_mean"] * 100
        report_lines.append(
            f"\n| {prob} | {int(tau)} | {ALGO_LABEL.get(bb['algorithm'], bb['algorithm'])} ({bb['MIGD_mean']:.4f}) | "
            f"{ld['MIGD_mean']:.4f} | {gap:+.2f}% |")

    report_lines.append("\n\n### 1.2 统计检验 (Wilcoxon 秩和 + Holm 校正, LDMT vs 基线, MIGD)\n")
    report_lines.append("| problem | τ_t | comparison | p | p_holm | sig |")
    report_lines.append("\n|---|---|---|---|---|---|")
    for _, t in tests.iterrows():
        report_lines.append(
            f"\n| {t['problem']} | {int(t['tau_t'])} | {t['comparison']} | {fmt(t['p_value'])} | {fmt(t['p_holm'])} | {t['sig']} |")

    # ---- figures ----
    # (1) boxplots of MIGD per problem per tau
    for prob in runs["problem"].unique():
        sub = runs[runs["problem"] == prob]
        fig, axes = plt.subplots(1, len(sub["tau_t"].unique()),
                                 figsize=(4 * len(sub["tau_t"].unique()), 4),
                                 squeeze=False)
        for j, tau in enumerate(sorted(sub["tau_t"].unique())):
            ax = axes[0][j]
            data, labels, cols = [], [], []
            for a in ALGO_ORDER:
                d = sub[(sub["tau_t"] == tau) & (sub["algorithm"] == a)]["MIGD"].values
                if len(d):
                    data.append(d); labels.append(ALGO_LABEL[a]); cols.append(COLORS[a])
            bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False)
            for patch, c in zip(bp["boxes"], cols):
                patch.set_facecolor(c); patch.set_alpha(0.7)
            ax.set_title(f"{prob} τ_t={tau}")
            ax.set_ylabel("MIGD" if j == 0 else "")
            ax.tick_params(axis="x", labelrotation=30)
        fig.suptitle(f"{prob}: MIGD 分布 (30 runs)")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, f"boxplot_MIGD_{prob}.png"), dpi=120)
        plt.close(fig)

    # (2) recovery curve: mean IGD vs gen-in-env, per problem at tau=25
    for prob in runs["problem"].unique():
        sub = comp[(comp["problem"] == prob) & (comp["tau_t"] == 25)]
        if sub.empty:
            continue
        sub = sub.copy(); sub["gen_in_env"] = sub["generation"] % sub["tau_t"]
        # only envs that actually changed (env_index>0)
        sub = sub[sub["environment_index"] > 0]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for a in ALGO_ORDER:
            d = sub[sub["algorithm"] == a].groupby("gen_in_env")["IGD"].mean()
            if len(d):
                ax.plot(d.index, d.values, marker="o", ms=3,
                        label=ALGO_LABEL[a], color=COLORS[a], lw=1.8)
        ax.set_xlabel("环境变化后代数 (0=刚变化)")
        ax.set_ylabel("IGD (均值, 30 runs)")
        ax.set_title(f"{prob} (τ_t=25): 环境变化后恢复曲线")
        ax.legend(); ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, f"recovery_{prob}.png"), dpi=120)
        plt.close(fig)

    # (3) migration-time bar chart per problem (tau=25)
    sub = runs[runs["tau_t"] == 25]
    if not sub.empty:
        probs = sorted(sub["problem"].unique())
        x = np.arange(len(probs)); w = 0.2
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for k, a in enumerate(ALGO_ORDER):
            vals = [sub[(sub["problem"] == p) & (sub["algorithm"] == a)]["total_mig_time_s"].mean() for p in probs]
            ax.bar(x + (k - 1.5) * w, vals, w, label=ALGO_LABEL[a], color=COLORS[a], alpha=0.85)
        ax.set_xticks(x); ax.set_xticklabels(probs)
        ax.set_ylabel("每次运行总迁移耗时 (s)")
        ax.set_title("迁移模块耗时对比 (τ_t=25, 30 runs 均值)")
        ax.legend(); ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "migration_time.png"), dpi=120)
        plt.close(fig)

    # (4) MIGD vs tau_t (trend with update frequency)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for j, prob in enumerate(sorted(runs["problem"].unique())):
        ax = axes[j]
        sub = runs[runs["problem"] == prob]
        for a in ALGO_ORDER:
            d = sub[sub["algorithm"] == a].groupby("tau_t")["MIGD"].mean()
            if len(d):
                ax.plot(d.index, d.values, marker="o", label=ALGO_LABEL[a], color=COLORS[a])
        ax.set_title(prob); ax.set_xlabel("τ_t"); ax.set_ylabel("MIGD")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.suptitle("MIGD 随环境更新频率 τ_t 变化")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "MIGD_vs_tau.png"), dpi=120)
    plt.close(fig)
else:
    report_lines.append("## 1. 主对比实验: 未找到 comparison_*.csv\n")

# ---------- ablation ----------
if abl is not None:
    runs_a = per_run_metrics(abl)
    runs_a.to_csv(os.path.join(SUMM, "ablation_runs.csv"), index=False)
    sa = summarize(runs_a)
    sa.to_csv(os.path.join(SUMM, "ablation_summary.csv"), index=False)
    report_lines.append("\n## 2. 消融实验 (FDA1/3/4, τ_t=25, 30次)\n")
    report_lines.append("变体: full(完整LDMT) / no_feature(改用GMM特征) / no_transfer(去迁移映射)。\n\n")
    report_lines.append("| problem | variant | MIGD | MHV | mig_t(s) |")
    report_lines.append("\n|---|---|---|---|---|")
    for _, r in sa.iterrows():
        report_lines.append(
            f"\n| {r['problem']} | {r['algorithm']} | {r['MIGD_mean']:.4f}±{r['MIGD_std']:.4f} | "
            f"{r['MHV_mean']:.4f} | {fmt(r['mig_time_mean'],5)} |")
    # ablation bar chart
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for j, prob in enumerate(sorted(runs_a["problem"].unique())):
        ax = axes[j]; sub = runs_a[runs_a["problem"] == prob]
        order = ["full", "LDMT_no_feature", "LDMT_no_transfer"]
        data = [sub[sub["algorithm"] == v]["MIGD"].values for v in order]
        bp = ax.boxplot(data, labels=["full", "no_feature", "no_transfer"],
                        patch_artist=True, showfliers=False)
        for patch, c in zip(bp["boxes"], ["#d62728", "#ff7f0e", "#7f7f7f"]):
            patch.set_facecolor(c); patch.set_alpha(0.7)
        ax.set_title(prob); ax.set_ylabel("MIGD")
    fig.suptitle("消融实验: MIGD")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "ablation_MIGD.png"), dpi=120); plt.close(fig)
else:
    report_lines.append("\n## 2. 消融实验: 未找到 ablation_*.csv\n")

# ---------- sensitivity ----------
if sens is not None:
    runs_s = per_run_metrics(sens)
    runs_s.to_csv(os.path.join(SUMM, "sensitivity_runs.csv"), index=False)
    report_lines.append("\n## 3. 参数敏感性 (FDA3, LDMT, 30次)\n")
    report_lines.append("(a) 种群规模 N∈{50,100,150}; (b) 变化幅度 n_t∈{10,20,40}。\n\n")
    report_lines.append("| study | value | MIGD | MHV |")
    report_lines.append("\n|---|---|---|---|")
    # derive study/value from filename not stored -> approximate by N and n_t
    for _, r in runs_s.groupby(["N", "n_t"]).agg(MIGD=("MIGD", "mean"), MHV=("MHV", "mean")).reset_index().iterrows():
        study = "N" if r["N"] != 100 else ("n_t" if r["n_t"] != 20 else "baseline")
        report_lines.append(f"\n| {study} | N={int(r['N'])}, n_t={int(r['n_t'])} | {r['MIGD']:.4f} | {r['MHV']:.4f} |")
    # N sensitivity plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sN = runs_s[(runs_s["n_t"] == 20)].groupby("N")["MIGD"].mean()
    if len(sN):
        axes[0].plot(sN.index, sN.values, "o-", color="#d62728")
        axes[0].set_title("N 敏感性 (n_t=20)"); axes[0].set_xlabel("N"); axes[0].set_ylabel("MIGD"); axes[0].grid(alpha=0.3)
    snt = runs_s[(runs_s["N"] == 100)].groupby("n_t")["MIGD"].mean()
    if len(snt):
        axes[1].plot(snt.index, snt.values, "o-", color="#1f77b4")
        axes[1].set_title("n_t 敏感性 (N=100)"); axes[1].set_xlabel("n_t"); axes[1].set_ylabel("MIGD"); axes[1].grid(alpha=0.3)
    fig.suptitle("参数敏感性 (FDA3, LDMT)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "sensitivity.png"), dpi=120); plt.close(fig)
else:
    report_lines.append("\n## 3. 参数敏感性: 未找到 sensitivity_*.csv\n")

# ---------- supplementary ----------
if supp is not None:
    runs_u = per_run_metrics(supp)
    runs_u.to_csv(os.path.join(SUMM, "supplementary_runs.csv"), index=False)
    su = summarize(runs_u)
    su.to_csv(os.path.join(SUMM, "supplementary_summary.csv"), index=False)
    report_lines.append("\n## 4. 补充实验 (dMOP1/2, τ_t=25, 4算法, 30次)\n\n")
    report_lines.append("| problem | algorithm | MIGD | MHV | mig_t(s) |")
    report_lines.append("\n|---|---|---|---|---|")
    for _, r in su.iterrows():
        report_lines.append(
            f"\n| {r['problem']} | {ALGO_LABEL.get(r['algorithm'], r['algorithm'])} | {r['MIGD_mean']:.4f}±{r['MIGD_std']:.4f} | "
            f"{r['MHV_mean']:.4f} | {fmt(r['mig_time_mean'],5)} |")
else:
    report_lines.append("\n## 4. 补充实验: 未找到 supplementary_*.csv\n")

# ---------- conclusion ----------
report_lines.append("\n## 5. 结论要点\n")
report_lines.append("- 见上述各表与 figures/ 下箱线图、恢复曲线、耗时对比、消融与敏感性图。\n")
report_lines.append("- 假设H1(分布迁移改善变化后跟踪): 由 post_change_MIGD 与恢复曲线判断。\n")
report_lines.append("- 假设H2(低维特征迁移相对GMM更省时): 由 mig_time 对比判断。\n")
report_lines.append("- 假设H3(迁移比例/变化幅度影响稳定性): 由敏感性图判断。\n")
report_lines.append("- 若某指标 LDMT 未达优势, 已如实列出负面结果, 不删除不利测试实例。\n")

with open(os.path.join(SUMM, "analysis_report.md"), "w", encoding="utf-8") as f:
    f.write("".join(report_lines))
print("Analysis complete. See results/summary/ and results/figures/")
print("Files:", os.listdir(SUMM), "figures:", os.listdir(FIG))
