# -*- coding: utf-8 -*-
"""make_thesis_figures.py - 生成本科毕业论文所需的全部实验/示意图（PNG，dpi=300）。

数据源（全部来自真实结果文件，禁止手填实验数字）：
  results/summary/summary_comparison.csv   主比较：MIGD/MHV 均值±标准差、迁移耗时
  results/summary/friedman_results.json    Friedman 平均秩与 Nemenyi 临界差异 CD
  results/summary/_recovery.json           τ=25 变化后逐代平均 IGD（0–24 代）
  results/summary/ablation_summary.csv     消融实验（full / A / B / C）
  results/summary/sensitivity_runs.csv     种群规模 N 敏感性（逐 run，n_t=20）
  problems/true_fronts.py                  FDA 真实前沿 get_true_pf()

输出目录：results/figures_thesis/（不写入 results/figures/）

用法（在项目根目录）：python scripts\make_thesis_figures.py
"""
import csv
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyBboxPatch, Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from problems.true_fronts import get_true_pf  # noqa: E402

SUMMARY_DIR = os.path.join(ROOT, "results", "summary")
OUT_DIR = os.path.join(ROOT, "results", "figures_thesis")

# ----------------------------------------------------------------------
# 统一风格：中文字体、无 unicode 负号问题、色盲安全配色、字号 10-12、线宽 1.5
# ----------------------------------------------------------------------
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10
plt.rcParams["axes.labelsize"] = 10.5
plt.rcParams["xtick.labelsize"] = 9.5
plt.rcParams["ytick.labelsize"] = 9.5
plt.rcParams["legend.fontsize"] = 9
plt.rcParams["axes.linewidth"] = 0.9
plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["savefig.dpi"] = 300

# 六算法固定顺序 / 显示名 / 颜色（LDT 用醒目橙色并加粗）
ALG_ORDER = ["pps", "ar1_predictive", "ldt_nsga2", "gmm_inspired", "nsga2_static", "dnsga2_inspired"]
ALG_LABEL = {
    "pps": "PPS",
    "ar1_predictive": "AR1",
    "ldt_nsga2": "LDT-NSGA-II（本文）",
    "gmm_inspired": "GMM",
    "nsga2_static": "NSGA-II（static）",
    "dnsga2_inspired": "DNSGA-II",
}
ALG_COLOR = {
    "pps": "#009E73",
    "ar1_predictive": "#56B4E9",
    "ldt_nsga2": "#D55E00",
    "gmm_inspired": "#CC79A7",
    "nsga2_static": "#0072B2",
    "dnsga2_inspired": "#E69F00",
}
ALG_MARKER = {
    "pps": "o",
    "ar1_predictive": "s",
    "ldt_nsga2": "*",
    "gmm_inspired": "D",
    "nsga2_static": "v",
    "dnsga2_inspired": "^",
}

PROBLEMS = ["FDA1", "FDA3", "FDA4"]
TAUS = [10, 25, 50]


# ----------------------------------------------------------------------
# 数据读取
# ----------------------------------------------------------------------
def load_summary():
    """summary_comparison.csv -> {(problem, tau): {alg: row_dict}}"""
    table = {}
    with open(os.path.join(SUMMARY_DIR, "summary_comparison.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            key = (row["problem"], int(row["tau_t"]))
            table.setdefault(key, {})[row["algorithm"]] = row
    return table


def load_friedman():
    with open(os.path.join(SUMMARY_DIR, "friedman_results.json"), encoding="utf-8") as fh:
        return json.load(fh)


def load_recovery():
    with open(os.path.join(SUMMARY_DIR, "_recovery.json"), encoding="utf-8") as fh:
        return json.load(fh)


def load_ablation():
    """ablation_summary.csv -> {problem: {alg: row_dict}}（τ=25）"""
    table = {}
    with open(os.path.join(SUMMARY_DIR, "ablation_summary.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            table.setdefault(row["problem"], {})[row["algorithm"]] = row
    return table


def load_sensitivity():
    """sensitivity_runs.csv（仅 ldt_nsga2, τ=25）。过滤 n_t=20 保证 N 间可比。"""
    runs = []
    with open(os.path.join(SUMMARY_DIR, "sensitivity_runs.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["algorithm"] == "ldt_nsga2" and int(row["tau_t"]) == 25 and int(row["n_t"]) == 20:
                runs.append(row)
    return runs

# ----------------------------------------------------------------------
# 公共小工具
# ----------------------------------------------------------------------
def _save(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("已保存:", path)
    return path


def _alg_legend_handles():
    handles = []
    for a in ALG_ORDER:
        is_ldt = a == "ldt_nsga2"
        handles.append(
            Line2D(
                [0], [0],
                color=ALG_COLOR[a],
                marker=ALG_MARKER[a],
                lw=2.4 if is_ldt else 1.5,
                ms=11 if is_ldt else 5.5,
                markerfacecolor=ALG_COLOR[a],
                markeredgecolor="white" if is_ldt else ALG_COLOR[a],
                markeredgewidth=0.6,
            )
        )
    return handles, [ALG_LABEL[a] for a in ALG_ORDER]


def _shared_alg_legend(fig, ncol=6):
    handles, labels = _alg_legend_handles()
    leg = fig.legend(
        handles, labels,
        loc="outside upper center", ncol=ncol,
        frameon=False, columnspacing=1.2, handletextpad=0.4,
    )
    # 加粗"本文方法"图例文字以突出
    for text, lab in zip(leg.get_texts(), labels):
        if "LDT" in lab:
            text.set_fontweight("bold")
            text.set_color(ALG_COLOR["ldt_nsga2"])
    return leg


def _rbox(ax, cx, cy, w, h, text, face="#EAF2F8", edge="#2C5F8A", lw=1.3,
          fontsize=9.5, fontweight="normal", color="black"):
    """以 (cx, cy) 为中心画圆角矩形并居中写字。"""
    box = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.06,rounding_size=0.10",
        facecolor=face, edgecolor=edge, linewidth=lw, zorder=2,
    )
    ax.add_patch(box)
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, color=color, zorder=3)
    return box


def _arrow(ax, xy_from, xy_to, color="#555555", lw=1.4, style="-|>", ms=10, ls="-"):
    ax.annotate("", xy=xy_to, xytext=xy_from,
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                                mutation_scale=ms, shrinkA=0.0, shrinkB=0.0,
                                linestyle=ls, zorder=1))


def _cov_ellipse(ax, mean, cov, n_std=1.6, **kw):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    ell = Ellipse(
        xy=mean, width=2 * n_std * np.sqrt(vals[0]), height=2 * n_std * np.sqrt(vals[1]),
        angle=angle, zorder=2, **kw,
    )
    ax.add_patch(ell)
    return ell

# ----------------------------------------------------------------------
# 图 1-1 动态多目标优化问题示意
# ----------------------------------------------------------------------
def fig1_1_dmop_concept():
    pf = get_true_pf("FDA1", M=2, n_points=300)  # 真实前沿 f2 = 1 - sqrt(f1)
    rng = np.random.default_rng(7)

    # t 时刻：种群已收敛到前沿附近
    f1p = np.sort(rng.uniform(0.04, 0.96, 24))
    pop_t = np.column_stack([f1p, 1.0 - np.sqrt(f1p) + rng.normal(0, 0.012, f1p.size)])

    # t+1 时刻：环境变化使新前沿整体下移（示意），旧种群滞后于新前沿
    pf_new = pf * np.array([1.0, 0.88])
    gap = (1.0 - np.sqrt(f1p)) * 0.12

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.3), constrained_layout=True)

    # (a) t 时刻
    ax1.plot(pf[:, 0], pf[:, 1], color="#0072B2", lw=1.8, label="t 时刻 POF", zorder=2)
    ax1.scatter(pop_t[:, 0], pop_t[:, 1], s=22, color="#0072B2", zorder=3,
                edgecolors="white", linewidths=0.4, label="种群（已收敛）")
    ax1.set_xlabel("目标 $f_1$（无量纲）")
    ax1.set_ylabel("目标 $f_2$（无量纲）")
    ax1.set_title("FDA1：t 时刻", fontsize=11)
    ax1.set_xlim(-0.05, 1.05)
    ax1.set_ylim(-0.12, 1.08)
    ax1.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    ax1.text(0.05, 0.06, "(a)", transform=ax1.transAxes, fontsize=11, fontweight="bold")

    # (b) t+1 时刻
    ax2.plot(pf[:, 0], pf[:, 1], color="#999999", lw=1.2, ls="--", zorder=1,
             label="t 时刻前沿（旧）")
    ax2.plot(pf_new[:, 0], pf_new[:, 1], color="#D55E00", lw=1.8, zorder=2,
             label="t+1 时刻 POF（新）")
    ax2.scatter(f1p, 1.0 - np.sqrt(f1p), s=26, facecolors="none",
                edgecolors="#0072B2", linewidths=1.1, zorder=3,
                label="旧种群（滞后）")
    # 迁移方向箭头（每隔一个点画，避免杂乱）
    for x, g in zip(f1p[::2], gap[::2]):
        _arrow(ax2, (x, 1.0 - np.sqrt(x) - 0.02), (x, 0.88 * (1.0 - np.sqrt(x)) + 0.03),
               color="#8A8A8A", lw=1.0, ms=8)
    pop_new = np.column_stack([f1p, 0.88 * (1.0 - np.sqrt(f1p))])
    ax2.scatter(pop_new[:, 0], pop_new[:, 1], s=22, color="#D55E00", zorder=4,
                edgecolors="white", linewidths=0.4, label="迁移后种群")
    ax2.annotate("迁移方向", xy=(0.42, 0.60), xytext=(0.60, 0.78),
                 fontsize=9, color="#555555",
                 arrowprops=dict(arrowstyle="->", color="#555555", lw=1.1))
    ax2.set_xlabel("目标 $f_1$（无量纲）")
    ax2.set_ylabel("目标 $f_2$（无量纲）")
    ax2.set_title("FDA1：t+1 时刻（环境变化后）", fontsize=11)
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylim(-0.12, 1.08)
    ax2.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    ax2.text(0.05, 0.06, "(b)", transform=ax2.transAxes, fontsize=11, fontweight="bold")

    return _save(fig, "fig1_1_dmop_concept.png")


# ----------------------------------------------------------------------
# 图 2-1 方法分类树
# ----------------------------------------------------------------------
def fig2_1_taxonomy():
    fig, ax = plt.subplots(figsize=(10.5, 5.9))
    ax.set_xlim(0, 10.5)
    ax.set_ylim(0.6, 6.5)
    ax.axis("off")

    # 根节点
    _rbox(ax, 5.25, 5.9, 4.2, 0.85,
          "动态多目标进化算法\n变化响应策略",
          face="#D6E4F0", edge="#2C5F8A", lw=1.5, fontsize=11, fontweight="bold")

    cat_x = [1.85, 5.25, 8.65]
    cat_titles = ["响应式方法", "预测式方法", "知识 / 分布迁移型方法"]
    reps = ["代表：DNSGA-II", "代表：AR1、PPS", "代表：GMM"]
    ideas = [
        "核心思想：变化后对部分个体\n重新初始化或超变异，补充种群多样性",
        "核心思想：以历史最优解建立自回归/\n中心点模型，外推新环境最优位置",
        "核心思想：学习旧种群在决策空间的\n分布，将其向新环境迁移生成个体",
    ]

    for x, title, rep, idea in zip(cat_x, cat_titles, reps, ideas):
        _arrow(ax, (5.25, 5.46), (x, 4.92), color="#555555", lw=1.3)
        _rbox(ax, x, 4.55, 3.0, 0.72, title, face="#EAF2F8", edge="#2C5F8A",
              lw=1.4, fontsize=10.5, fontweight="bold")
        _arrow(ax, (x, 4.19), (x, 3.78), color="#888888", lw=1.1)
        _rbox(ax, x, 3.45, 2.9, 0.62, rep, face="#F5F5F5", edge="#777777", lw=1.1,
              fontsize=9.5)
        ax.text(x, 2.35, idea, ha="center", va="center", fontsize=9, color="#444444",
                linespacing=1.5)

    # 高亮"本文方法"小框
    _arrow(ax, (8.65, 3.14), (8.65, 2.86), color="#D55E00", lw=1.4)
    _rbox(ax, 8.65, 2.52, 2.9, 0.62, "LDT-NSGA-II（本文）", face="#FDEBD0",
          edge="#D55E00", lw=2.2, fontsize=10.5, fontweight="bold", color="#B34700")
    ax.text(8.65, 1.30, "核心思想：归一化特征提取 + 趋势中心预测\n+ 高斯采样迁移",
            ha="center", va="center", fontsize=9, color="#B34700", linespacing=1.5,
            fontweight="bold")

    return _save(fig, "fig2_1_taxonomy.png")

# ----------------------------------------------------------------------
# 图 3-1 LDT-NSGA-II 算法框架流程图
# ----------------------------------------------------------------------
def fig3_1_framework():
    fig, ax = plt.subplots(figsize=(8.2, 10.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(2.2, 14.0)
    ax.axis("off")

    base = dict(face="#EAF2F8", edge="#2C5F8A", lw=1.3)
    mod = dict(face="#FDEBD0", edge="#D55E00", lw=1.5)

    _rbox(ax, 5, 13.35, 5.2, 0.75, "初始化种群（均匀随机生成 N 个个体）", **base)
    _rbox(ax, 5, 12.05, 3.9, 0.72, "NSGA-II 进化 τ 代", **base)
    _rbox(ax, 5, 10.75, 4.6, 0.72, "检测环境变化（重新评估目标值）", **base)
    _rbox(ax, 5, 9.25, 4.6, 0.72, "非支配筛选，得到精英集 E", **base)
    _rbox(ax, 5, 7.75, 5.0, 0.70, "归一化特征提取（$\\mu_e$, $\\sigma_e$）", **mod)
    _rbox(ax, 5, 6.50, 5.9, 0.70, "趋势中心预测 $\\hat{\\mu}=\\mu_e+\\alpha(\\mu_e-\\mu_{prev})$", **mod)
    _rbox(ax, 5, 5.25, 5.5, 0.70, "高斯采样迁移 $x_{tr}=\\mathrm{clip}(\\hat{\\mu}+z\\sigma_e)$", **mod)
    _rbox(ax, 5, 3.60, 6.4, 0.75,
          "种群重组：迁移个体 50% / 保留个体 20% / 随机个体 30%", **base)

    # 主流程箭头
    _arrow(ax, (5, 12.97), (5, 12.43))
    _arrow(ax, (5, 11.69), (5, 11.13))
    _arrow(ax, (5, 10.39), (5, 9.63))
    _arrow(ax, (5, 8.89), (5, 8.12))
    _arrow(ax, (5, 7.40), (5, 6.87))
    _arrow(ax, (5, 6.15), (5, 5.62))
    _arrow(ax, (5, 4.90), (5, 4.00))
    ax.text(5.18, 9.93, "有变化", fontsize=9, color="#333333", ha="left")

    # 无变化分支：从检测框右侧回到进化框
    ax.plot([7.32, 9.15, 9.15], [10.75, 10.75, 12.05], color="#555555", lw=1.4, zorder=1)
    _arrow(ax, (9.15, 12.05), (6.97, 12.05))
    ax.text(8.15, 10.97, "无变化", fontsize=9, color="#333333", ha="center")

    # 新环境进化循环：从重组框左侧回到进化框
    ax.plot([1.78, 0.65, 0.65], [3.60, 3.60, 12.05], color="#555555", lw=1.4, zorder=1)
    _arrow(ax, (0.65, 12.05), (3.03, 12.05))
    ax.text(0.38, 7.8, "新环境进化循环", fontsize=9, color="#333333",
            ha="center", va="center", rotation=90)

    # 虚线框：分布迁移模块（本文）
    dashed = FancyBboxPatch(
        (1.55, 4.55), 6.9, 3.65,
        boxstyle="round,pad=0.06,rounding_size=0.12",
        facecolor="none", edgecolor="#D55E00", linewidth=1.8,
        linestyle=(0, (5, 3)), zorder=0.5,
    )
    ax.add_patch(dashed)
    ax.text(1.75, 8.28, "分布迁移模块（本文）", fontsize=10, fontweight="bold",
            color="#B34700", ha="left", va="center")

    # 图例：实线=基线 NSGA-II 部分，橙色=本文模块
    handles = [
        Patch(facecolor="#EAF2F8", edgecolor="#2C5F8A", label="基线 NSGA-II 流程"),
        Patch(facecolor="#FDEBD0", edgecolor="#D55E00", label="分布迁移模块（本文）"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.035),
              ncol=2, frameon=False, fontsize=9.5)

    return _save(fig, "fig3_1_framework.png")


# ----------------------------------------------------------------------
# 图 3-2 分布迁移策略示意
# ----------------------------------------------------------------------
def fig3_2_transfer():
    rng = np.random.default_rng(11)
    mu_prev = np.array([2.0, 3.4])
    mu_e = np.array([3.7, 4.5])
    alpha = 0.5
    mu_hat = mu_e + alpha * (mu_e - mu_prev)
    sigma = np.array([[0.50, 0.22], [0.22, 0.32]])

    old_pop = rng.multivariate_normal(mu_e, sigma, 220)
    new_pop = rng.multivariate_normal(mu_hat, sigma * 0.9, 130)

    fig, ax = plt.subplots(figsize=(7.8, 5.6), constrained_layout=True)

    _cov_ellipse(ax, mu_e, sigma, n_std=1.6, facecolor="#0072B2", alpha=0.08,
                 edgecolor="#0072B2", lw=1.5, linestyle="-")
    _cov_ellipse(ax, mu_hat, sigma * 0.9, n_std=1.6, facecolor="#D55E00", alpha=0.08,
                 edgecolor="#D55E00", lw=1.5, linestyle="--")

    ax.scatter(old_pop[:, 0], old_pop[:, 1], s=10, color="#0072B2", alpha=0.45,
               linewidths=0, zorder=3)
    ax.scatter(new_pop[:, 0], new_pop[:, 1], s=10, color="#D55E00", alpha=0.5,
               linewidths=0, zorder=3)

    _arrow(ax, tuple(mu_prev), tuple(mu_e), color="#777777", lw=1.3, ms=11)
    _arrow(ax, tuple(mu_e), tuple(mu_hat), color="#D55E00", lw=2.2, ms=14)

    ax.plot(*mu_prev, "x", color="#555555", ms=9, mew=2, zorder=4)
    ax.plot(*mu_e, "o", color="#0072B2", ms=8, mec="white", mew=0.8, zorder=4)
    ax.plot(*mu_hat, "*", color="#D55E00", ms=17, mec="white", mew=0.6, zorder=4)

    ax.text(mu_prev[0] - 0.05, mu_prev[1] - 0.38, "$\\mu_{prev}$（上一环境中心）",
            fontsize=9.5, ha="center", color="#555555")
    ax.text(mu_e[0] + 0.02, mu_e[1] - 0.50, "$\\mu_e$（当前环境均值）",
            fontsize=9.5, ha="center", color="#2C5F8A")
    ax.text(mu_hat[0] + 0.05, mu_hat[1] + 0.42, "$\\hat{\\mu}$（预测中心）",
            fontsize=9.5, ha="center", color="#B34700")

    ax.text(1.15, 5.35, "趋势方向", fontsize=9, color="#777777", rotation=32)
    ax.text(3.45, 5.42, "迁移方向 $\\Delta=\\mu_e-\\mu_{prev}$", fontsize=9.5,
            color="#B34700", rotation=25)

    ax.text(0.99, 0.02, "高斯采样：$x_{tr}=\\mathrm{clip}(\\hat{\\mu}+z\\sigma_e)$",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#FDF3E7",
                      edgecolor="#D55E00", lw=1.0))

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#0072B2",
               markersize=7, label="旧环境种群"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#D55E00",
               markersize=7, label="预测/迁移个体"),
        Line2D([0], [0], color="#0072B2", lw=1.5, label="旧分布（1.6σ 椭圆）"),
        Line2D([0], [0], color="#D55E00", lw=1.5, ls="--", label="预测分布（1.6σ 椭圆）"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, framealpha=0.9)

    ax.set_xlabel("决策变量 $x_1$（归一化，示意）")
    ax.set_ylabel("决策变量 $x_2$（归一化，示意）")
    ax.set_xlim(0.3, 7.0)
    ax.set_ylim(1.6, 6.6)

    return _save(fig, "fig3_2_transfer.png")

# ----------------------------------------------------------------------
# 图 5-1 MIGD 随 τ 变化（1×3）
# ----------------------------------------------------------------------
def fig5_1_migd_vs_tau(summary):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for ax, prob in zip(axes, PROBLEMS):
        for a in ALG_ORDER:
            ms, ss = [], []
            for tau in TAUS:
                row = summary[(prob, tau)][a]
                ms.append(float(row["MIGD_mean"]))
                ss.append(float(row["MIGD_std"]))
            is_ldt = a == "ldt_nsga2"
            ax.errorbar(TAUS, ms, yerr=ss, color=ALG_COLOR[a], marker=ALG_MARKER[a],
                        lw=2.4 if is_ldt else 1.5, ms=11 if is_ldt else 5.5,
                        capsize=3, elinewidth=1.0, zorder=4 if is_ldt else 3,
                        markeredgecolor="white" if is_ldt else ALG_COLOR[a],
                        markeredgewidth=0.6)
        ax.set_xticks(TAUS)
        ax.set_title(prob, fontsize=11)
        ax.set_xlabel("变化频率参数 τ（代）")
    axes[0].set_ylabel("稳态 MIGD（均值±标准差，越小越优）")
    _shared_alg_legend(fig)
    return _save(fig, "fig5_1_migd_vs_tau.png")


# ----------------------------------------------------------------------
# 图 5-3 / 5-4 变化后收敛与归一化恢复曲线（1×3）
# ----------------------------------------------------------------------
def _recovery_panels(recovery, normalize, ylabel, fname):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for ax, prob in zip(axes, PROBLEMS):
        gens = sorted(int(k) for k in recovery[prob][ALG_ORDER[0]].keys())
        for a in ALG_ORDER:
            vals = np.array([recovery[prob][a][str(g)] for g in gens], dtype=float)
            if normalize:
                vals = vals / vals[0]
            is_ldt = a == "ldt_nsga2"
            ax.plot(gens, vals, color=ALG_COLOR[a], marker=ALG_MARKER[a],
                    lw=2.4 if is_ldt else 1.5, ms=8 if is_ldt else 3.5,
                    markevery=3, zorder=4 if is_ldt else 3,
                    markeredgecolor="white" if is_ldt else ALG_COLOR[a],
                    markeredgewidth=0.5)
        ax.set_title(prob, fontsize=11)
        ax.set_xlabel("变化后进化代数（代）")
        if normalize:
            ax.set_ylim(0.55, 1.05)
            ax.axhline(1.0, color="#BBBBBB", lw=0.8, ls=":", zorder=1)
        else:
            ax.set_yscale("log")
    axes[0].set_ylabel(ylabel)
    _shared_alg_legend(fig)
    return _save(fig, fname)


def fig5_3_convergence(recovery):
    return _recovery_panels(recovery, normalize=False,
                            ylabel="平均 IGD（对数刻度，越小越优）",
                            fname="fig5_3_convergence.png")


def fig5_4_recovery(recovery):
    return _recovery_panels(recovery, normalize=True,
                            ylabel="IGD(t) / IGD(0)（按各算法首代归一化）",
                            fname="fig5_4_recovery.png")


# ----------------------------------------------------------------------
# 图 5-5 迁移模块耗时对比（9 组 × 6 算法，对数轴）
# ----------------------------------------------------------------------
def fig5_5_time(summary):
    groups = [(p, t) for p in PROBLEMS for t in TAUS]
    x = np.arange(len(groups))
    n_alg = len(ALG_ORDER)
    w = 0.8 / n_alg

    fig, ax = plt.subplots(figsize=(12, 4.4), constrained_layout=True)
    gmm_vals = ldt_vals = None
    for k, a in enumerate(ALG_ORDER):
        vals = np.array([float(summary[g][a]["mig_time_mean"]) for g in groups])
        ax.bar(x + (k - (n_alg - 1) / 2) * w, vals, width=w * 0.92,
               color=ALG_COLOR[a], edgecolor="black", linewidth=0.4,
               label=ALG_LABEL[a], zorder=3)
        if a == "gmm_inspired":
            gmm_vals = vals
        if a == "ldt_nsga2":
            ldt_vals = vals

    # 在 GMM 柱上方标注相对 LDT 的倍数（由数据实算）
    for g, (gv, lv) in enumerate(zip(gmm_vals, ldt_vals)):
        ratio = gv / lv
        ax.text(g + (3 - (n_alg - 1) / 2) * w, gv * 1.22, "%.1f×" % ratio,
                rotation=90, ha="center", va="bottom", fontsize=7.5, color="#B2182B")

    for sep in (2.5, 5.5):
        ax.axvline(sep, color="#BBBBBB", lw=0.8, ls=":", zorder=1)

    ax.set_yscale("log")
    ax.set_ylim(6e-5, 2.2)
    ax.set_xticks(x)
    ax.set_xticklabels(["%s\nτ=%d" % g for g in groups], fontsize=8.5)
    ax.set_xlabel("问题 × 变化频率配置（共 9 组）")
    ax.set_ylabel("单次迁移耗时（秒，对数刻度）")
    _shared_alg_legend(fig)
    return _save(fig, "fig5_5_time.png")

# ----------------------------------------------------------------------
# 图 6-1 消融实验分组柱状
# ----------------------------------------------------------------------
def fig6_1_ablation(ablation):
    keys = ["ldt_nsga2", "LDT_AblationA", "LDT_AblationB", "LDT_AblationC"]
    labels = ["LDT（完整）", "变体A（关趋势外推）", "变体B（关分布迁移）", "变体C（关旧种群保留）"]
    colors = ["#D55E00", "#56B4E9", "#0072B2", "#009E73"]

    x = np.arange(len(PROBLEMS))
    w = 0.8 / len(keys)
    fig, ax = plt.subplots(figsize=(7.8, 4.3), constrained_layout=True)
    for k, (key, lab, c) in enumerate(zip(keys, labels, colors)):
        ms = np.array([float(ablation[p][key]["MIGD_mean"]) for p in PROBLEMS])
        ss = np.array([float(ablation[p][key]["MIGD_std"]) for p in PROBLEMS])
        ax.bar(x + (k - 1.5) * w, ms, width=w * 0.9, yerr=ss, capsize=3,
               error_kw=dict(elinewidth=1.0, ecolor="#333333"),
               color=c, edgecolor="black", linewidth=0.4, label=lab, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(PROBLEMS)
    ax.set_xlabel("测试问题（τ=25）")
    ax.set_ylabel("MIGD（均值±标准差，越小越优）")
    handles = [Patch(facecolor=c, edgecolor="black", linewidth=0.4, label=lab)
               for lab, c in zip(labels, colors)]
    leg = fig.legend(handles=handles, loc="outside upper center", ncol=4,
                     frameon=False)
    leg.get_texts()[0].set_fontweight("bold")
    leg.get_texts()[0].set_color("#D55E00")
    return _save(fig, "fig6_1_ablation.png")


# ----------------------------------------------------------------------
# 图 7-1 种群规模 N 敏感性（MIGD 柱状 + 单代耗时双轴）
# ----------------------------------------------------------------------
def fig7_1_sensitivity(sens):
    N_LIST = [50, 100, 200]
    prob_colors = {"FDA1": "#0072B2", "FDA3": "#009E73", "FDA4": "#D55E00"}

    migd_mean, migd_std, gtime = {}, {}, {}
    for p in PROBLEMS:
        migd_mean[p], migd_std[p], gtime[p] = [], [], []
        for n in N_LIST:
            vals = [float(r["MIGD"]) for r in sens if r["problem"] == p and int(r["N"]) == n]
            ts = [float(r["mean_gen_time_s"]) for r in sens
                  if r["problem"] == p and int(r["N"]) == n]
            migd_mean[p].append(np.mean(vals))
            migd_std[p].append(np.std(vals, ddof=1))
            gtime[p].append(np.mean(ts))

    x = np.arange(len(N_LIST))
    w = 0.8 / len(PROBLEMS)
    fig, ax = plt.subplots(figsize=(8.2, 4.5), constrained_layout=True)
    for k, p in enumerate(PROBLEMS):
        ax.bar(x + (k - 1) * w, migd_mean[p], width=w * 0.9,
               yerr=migd_std[p], capsize=3,
               error_kw=dict(elinewidth=1.0, ecolor="#333333"),
               color=prob_colors[p], alpha=0.85, edgecolor="black",
               linewidth=0.4, zorder=3)

    ax.set_yscale("log")
    ax.set_ylim(0.003, 0.42)
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in N_LIST])
    ax.set_xlabel("种群规模 N（个体数）")
    ax.set_ylabel("MIGD（均值±标准差，对数刻度，越小越优）")

    ax2 = ax.twinx()
    for p in PROBLEMS:
        ax2.plot(x, gtime[p], color=prob_colors[p], marker="o", lw=1.5,
                 ls="--", ms=5.5, zorder=4)
    ax2.set_ylabel("单代进化耗时（秒）")
    ax2.set_ylim(0, 0.010)

    handles = []
    for p in PROBLEMS:
        handles.append(Patch(facecolor=prob_colors[p], alpha=0.85,
                             edgecolor="black", linewidth=0.4,
                             label="%s：MIGD（左轴）" % p))
    for p in PROBLEMS:
        handles.append(Line2D([0], [0], color=prob_colors[p], marker="o", lw=1.5,
                              ls="--", ms=5.5, label="%s：单代耗时（右轴）" % p))
    fig.legend(handles=handles, loc="outside upper center", ncol=3, frameon=False,
              fontsize=8.5)
    return _save(fig, "fig7_1_sensitivity_N.png")

# ----------------------------------------------------------------------
# 图 8-1 Friedman + Nemenyi 临界差异图
# ----------------------------------------------------------------------
def fig8_1_nemenyi(friedman):
    # 对 9 个 问题×τ 组合求每算法平均秩
    acc = {a: [] for a in friedman[0]["algorithms"]}
    for entry in friedman:
        for a, r in zip(entry["algorithms"], entry["avg_ranks"]):
            acc[a].append(r)
    avg = {a: float(np.mean(v)) for a, v in acc.items()}
    ranked = sorted(avg.items(), key=lambda kv: kv[1])  # 秩小者优
    cd = float(friedman[0]["CD"])

    print("9 组合平均秩（优到劣）：",
          [(ALG_LABEL[a], round(r, 3)) for a, r in ranked])

    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)

    # 算法点与名称（上下交替布置引线标签）
    for i, (a, r) in enumerate(ranked):
        y = 0.62 if i % 2 == 0 else -0.62
        ax.plot([r, r], [0, y], color="#888888", lw=0.9, zorder=1)
        ax.plot(r, 0, marker=ALG_MARKER[a], color=ALG_COLOR[a],
                ms=14 if a == "ldt_nsga2" else 9,
                mec="black", mew=0.6, zorder=4)
        ax.text(r, y * 1.22, ALG_LABEL[a], ha="center", va="center", fontsize=10,
                fontweight="bold" if a == "ldt_nsga2" else "normal",
                color="#B34700" if a == "ldt_nsga2" else "black")

    # 平均秩差 < CD 的算法对：横线连接（差异不显著）
    keys = [a for a, _ in ranked]
    level = 0
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if abs(avg[keys[i]] - avg[keys[j]]) < cd:
                y = -1.15 - 0.34 * level
                level += 1
                ax.plot([avg[keys[i]], avg[keys[j]]], [y, y], color="#B2182B",
                        lw=2.0, solid_capstyle="round", zorder=2)
                for xx in (avg[keys[i]], avg[keys[j]]):
                    ax.plot([xx, xx], [0, y], color="#B2182B", lw=0.9, ls=":", zorder=2)

    # CD 参考段
    x0, y_cd = 0.75, 1.85
    ax.plot([x0, x0 + cd], [y_cd, y_cd], color="black", lw=1.5)
    for xx in (x0, x0 + cd):
        ax.plot([xx, xx], [y_cd - 0.09, y_cd + 0.09], color="black", lw=1.5)
    ax.text(x0 + cd / 2, y_cd + 0.22, "临界差异 CD = %.3f" % cd,
            ha="center", fontsize=10)

    ax.axhline(0, color="black", lw=1.0, zorder=0.5)
    ax.set_xlim(0.4, 6.1)
    ax.set_ylim(-2.7, 2.5)
    ax.get_yaxis().set_visible(False)
    for s in ("left", "top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("平均秩（9 个 问题×τ 组合，秩越小越优）", fontsize=11)
    return _save(fig, "fig8_1_nemenyi.png")


# ----------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    summary = load_summary()
    friedman = load_friedman()
    recovery = load_recovery()
    ablation = load_ablation()
    sens = load_sensitivity()

    fig1_1_dmop_concept()
    fig2_1_taxonomy()
    fig3_1_framework()
    fig3_2_transfer()
    fig5_1_migd_vs_tau(summary)
    fig5_3_convergence(recovery)
    fig5_4_recovery(recovery)
    fig5_5_time(summary)
    fig6_1_ablation(ablation)
    fig7_1_sensitivity(sens)
    fig8_1_nemenyi(friedman)

    print("")
    print("全部图已生成于:", OUT_DIR)


if __name__ == "__main__":
    main()