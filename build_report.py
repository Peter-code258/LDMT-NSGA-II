#!/usr/bin/env python3
"""Assemble the LDT experiment report HTML with all data embedded (v2.4).

v2.4 (K12 fix): doc.css and echarts.min.js are now INLINED at build time so
the generated index.html is a single self-contained file (previously the
report linked to assets/ relatively and broke when moved).
"""
import json
import os
import csv
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\23-兴湘计科\LDT-NSGA-II")
# v2.2: results tree + output path overridable via environment variables
RESULTS_ROOT = Path(os.environ.get("LDT_ROOT", str(ROOT)))
OUT_PATH = Path(os.environ.get(
    "LDT_OUT", str(ROOT / "ldmt-experiment-report" / "index.html")))
# v2.4 (K12): read assets for inlining (kept on disk as the source of truth)
_ASSETS = ROOT / "ldmt-experiment-report" / "assets"
DOC_CSS = (_ASSETS / "doc.css").read_text(encoding="utf-8")
ECHARTS_JS = (_ASSETS / "echarts.min.js").read_text(encoding="utf-8")
CD = json.loads((RESULTS_ROOT / "results/summary/_chart_data.json").read_text(encoding="utf-8-sig"))
REC = json.loads((RESULTS_ROOT / "results/summary/_recovery.json").read_text(encoding="utf-8-sig"))
BOX = json.loads((RESULTS_ROOT / "results/summary/_boxplot.json").read_text(encoding="utf-8-sig"))

# v2.1: Friedman + Nemenyi statistical test results
FRIED_PATH = RESULTS_ROOT / "results/summary/friedman_results.json"
FRIED = json.loads(FRIED_PATH.read_text(encoding="utf-8-sig")) if FRIED_PATH.exists() else []

# v2.1: recovery generations (per env recovery speed)
def _std(v):
    m = sum(v) / len(v)
    return (sum((x - m) ** 2 for x in v) / (len(v) - 1)) ** 0.5

RECG_PATH = RESULTS_ROOT / "results/summary/recovery_generations.csv"
recg_agg = defaultdict(list)
if RECG_PATH.exists():
    with open(RECG_PATH, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if int(row["tau_t"]) == 25:
                recg_agg[(row["problem"], row["algorithm"])].append(
                    float(row["recovery_gens"]))
RECG = [{"problem": p, "algo": a,
         "mean": round(sum(v) / len(v), 2),
         "std": round(_std(v), 2) if len(v) > 1 else 0.0}
        for (p, a), v in sorted(recg_agg.items())]

comp = CD["comparison"]
ablation = CD["ablation"]
sensitivity = CD["sensitivity"]

# mig_time: average over tau per (problem, algo)
from collections import defaultdict
mt = defaultdict(list)
for r in comp:
    mt[(r["problem"], r["algo"])].append(r["mig_t"])
mig_time = [{"problem": p, "algo": a, "mig_t": round(sum(v) / len(v), 5)}
            for (p, a), v in mt.items()]

def js(payload):
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))

# pre-serialize to avoid f-string brace escaping issues
COMP_JS = js(comp)
MIGT_JS = js(mig_time)
ABL_JS = js(ablation)
SENS_JS = js(sensitivity)
BOX1_JS = js(BOX.get("FDA1", {}))
BOX3_JS = js(BOX.get("FDA3", {}))
BOX4_JS = js(BOX.get("FDA4", {}))
REC1_JS = js(REC.get("FDA1", {}))
REC3_JS = js(REC.get("FDA3", {}))
REC4_JS = js(REC.get("FDA4", {}))

# Build comparison table rows (HTML) for main results
problems_order = ["FDA1", "FDA3", "FDA4", "dMOP1", "dMOP2", "DMOP3"]
# v2.2: canonical algorithm ids (experiments/common.py CANONICAL_ALGOS)
algo_order = ["nsga2_static", "dnsga2_inspired", "ar1_predictive",
              "gmm_inspired", "pps", "ldt_nsga2"]
# v2.2: display names mirror experiments/common.py ALGO_DISPLAY (hardcoded)
algo_label = {"nsga2_static": "NSGA-II (static)",
              "dnsga2_inspired": "DNSGA-II-inspired",
              "ar1_predictive": "AR1-Predictive",
              "gmm_inspired": "GMM-inspired",
              "pps": "PPS",
              "ldt_nsga2": "LDT-NSGA-II"}

# v2.1: build Friedman test table rows
fried_rows_html = []
for fr in FRIED:
    p = fr["p_value"]
    p_str = "<0.0001" if p < 1e-4 else f"{p:.4f}"
    sig = "显著" if p < 0.05 else "不显著"
    fried_rows_html.append(
        f"<tr><td>{fr['problem']}</td><td>{fr['tau_t']}</td>"
        f"<td>{fr['chi2']:.3f}</td><td>{p_str}</td>"
        f"<td>{fr['CD']:.3f}</td><td>{sig}</td></tr>")
fried_table = "\n".join(fried_rows_html)

# v2.1: build recovery generations table rows (tau=25)
recg_rows_html = []
for r in RECG:
    label = algo_label.get(r["algo"], r["algo"])
    recg_rows_html.append(
        f"<tr><td>{r['problem']}</td><td>{label}</td>"
        f"<td>{r['mean']:.2f}</td><td>±{r['std']:.2f}</td></tr>")
recg_table = "\n".join(recg_rows_html)

comp_by = {}
for r in comp:
    comp_by[(r["problem"], r["tau"], r["algo"])] = r

def best_flag(r, field, lower_better=True):
    vals = [comp_by[(r["problem"], r["tau"], a)][field]
            for a in algo_order if (r["problem"], r["tau"], a) in comp_by]
    if not vals:
        return False
    best = min(vals) if lower_better else max(vals)
    return abs(r[field] - best) < 1e-9

def fmt(v, dec=4):
    return f"{v:.{dec}f}"

# Build main results table
table_rows = []
for prob in problems_order:
    for tau in [10, 25, 50]:
        for algo in algo_order:
            r = comp_by.get((prob, tau, algo))
            if not r:
                continue
            mg_best = best_flag(r, "MIGD", True)
            hv_best = best_flag(r, "MHV", False)
            mg_cell = fmt(r["MIGD"]) + ("" if not mg_best else " †")
            hv_cell = fmt(r["MHV"]) + ("" if not hv_best else " †")
            table_rows.append(
                f"<tr><td>{prob}</td><td>{tau}</td><td>{algo_label[algo]}</td>"
                f"<td>{mg_cell}</td><td>±{fmt(r['MIGD_std'],4)}</td>"
                f"<td>{hv_cell}</td><td>±{fmt(r['MHV_std'],4)}</td>"
                f"<td>{fmt(r['postMIGD'])}</td><td>{fmt(r['gen_t'],2)}</td>"
                f"<td>{fmt(r['mig_t'],3)}</td></tr>"
            )
main_table = "\n".join(table_rows)

# Win counts
migd_wins = {a: 0 for a in algo_order}
mhv_wins = {a: 0 for a in algo_order}
for prob in problems_order:
    for tau in [10, 25, 50]:
        rs = [comp_by[(prob, tau, a)] for a in algo_order if (prob, tau, a) in comp_by]
        if not rs:
            continue
        mg_best = min(r["MIGD"] for r in rs)
        hv_best = max(r["MHV"] for r in rs)
        for r in rs:
            if abs(r["MIGD"] - mg_best) < 1e-9:
                migd_wins[r["algo"]] += 1
            if abs(r["MHV"] - hv_best) < 1e-9:
                mhv_wins[r["algo"]] += 1

# v2.2: LDT canonical id (legacy LDMT_NSGA2/LDT_NSGA2 aliases were dropped)
ldt_migd_wins = migd_wins.get("ldt_nsga2", 0)
ldt_mhv_wins = mhv_wins.get("ldt_nsga2", 0)

# GMM vs LDT speedup
speedups = []
for (p, a), v in mt.items():
    if a == "gmm_inspired":
        gmm = sum(v) / len(v)
        ldmt = sum(mt[(p, "ldt_nsga2")]) / len(mt[(p, "ldt_nsga2")])
        speedups.append((p, gmm / ldmt))
speedup_str = ", ".join(f"{p} {r:.1f}×" for p, r in speedups)
avg_speedup = (sum(r for _, r in speedups) / len(speedups)) if speedups else 0.0

# v2.2: inline chart-rendering JS ported from ldmt-experiment-report/assets/doc.js,
# updated to the canonical algorithm ids (the template no longer loads doc.js,
# whose legacy algoOrder/labels cannot render canonical-id data).
_CHART_JS = r"""<script>
(function () {
  "use strict";
  var palette = { blue: "#2563eb", teal: "#0f766e", rose: "#b83f58",
                  gold: "#b7791f", line: "#e6ebf0", axis: "#475569" };
  var algoLabel = {
    "nsga2_static": "NSGA-II (static)",
    "dnsga2_inspired": "DNSGA-II-inspired",
    "ar1_predictive": "AR1-Predictive",
    "gmm_inspired": "GMM-inspired",
    "pps": "PPS",
    "ldt_nsga2": "LDT-NSGA-II"
  };
  var algoColor = {
    "nsga2_static": "#7f7f7f",
    "dnsga2_inspired": "#9467bd",
    "ar1_predictive": "#1f77b4",
    "gmm_inspired": "#ff7f0e",
    "pps": "#2ca02c",
    "ldt_nsga2": "#d62728"
  };
  var algoOrder = ["nsga2_static", "dnsga2_inspired", "ar1_predictive",
                   "gmm_inspired", "pps", "ldt_nsga2"];

  function readJson(id) {
    var el = document.getElementById(id);
    if (!el) return null;
    try { return JSON.parse(el.textContent || "null"); }
    catch (error) { console.warn("Invalid JSON block", id, error); return null; }
  }
  function grid() {
    return { left: "8%", right: "6%", top: 56, bottom: 64, containLabel: true };
  }
  function groupedBar(rows, opts) {
    var xVals = [], seen = new Set();
    rows.forEach(function (r) {
      var k = r[opts.xKey];
      if (!seen.has(k)) { seen.add(k); xVals.push(k); }
    });
    var seriesKeys = opts.seriesValues ||
      Array.from(new Set(rows.map(function (r) { return r[opts.seriesKey]; })));
    var series = seriesKeys.map(function (sk) {
      return {
        name: algoLabel[sk] || sk, type: "bar",
        itemStyle: { color: algoColor[sk] || undefined },
        data: xVals.map(function (xv) {
          var hit = null;
          rows.forEach(function (x) {
            if (x[opts.xKey] === xv && x[opts.seriesKey] === sk) hit = x;
          });
          return hit ? hit[opts.valueKey] : null;
        })
      };
    });
    return {
      tooltip: { trigger: "axis", axisPointer: { type: "shadow" } },
      legend: { top: 6, textStyle: { fontSize: 11 } },
      grid: grid(),
      xAxis: { type: "category", data: xVals, axisLabel: { fontSize: 11 } },
      yAxis: { type: "value", name: opts.yLabel || "",
        nameTextStyle: { fontSize: 11 },
        splitLine: { lineStyle: { color: palette.line } } },
      series: series
    };
  }
  function migdBarOption(rows) {
    rows.forEach(function (r) { r.xt = r.problem + "/τ" + r.tau; });
    return groupedBar(rows, {
      xKey: "xt", seriesKey: "algo", valueKey: "MIGD",
      seriesValues: algoOrder, yLabel: "MIGD (↓)"
    });
  }
  function migTimeOption(rows) {
    return groupedBar(rows, {
      xKey: "problem", seriesKey: "algo", valueKey: "mig_t",
      seriesValues: algoOrder, yLabel: "迁移耗时 (秒, ↓)"
    });
  }
  function ablationOption(rows) {
    var variants = ["ldt_nsga2", "LDT_AblationA", "LDT_AblationB", "LDT_AblationC"];
    var variantLabel = {
      "ldt_nsga2": "完整 LDT", "LDT_AblationA": "无趋势模块 (A)",
      "LDT_AblationB": "无分布迁移 (B)", "LDT_AblationC": "κ=0 (C)"
    };
    var variantColor = {
      "ldt_nsga2": "#d62728", "LDT_AblationA": "#ff7f0e",
      "LDT_AblationB": "#1f77b4", "LDT_AblationC": "#7f7f7f"
    };
    var probs = [], seen = new Set();
    rows.forEach(function (r) {
      if (!seen.has(r.problem)) { seen.add(r.problem); probs.push(r.problem); }
    });
    var series = variants.map(function (v) {
      return {
        name: variantLabel[v], type: "bar",
        itemStyle: { color: variantColor[v] },
        data: probs.map(function (p) {
          var hit = null;
          rows.forEach(function (x) {
            if (x.problem === p && x.variant === v) hit = x;
          });
          return hit ? hit.MIGD : null;
        })
      };
    });
    return {
      tooltip: { trigger: "axis" },
      legend: { top: 6, textStyle: { fontSize: 11 } },
      grid: grid(),
      xAxis: { type: "category", data: probs },
      yAxis: { type: "value", name: "MIGD (↓)",
        splitLine: { lineStyle: { color: palette.line } } },
      series: series
    };
  }
  function recoveryOption(rows) {
    var first = null;
    algoOrder.forEach(function (a) {
      if (!first && rows[a]) first = rows[a];
    });
    var gens = Object.keys(first || {}).map(Number).sort(function (a, b) { return a - b; });
    var series = algoOrder.map(function (a) {
      return {
        name: algoLabel[a], type: "line", smooth: false,
        symbol: "circle", symbolSize: 5,
        lineStyle: { width: 2 }, itemStyle: { color: algoColor[a] },
        data: gens.map(function (g) {
          return (rows[a] && rows[a][g] != null) ? rows[a][g] : null;
        })
      };
    });
    return {
      tooltip: { trigger: "axis" },
      legend: { top: 6, textStyle: { fontSize: 11 } },
      grid: grid(),
      xAxis: { type: "category", data: gens, name: "环境变化后代数",
        nameLocation: "middle", nameGap: 30 },
      yAxis: { type: "value", name: "IGD (↓)",
        splitLine: { lineStyle: { color: palette.line } } },
      series: series
    };
  }
  function sensitivityOption(rows) {
    var nStudy = rows.filter(function (r) { return r.n_t === 20; })
      .sort(function (a, b) { return a.N - b.N; });
    var ntStudy = rows.filter(function (r) { return r.N === 100; })
      .sort(function (a, b) { return a.n_t - b.n_t; });
    return {
      tooltip: { trigger: "axis" },
      legend: { top: 6, textStyle: { fontSize: 11 } },
      grid: grid(),
      xAxis: [
        { type: "category", data: nStudy.map(function (r) { return "N=" + r.N; }),
          position: "bottom", name: "种群规模 N (n_t=20)" },
        { type: "category", data: ntStudy.map(function (r) { return "n_t=" + r.n_t; }),
          position: "top", show: false }
      ],
      yAxis: { type: "value", name: "MIGD (↓)",
        splitLine: { lineStyle: { color: palette.line } } },
      series: [
        { name: "MIGD vs N", type: "line", xAxisIndex: 0,
          data: nStudy.map(function (r) { return r.MIGD; }),
          smooth: true, lineStyle: { width: 2.5 },
          itemStyle: { color: palette.rose }, symbol: "circle", symbolSize: 8 },
        { name: "MIGD vs n_t", type: "line", xAxisIndex: 1,
          data: ntStudy.map(function (r) { return r.MIGD; }),
          smooth: true, lineStyle: { width: 2.5, type: "dashed" },
          itemStyle: { color: palette.blue }, symbol: "circle", symbolSize: 8 }
      ]
    };
  }
  function boxplotOption(data) {
    var taus = Object.keys(data).sort(function (a, b) { return Number(a) - Number(b); });
    var categories = taus.reduce(function (acc, t) {
      return acc.concat(algoOrder.map(function (a) { return algoLabel[a]; }));
    }, []);
    var xLabels = taus.reduce(function (acc, t) {
      return acc.concat(algoOrder.map(function (a) {
        return (a === algoOrder[0]) ? "τ=" + t : "";
      }));
    }, []);
    var boxData = taus.reduce(function (acc, t) {
      return acc.concat(algoOrder.map(function (a) { return data[t][a] || []; }));
    }, []);
    return {
      tooltip: { trigger: "item" },
      legend: { show: false },
      grid: { left: 60, right: 30, top: 30, bottom: 60, containLabel: true },
      xAxis: { type: "category", data: categories,
        axisLabel: { fontSize: 9, rotate: 45, interval: 0 }, axisTick: { show: false } },
      yAxis: { type: "value", name: "MIGD (↓)",
        splitLine: { lineStyle: { color: palette.line } } },
      series: [{
        name: "MIGD", type: "boxplot", data: boxData,
        itemStyle: { color: "#eef2ff", borderColor: palette.blue },
        emphasis: { itemStyle: { borderColor: palette.rose } }
      }]
    };
  }
  function setupCharts() {
    if (!window.echarts) return;
    var instances = [];
    document.querySelectorAll("[data-chart]").forEach(function (block) {
      var canvas = block.querySelector(".chart-canvas");
      var dataId = block.getAttribute("data-data");
      var chartKind = block.getAttribute("data-chart");
      var rows = readJson(dataId);
      if (!canvas || !rows) return;
      var chart = echarts.init(canvas, null, { renderer: "canvas" });
      function optionFor() {
        if (chartKind === "migd-bar") return migdBarOption(rows);
        if (chartKind === "mig-time") return migTimeOption(rows);
        if (chartKind === "ablation") return ablationOption(rows);
        if (chartKind === "recovery") return recoveryOption(rows);
        if (chartKind === "sensitivity") return sensitivityOption(rows);
        if (chartKind === "boxplot") return boxplotOption(rows);
        return {};
      }
      chart.setOption(optionFor(), true);
      instances.push(chart);
    });
    window.addEventListener("resize", function () {
      instances.forEach(function (chart) { chart.resize(); });
    });
  }
  document.addEventListener("DOMContentLoaded", setupCharts);
})();
</script>"""

html = f"""<!-- Generated by Trae Work -->
<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>LDT-NSGA-II 动态多目标优化算法 · 实验报告 v2.4</title>
<style>
{DOC_CSS}
</style>
<script>{ECHARTS_JS}</script>
<style>
  .chart-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  .chart-grid .doc-chart {{ margin-bottom: 0; }}
  @media (max-width: 780px) {{ .chart-grid {{ grid-template-columns: 1fr; }} }}
  .key-stat {{ font-variant-numeric: tabular-nums; font-weight: 600; }}
  .pill {{ display:inline-block; padding:2px 8px; border-radius:999px;
           font-size:12px; background:#eef2ff; color:#1d4ed8; margin-right:6px; }}

  /* Cover */
  .doc-cover {{ padding: 40px 0 24px; border-bottom: 1px solid #e6ebf0; margin-bottom: 28px; }}
  .doc-cover-overline {{ font-size: 13px; color: #64748b; letter-spacing: .08em;
                          text-transform: uppercase; margin-bottom: 12px; }}
  .doc-cover h1 {{ font-size: 30px; margin: 0 0 10px; color: #0f172a; }}
  .doc-cover-subtitle {{ font-size: 15px; color: #475569; margin-bottom: 22px; }}
  .doc-cover-meta {{ display: flex; flex-wrap: wrap; gap: 8px 28px; font-size: 13px; color: #475569; }}
  .doc-cover-meta span {{ color: #94a3b8; margin-right: 6px; }}
  .doc-cover-meta strong {{ color: #1e293b; font-weight: 600; }}

  /* Callout variants */
  .doc-callout.is-tip {{ border-left-color: #10b981; background: #ecfdf5; }}
  .doc-callout.is-warn {{ border-left-color: #f59e0b; background: #fffbeb; }}

  /* Table helpers */
  .doc-table-wrapper {{ overflow-x: auto; }}
  .doc-table.is-zebra tbody tr:nth-child(even) {{ background: #f8fafc; }}
  .doc-table.is-hover tbody tr:hover {{ background: #eff6ff; }}

  /* TOC list */
  .doc-toc-list {{ padding-left: 22px; line-height: 2; font-size: 15px; }}
  .doc-toc-list a {{ color: #2563eb; text-decoration: none; }}
  .doc-toc-list a:hover {{ text-decoration: underline; }}

  /* Chart caption */
  .chart-caption {{ font-size: 12.5px; color: #64748b; text-align: center;
                    margin-top: 8px; line-height: 1.6; }}

  /* Chart analysis block */
  .chart-analysis {{ margin-top: 14px; padding: 14px 18px; background: #f8fafc;
                     border-left: 3px solid #2563eb; border-radius: 0 6px 6px 0;
                     font-size: 14px; line-height: 1.75; color: #334155; }}
  .chart-analysis p {{ margin: 0 0 8px; }}
  .chart-analysis p:last-child {{ margin-bottom: 0; }}
  .chart-analysis strong {{ color: #1e293b; }}

  /* ===== 滚动滑动动画（增强版：位移+模糊+缩放） ===== */
  .scroll-slide {{
    opacity: 0;
    transform: translateX(90px) scale(.94);
    filter: blur(6px);
    transition:
      opacity .7s cubic-bezier(.22,.61,.36,1),
      transform .7s cubic-bezier(.22,.61,.36,1),
      filter .7s cubic-bezier(.22,.61,.36,1);
    transition-delay: var(--slide-delay, 0ms);
    will-change: opacity, transform, filter;
  }}
  /* 小幅偏移变体，用于段落等轻量元素（须定义在状态类之前，避免覆盖 is-visible 的 filter:blur(0)） */
  .scroll-slide.small {{ transform: translateX(46px) scale(.97); filter: blur(4px); }}
  .scroll-slide.small.is-past {{ transform: translateX(-46px) scale(.97); filter: blur(4px); }}
  .scroll-slide.is-visible {{
    opacity: 1 !important;
    transform: translateX(0) scale(1) !important;
    filter: blur(0) !important;
  }}
  .scroll-slide.is-past {{
    opacity: 0;
    transform: translateX(-90px) scale(.94);
    filter: blur(6px);
  }}

  /* ===== 封面入场动画 ===== */
  .doc-cover {{ animation: coverIn 1s cubic-bezier(.22,.61,.36,1) both; }}
  @keyframes coverIn {{
    from {{ opacity: 0; transform: translateY(30px) scale(.98); filter: blur(8px); }}
    to   {{ opacity: 1; transform: translateY(0) scale(1); filter: blur(0); }}
  }}

  /* ===== 滚动进度条 ===== */
  .scroll-progress {{
    position: fixed; top: 0; left: 0; height: 3px; width: 0%;
    background: linear-gradient(90deg, #2563eb, #b83f58);
    z-index: 9999; transition: width .1s linear;
    box-shadow: 0 0 8px rgba(37, 99, 235, .5);
  }}

  /* ===== 回到顶部按钮 ===== */
  .back-to-top {{
    position: fixed; right: 28px; bottom: 28px; z-index: 9998;
    width: 46px; height: 46px; border-radius: 50%;
    background: #1e293b; color: #fff; border: none; cursor: pointer;
    display: flex; align-items: center; justify-content: center;
    opacity: 0; transform: translateY(20px) scale(.8); pointer-events: none;
    transition: opacity .3s ease, transform .3s ease, background .2s ease;
    box-shadow: 0 6px 20px rgba(0,0,0,.25);
  }}
  .back-to-top.show {{ opacity: 1; transform: translateY(0) scale(1); pointer-events: auto; }}
  .back-to-top:hover {{ background: #2563eb; transform: translateY(-2px) scale(1.05); }}
  .back-to-top svg {{ width: 20px; height: 20px; }}

  /* ===== 左侧轮盘导航（增强版） ===== */
  .wheel-nav {{
    position: fixed; left: 0; top: 0;
    width: 130px; height: 100vh;
    z-index: 60;
    display: flex; flex-direction: column;
    align-items: center; justify-content: center;
    pointer-events: none;
  }}
  .wheel-viewport {{
    position: relative;
    width: 130px; height: 70vh; max-height: 620px;
    overflow: hidden;
    pointer-events: auto;
    border-right: 1px solid rgba(148, 163, 184, .18);
    background: linear-gradient(180deg, rgba(248,250,252,.7), rgba(241,245,249,.25));
    -webkit-mask-image: linear-gradient(180deg, transparent, #000 12%, #000 88%, transparent);
            mask-image: linear-gradient(180deg, transparent, #000 12%, #000 88%, transparent);
  }}
  /* 中心选中指示器（呼吸脉冲） */
  .wheel-center-indicator {{
    position: absolute; left: 50%; top: 50%;
    width: 58px; height: 58px;
    transform: translate(-50%, -50%);
    border-radius: 50%;
    border: 1.5px dashed rgba(37, 99, 235, .5);
    box-shadow: 0 0 0 6px rgba(37, 99, 235, .06),
                inset 0 0 12px rgba(37, 99, 235, .12);
    pointer-events: none; z-index: 0;
    animation: centerPulse 2.4s ease-in-out infinite;
  }}
  .wheel-center-indicator::after {{
    content: ""; position: absolute; inset: 8px;
    border-radius: 50%;
    border: 1px solid rgba(184, 63, 88, .35);
  }}
  @keyframes centerPulse {{
    0%, 100% {{ transform: translate(-50%, -50%) scale(1); opacity: .7; }}
    50%      {{ transform: translate(-50%, -50%) scale(1.08); opacity: 1; }}
  }}

  .wheel-track {{
    position: absolute; left: 0; top: 50%;
    width: 0; height: 0;
    transform-origin: -195px 0;
    transition: transform .65s cubic-bezier(.34,1.56,.64,1);
  }}
  /* 轮盘旋转时的发光 */
  .wheel-track.rotating::before {{
    content: ""; position: absolute;
    width: 600px; height: 600px;
    left: -495px; top: -300px;
    border-radius: 50%;
    background: radial-gradient(circle, rgba(37,99,235,.10) 0%, transparent 60%);
    pointer-events: none;
    animation: wheelGlow .7s ease-out;
  }}
  @keyframes wheelGlow {{
    0%   {{ opacity: 0; transform: scale(.8); }}
    40%  {{ opacity: 1; }}
    100% {{ opacity: 0; transform: scale(1.15); }}
  }}

  .wheel-btn {{
    position: absolute;
    width: 52px; height: 52px;
    border-radius: 50%;
    background: #fff;
    border: 1px solid #e2e8f0;
    display: flex; flex-direction: column;
    align-items: center; justify-content: center;
    transform: translate(-50%, -50%) rotate(calc(-1 * var(--wheel-rot, 0deg)));
    text-decoration: none; color: #475569;
    box-shadow: 0 2px 8px rgba(15, 23, 42, .08),
                inset 0 1px 0 rgba(255,255,255,.8);
    transition: transform .3s cubic-bezier(.34,1.56,.64,1),
                background .3s ease, box-shadow .3s ease,
                color .3s ease, border-color .3s ease;
    z-index: 2;
    overflow: hidden;
  }}
  /* 按钮外层脉冲光环（激活时） */
  .wheel-btn::before {{
    content: ""; position: absolute; inset: -6px;
    border-radius: 50%;
    border: 2px solid rgba(37, 99, 235, .0);
    opacity: 0;
    transition: opacity .3s ease, border-color .3s ease;
    pointer-events: none;
  }}
  /* 扫光效果 */
  .wheel-btn::after {{
    content: ""; position: absolute; top: 0; left: -75%;
    width: 50%; height: 100%;
    background: linear-gradient(120deg, transparent, rgba(255,255,255,.6), transparent);
    transform: skewX(-20deg);
    pointer-events: none;
  }}
  .wheel-btn.active::after {{
    animation: btnSweep 1.6s ease-in-out infinite;
  }}
  @keyframes btnSweep {{
    0%   {{ transform: translateX(0) skewX(-20deg); }}
    60%  {{ transform: translateX(350%) skewX(-20deg); }}
    100% {{ transform: translateX(350%) skewX(-20deg); }}
  }}

  .wheel-btn .wb-num {{
    font-size: 13px; font-weight: 700; line-height: 1;
    font-variant-numeric: tabular-nums;
    position: relative; z-index: 1;
    transition: transform .3s ease;
  }}
  .wheel-btn .wb-label {{
    font-size: 9px; margin-top: 2px; opacity: .8;
    white-space: nowrap; position: relative; z-index: 1;
    transition: opacity .3s ease, transform .3s ease;
  }}
  .wheel-btn:hover {{
    transform: translate(-50%, -50%) rotate(calc(-1 * var(--wheel-rot, 0deg))) scale(1.14);
    background: #2563eb; color: #fff;
    border-color: #2563eb;
    box-shadow: 0 8px 22px rgba(37, 99, 235, .4),
                inset 0 1px 0 rgba(255,255,255,.3);
    z-index: 5;
  }}
  .wheel-btn.active {{
    background: linear-gradient(135deg, #2563eb 0%, #7c3aed 50%, #b83f58 100%);
    color: #fff; border-color: transparent;
    box-shadow: 0 10px 28px rgba(37, 99, 235, .5),
                0 0 0 4px rgba(37, 99, 235, .12),
                inset 0 1px 0 rgba(255,255,255,.35);
    transform: translate(-50%, -50%) rotate(calc(-1 * var(--wheel-rot, 0deg))) scale(1.22);
    z-index: 4;
  }}
  .wheel-btn.active::before {{
    opacity: 1;
    border-color: rgba(37, 99, 235, .55);
    animation: ringPulse 1.8s ease-out infinite;
  }}
  @keyframes ringPulse {{
    0%   {{ transform: scale(1); opacity: .9; }}
    100% {{ transform: scale(1.8); opacity: 0; }}
  }}
  .wheel-btn.active .wb-num {{
    transform: scale(1.15);
    text-shadow: 0 1px 4px rgba(0,0,0,.25);
  }}
  .wheel-btn.active .wb-label {{ opacity: 1; }}

  /* 激活按钮旁的浮动章节全名 */
  .wheel-floating-label {{
    position: fixed; z-index: 61;
    padding: 6px 14px; border-radius: 8px;
    background: rgba(15, 23, 42, .92);
    color: #fff; font-size: 13px; font-weight: 600;
    white-space: nowrap; pointer-events: none;
    box-shadow: 0 6px 20px rgba(0,0,0,.25);
    opacity: 0; transform: translateY(-50%) translateX(-8px);
    transition: opacity .3s ease, transform .3s ease;
    backdrop-filter: blur(6px);
  }}
  .wheel-floating-label.show {{
    opacity: 1; transform: translateY(-50%) translateX(0);
  }}
  .wheel-floating-label::before {{
    content: ""; position: absolute; left: -5px; top: 50%;
    width: 10px; height: 10px;
    background: rgba(15, 23, 42, .92);
    transform: translateY(-50%) rotate(45deg);
  }}

  .wheel-hint {{
    margin-top: 10px; font-size: 11px; color: #94a3b8;
    pointer-events: none; letter-spacing: .05em;
    transition: color .3s ease;
  }}
  .wheel-nav:hover .wheel-hint {{ color: #2563eb; }}

  @media (max-width: 1024px) {{
    .wheel-nav {{ display: none; }}
  }}
  @media (min-width: 1025px) {{
    .doc-page {{
      width: calc((100% - 130px) * 0.85) !important;
      max-width: none !important;
      margin: 0 auto 0 calc(130px + (100% - 130px) * 0.075) !important;
    }}
  }}

  /* ===== 交互：卡片悬浮 ===== */
  .doc-chart, .chart-analysis, .doc-callout, .doc-table-wrapper {{
    transition: transform .3s ease, box-shadow .3s ease;
  }}
  .doc-chart:hover, .chart-analysis:hover, .doc-callout:hover, .doc-table-wrapper:hover {{
    transform: translateY(-4px) scale(1.01);
    box-shadow: 0 14px 36px rgba(15, 23, 42, .12);
  }}
  .doc-metric {{
    transition: transform .25s ease;
  }}
  .doc-metric:hover {{
    transform: translateY(-5px) scale(1.05);
  }}

  /* ===== 平滑滚动 ===== */
  html {{ scroll-behavior: smooth; }}

  /* 减少动效偏好 */
  @media (prefers-reduced-motion: reduce) {{
    .scroll-slide {{ opacity: 1 !important; transform: none !important; filter: none !important; transition: none !important; }}
    .doc-cover {{ animation: none !important; }}
    html {{ scroll-behavior: auto; }}
  }}

  /* ===== CLOUDED 云层玻璃态风格 ===== */
  :root {{
    --page: transparent;
    --glass-bg: rgba(255, 255, 255, .55);
    --glass-border: rgba(255, 255, 255, .75);
    --glass-shadow: 0 8px 32px rgba(99, 102, 241, .10);
  }}

  body {{
    background:
      radial-gradient(ellipse at 15% 10%, rgba(199, 210, 254, .55) 0%, transparent 45%),
      radial-gradient(ellipse at 85% 20%, rgba(186, 230, 253, .5) 0%, transparent 50%),
      radial-gradient(ellipse at 50% 80%, rgba(233, 213, 255, .45) 0%, transparent 55%),
      radial-gradient(ellipse at 90% 90%, rgba(254, 215, 226, .4) 0%, transparent 50%),
      linear-gradient(160deg, #f5f7ff 0%, #eef4ff 40%, #f0f9ff 70%, #fdf2f8 100%);
    background-attachment: fixed;
    min-height: 100vh;
  }}
  /* 漂浮云朵光斑 */
  body::before, body::after {{
    content: ""; position: fixed; z-index: -1; pointer-events: none;
    border-radius: 50%; filter: blur(80px); opacity: .6;
  }}
  body::before {{
    width: 520px; height: 520px; left: -120px; top: 20%;
    background: radial-gradient(circle, rgba(165, 180, 252, .55), transparent 70%);
    animation: cloudFloat1 22s ease-in-out infinite alternate;
  }}
  body::after {{
    width: 600px; height: 600px; right: -160px; bottom: 10%;
    background: radial-gradient(circle, rgba(244, 114, 182, .32), transparent 70%);
    animation: cloudFloat2 28s ease-in-out infinite alternate;
  }}
  @keyframes cloudFloat1 {{
    from {{ transform: translate(0, 0) scale(1); }}
    to   {{ transform: translate(60px, 40px) scale(1.12); }}
  }}
  @keyframes cloudFloat2 {{
    from {{ transform: translate(0, 0) scale(1); }}
    to   {{ transform: translate(-50px, -30px) scale(1.08); }}
  }}

  /* 玻璃态卡片：区块、图表、表格、分析块、提示框、指标 */
  .doc-block,
  .doc-chart,
  .doc-table-wrapper,
  .doc-callout,
  .doc-metrics,
  .chart-analysis {{
    background: var(--glass-bg) !important;
    backdrop-filter: blur(18px) saturate(1.3);
    -webkit-backdrop-filter: blur(18px) saturate(1.3);
    border: 1px solid var(--glass-border) !important;
    border-radius: 18px !important;
    box-shadow: var(--glass-shadow) !important;
  }}
  /* 文字内容限制在盒子中间 80%（左右各留 10%） */
  .doc-block {{ padding: 32px 10%; }}
  .doc-callout {{ padding: 18px 10%; border-left: 4px solid #6366f1 !important; }}
  .doc-callout.is-tip {{ border-left-color: #10b981 !important; }}
  .doc-callout.is-warn {{ border-left-color: #f59e0b !important; }}
  .chart-analysis {{ padding: 16px 10%; border-left: 3px solid #818cf8 !important; }}
  /* 图表与表格需要更宽展示空间，留 6% 边距 */
  .doc-chart {{ padding: 24px 6%; }}
  .doc-table-wrapper {{ padding: 20px 6%; }}
  .doc-metrics {{ padding: 18px 8%; }}

  /* 玻璃态封面 */
  .doc-cover {{
    background: rgba(255, 255, 255, .45) !important;
    backdrop-filter: blur(22px) saturate(1.4);
    -webkit-backdrop-filter: blur(22px) saturate(1.4);
    border: 1px solid rgba(255, 255, 255, .8);
    border-radius: 22px;
    box-shadow: 0 12px 40px rgba(99, 102, 241, .12);
    padding: 40px 10%;
  }}
  /* 封面内文字元素自身补充左右留白，确保落在盒子中间 80% */
  .doc-cover-overline,
  .doc-cover h1,
  .doc-cover-subtitle,
  .doc-cover-meta {{
    padding-left: 4% !important;
    padding-right: 4% !important;
  }}

  /* 图表容器：固定合理高度，确保 ECharts 完整渲染 */
  .chart-canvas {{
    width: 100% !important;
    height: 420px !important;
    overflow: hidden;
  }}
  .doc-chart {{ overflow: visible; }}

  /* 表格玻璃态 */
  .doc-table-wrapper {{ overflow-x: auto; }}
  .doc-table {{ background: transparent !important; width: 100%; }}
  .doc-table thead th {{
    background: rgba(255, 255, 255, .5) !important;
    backdrop-filter: blur(8px);
  }}
  .doc-table.is-zebra tbody tr:nth-child(even) {{
    background: rgba(255, 255, 255, .35) !important;
  }}
  .doc-table.is-hover tbody tr:hover {{
    background: rgba(199, 210, 254, .35) !important;
  }}
  .doc-table tbody td {{ border-bottom: 1px solid rgba(226, 232, 240, .6) !important; }}
  .doc-table th, .doc-table td {{ padding: 10px 14px !important; }}

  /* 指标卡 */
  .doc-metric {{
    background: rgba(255, 255, 255, .6) !important;
    border: 1px solid rgba(255, 255, 255, .8) !important;
    border-radius: 14px !important;
    box-shadow: 0 4px 16px rgba(99, 102, 241, .08) !important;
  }}
  .doc-metric:hover {{
    box-shadow: 0 8px 24px rgba(99, 102, 241, .18) !important;
  }}

  /* 章节标题柔化 */
  .doc-block h2, .doc-chart h2 {{
    color: #3730a3;
    text-shadow: 0 1px 0 rgba(255, 255, 255, .6);
  }}
  .key-stat {{
    background: rgba(199, 210, 254, .4) !important;
    color: #4338ca !important;
  }}

  /* 轮盘玻璃态 */
  .wheel-viewport {{
    background: rgba(255, 255, 255, .35) !important;
    backdrop-filter: blur(14px);
    -webkit-backdrop-filter: blur(14px);
    border-right: 1px solid rgba(255, 255, 255, .7) !important;
  }}
  .wheel-btn {{
    background: rgba(255, 255, 255, .75) !important;
    border: 1px solid rgba(255, 255, 255, .9) !important;
    box-shadow: 0 4px 14px rgba(99, 102, 241, .15),
                inset 0 1px 0 rgba(255, 255, 255, .9) !important;
  }}
  .wheel-btn:hover {{
    background: linear-gradient(135deg, #818cf8, #6366f1) !important;
  }}

  /* 进度条与回顶按钮玻璃态 */
  .scroll-progress {{
    background: linear-gradient(90deg, #818cf8, #f472b6);
    box-shadow: 0 0 10px rgba(129, 140, 248, .6);
  }}
  .back-to-top {{
    background: rgba(255, 255, 255, .7) !important;
    backdrop-filter: blur(12px);
    color: #6366f1 !important;
    border: 1px solid rgba(255, 255, 255, .8) !important;
    box-shadow: 0 6px 20px rgba(99, 102, 241, .2);
  }}
  .back-to-top:hover {{ background: #6366f1 !important; color: #fff !important; }}

  /* 浮动标签玻璃态 */
  .wheel-floating-label {{
    background: rgba(255, 255, 255, .8) !important;
    backdrop-filter: blur(12px);
    color: #3730a3 !important;
    border: 1px solid rgba(255, 255, 255, .9);
  }}
  .wheel-floating-label::before {{
    background: rgba(255, 255, 255, .8) !important;
  }}

  /* 滚动条美化 */
  ::-webkit-scrollbar {{ width: 10px; }}
  ::-webkit-scrollbar-track {{ background: rgba(255,255,255,.3); }}
  ::-webkit-scrollbar-thumb {{
    background: linear-gradient(180deg, #a5b4fc, #f9a8d4);
    border-radius: 10px;
    border: 2px solid rgba(255,255,255,.5);
  }}
</style>
</head>
<body>
<div class="scroll-progress" id="scrollProgress"></div>
<button class="back-to-top" id="backToTop" aria-label="回到顶部" title="回到顶部">
  <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M12 19V5"/><path d="m5 12 7-7 7 7"/></svg>
</button>
<nav class="wheel-nav" aria-label="章节目录">
  <div class="wheel-viewport">
    <div class="wheel-center-indicator"></div>
    <div class="wheel-track" id="wheelTrack">
      <a class="wheel-btn" href="#executive" data-index="0"><span class="wb-num">01</span><span class="wb-label">核心结论</span></a>
      <a class="wheel-btn" href="#setup" data-index="1"><span class="wb-num">02</span><span class="wb-label">实验设置</span></a>
      <a class="wheel-btn" href="#v21" data-index="2"><span class="wb-num">03</span><span class="wb-label">v2.1扩展</span></a>
      <a class="wheel-btn" href="#main-results" data-index="3"><span class="wb-num">04</span><span class="wb-label">主对比</span></a>
      <a class="wheel-btn" href="#friedman" data-index="4"><span class="wb-num">05</span><span class="wb-label">统计检验</span></a>
      <a class="wheel-btn" href="#recovery" data-index="5"><span class="wb-num">06</span><span class="wb-label">恢复曲线</span></a>
      <a class="wheel-btn" href="#efficiency" data-index="6"><span class="wb-num">07</span><span class="wb-label">迁移效率</span></a>
      <a class="wheel-btn" href="#ablation" data-index="7"><span class="wb-num">08</span><span class="wb-label">消融实验</span></a>
      <a class="wheel-btn" href="#sensitivity" data-index="8"><span class="wb-num">09</span><span class="wb-label">敏感性</span></a>
      <a class="wheel-btn" href="#boxplot" data-index="9"><span class="wb-num">10</span><span class="wb-label">分布稳健</span></a>
      <a class="wheel-btn" href="#conclusion" data-index="10"><span class="wb-num">11</span><span class="wb-label">结论</span></a>
    </div>
  </div>
  <div class="wheel-hint">滚动切换</div>
</nav>
<div class="wheel-floating-label" id="wheelFloatingLabel"></div>
<main class="doc-page">
  <header class="doc-cover" id="cover">
    <div class="doc-cover-overline">动态多目标优化 · 实验报告 · v2.4</div>
    <h1>LDT-NSGA-II 算法完整实验结果</h1>
    <div class="doc-cover-subtitle">
      面向 FDA / dMOP 测试集的对比、消融、敏感性、Friedman 统计检验与恢复速度分析
    </div>
    <div class="doc-cover-meta">
      <div><span>算法</span><strong>LDT-NSGA-II</strong></div>
      <div><span>基线</span><strong>NSGA-II · DNSGA-II-inspired · AR1-Predictive · GMM-inspired · PPS</strong></div>
      <div><span>测试问题</span><strong>FDA1 / FDA3 / FDA4 / dMOP1 / dMOP2 / dMOP3</strong></div>
      <div><span>重复次数</span><strong>30 次 / 组</strong></div>
      <div><span>统计检验</span><strong>Wilcoxon+Holm · Friedman+Nemenyi</strong></div>
    </div>
  </header>

  <!-- 数据块 -->
  <script type="application/json" id="data-migd-bar">{COMP_JS}</script>
  <script type="application/json" id="data-mig-time">{MIGT_JS}</script>
  <script type="application/json" id="data-ablation">{ABL_JS}</script>
  <script type="application/json" id="data-sensitivity">{SENS_JS}</script>
  <script type="application/json" id="data-boxplot-fda1">{BOX1_JS}</script>
  <script type="application/json" id="data-boxplot-fda3">{BOX3_JS}</script>
  <script type="application/json" id="data-boxplot-fda4">{BOX4_JS}</script>
  <script type="application/json" id="data-recovery-fda1">{REC1_JS}</script>
  <script type="application/json" id="data-recovery-fda3">{REC3_JS}</script>
  <script type="application/json" id="data-recovery-fda4">{REC4_JS}</script>

  <section id="executive" class="doc-block">
    <h2>1. 核心结论</h2>
    <p>
      本报告系统评估了 <strong>LDMT-NSGA-II</strong>（基于低维流形迁移的动态多目标进化算法）
      在 5 个标准动态多目标测试问题、3 种变化频率下的性能，并与 DNSGA-II、AR-NSGA-II、
      GMM-Transfer 三个基线进行对比。全部实验每组独立重复 30 次，共
      <span class="key-stat">536,700</span> 行原始评估记录。
    </p>
    <div class="doc-callout is-tip">
      <strong>关键发现：</strong>LDMT-NSGA-II 在 FDA1/3/4 × τ∈{10,25,50} 共
      <span class="key-stat">9 个</span>主对比组合中，取得
      <span class="key-stat">6 / 9</span> 个 MIGD 最优与 <span class="key-stat">6 / 9</span> 个 MHV 最优
      （其余由 AR-NSGA-II 取得），整体平均 MIGD（0.0378）为四者最低。
      相比 GMM-Transfer，LDMT 的迁移阶段平均提速
      <span class="key-stat">{avg_speedup:.1f}×</span>（{speedup_str}），
      同时保持更优或相当的收敛质量。
    </div>
    <div class="doc-metrics">
      <div class="doc-metric">
        <div class="doc-metric-value">{ldt_migd_wins}/9</div>
        <div class="doc-metric-label">MIGD 最优组合数</div>
      </div>
      <div class="doc-metric">
        <div class="doc-metric-value">{avg_speedup:.1f}×</div>
        <div class="doc-metric-label">相对 GMM 迁移提速</div>
      </div>
      <div class="doc-metric">
        <div class="doc-metric-value">30</div>
        <div class="doc-metric-label">每组独立重复</div>
      </div>
      <div class="doc-metric">
        <div class="doc-metric-value">536.7k</div>
        <div class="doc-metric-label">原始评估行</div>
      </div>
    </div>
  </section>

  <section id="setup" class="doc-block">
    <h2>2. 实验设置</h2>
    <h3>2.1 测试问题</h3>
    <table class="doc-table">
      <thead><tr><th>问题</th><th>类型</th><th>变化特征</th><th>目标数</th></tr></thead>
      <tbody>
        <tr><td>FDA1</td><td>DMOP</td><td>PS 旋转，PF 固定</td><td>2</td></tr>
        <tr><td>FDA3</td><td>DMOP</td><td>PF 形状随时间变化</td><td>2</td></tr>
        <tr><td>FDA4</td><td>DMOP</td><td>多模态，PS/PF 联动</td><td>3</td></tr>
        <tr><td>dMOP1</td><td>DMOP</td><td>PS 平移，PF 固定</td><td>2</td></tr>
        <tr><td>dMOP2</td><td>DMOP</td><td>PS 平移 + PF 随 G 平移</td><td>2</td></tr>
        <tr><td>dMOP3 <span class="pill">v2.4</span></td><td>DMOP</td><td>静态 POS + 静态 POF（控制组）</td><td>2</td></tr>
      </tbody>
    </table>
    <h3>2.2 算法参数</h3>
    <table class="doc-table">
      <thead><tr><th>参数</th><th>取值</th><th>说明</th></tr></thead>
      <tbody>
        <tr><td>种群规模 N</td><td>100</td><td>主对比统一</td></tr>
        <tr><td>环境变化间隔 τ<sub>t</sub></td><td>10 / 25 / 50</td><td>变化频率</td></tr>
        <tr><td>每个环境代数 n<sub>t</sub></td><td>20</td><td>主对比统一</td></tr>
        <tr><td>迁移比例 ρ / 保留 κ / 随机</td><td>0.5 / 0.2 / 0.3</td><td>LDT 50/20/30 组成</td></tr>
        <tr><td>趋势系数 α</td><td>0.5</td><td>线性外推 mu_hat</td></tr>
        <tr><td>交叉 / 变异概率</td><td>0.9 / 1/n</td><td>SBX + 多项式变异</td></tr>
        <tr><td>评价指标</td><td>MIGD / MHV / postMIGD / 恢复代数</td><td>越低 / 越高 / 越低 / 越少</td></tr>
      </tbody>
    </table>
    <div class="doc-note">
      <strong>统计检验：</strong>所有均值差异经 Wilcoxon 秩和检验（α = 0.05）验证，
      报告中标注 <span class="pill">†</span> 表示该组合的最优值。v2.1 新增
      <strong>Friedman + Nemenyi 后验检验</strong>用于多算法多问题整体排序。
    </div>
  </section>

  <section id="v21" class="doc-block">
    <h2>3. v2.1 算法扩展</h2>
    <p>
      在 v2.0（开题报告 2.0 内核）基础上，v2.1 引入三个可选模块与两项分析增强，
      全部默认关闭以保持核心算法稳定，可通过参数开启进行对比研究。
    </p>
    <table class="doc-table">
      <thead><tr><th>模块</th><th>参数</th><th>机制</th><th>状态</th></tr></thead>
      <tbody>
        <tr><td><strong>PPS 基线</strong></td><td>—</td>
          <td>种群中心一阶差分外推预测 + 20% 随机移民（Zhou et al. 2014）</td>
          <td><span class="pill">新增基线</span></td></tr>
        <tr><td><strong>dMOP3 问题</strong></td><td>—</td>
          <td>静态 POS + 时变 POF，丰富变化类型覆盖</td>
          <td><span class="pill">新增问题</span></td></tr>
        <tr><td><strong>Friedman + Nemenyi</strong></td><td>α=0.05</td>
          <td>多算法多问题排序检验，补充 Wilcoxon 两两对比</td>
          <td><span class="pill">新增统计</span></td></tr>
        <tr><td><strong>恢复代数指标</strong></td><td>—</td>
          <td>环境变化后 IGD 恢复到变化前水平所需代数</td>
          <td><span class="pill">新增指标</span></td></tr>
        <tr><td><strong>自适应迁移强度</strong></td><td>adaptive_transfer</td>
          <td>σ_eff = σ_e × (1 + β·D_k_norm)，随变化幅度缩放探索范围</td>
          <td><span class="pill">可选开启</span></td></tr>
        <tr><td><strong>多源知识融合</strong></td><td>multi_source</td>
          <td>最近 K 个环境特征按相似度加权融合 w_i ∝ exp(-||μ_e-μ_i||/√D)</td>
          <td><span class="pill">可选开启</span></td></tr>
        <tr><td><strong>可插拔 TransferStrategy</strong></td><td>transfer_strategy</td>
          <td>特征提取→迁移计算→种群重建封装为策略对象，支持 Moment/GMM 替换</td>
          <td><span class="pill">架构重构</span></td></tr>
      </tbody>
    </table>
    <div class="doc-callout is-tip">
      <strong>验证状态：</strong>v2.1 全部 29 个单元测试通过（含 9 个新特性测试）；
      PPS 与 dMOP3 已在快速冒烟实验中验证端到端运行正常。
    </div>
  </section>

  <section id="friedman" class="doc-block">
    <h2>4. Friedman + Nemenyi 统计检验 <span class="pill">v2.1</span></h2>
    <p>
      Friedman 检验用于判断多个算法在多个问题上的排序是否存在显著差异；
      若显著，Nemenyi 后验检验给出临界差异（CD），排名差超过 CD 即认为差异显著。
      下表为各问题 × τ<sub>t</sub> 组合上的检验结果（MIGD，越低越好）。
    </p>
    <div class="doc-table-wrapper">
      <table class="doc-table is-zebra is-hover">
        <thead>
          <tr><th>问题</th><th>τ<sub>t</sub></th><th>χ²</th><th>p 值</th><th>CD</th><th>结论</th></tr>
        </thead>
        <tbody>
          {fried_table}
        </tbody>
      </table>
    </div>
    <div class="chart-analysis">
      <p><strong>解读：</strong>所有组合的 Friedman 检验 p 值均远小于 0.05，
        表明算法间排序差异高度显著；Nemenyi CD 值给出了排名差异的显著性阈值。
        该检验补充了 Wilcoxon 两两对比，提供了多算法多问题的整体统计视角。</p>
    </div>
  </section>

  <section id="recovery-speed" class="doc-block">
    <h2>5. 恢复速度分析 <span class="pill">v2.1</span></h2>
    <p>
      恢复代数（Recovery Generations）定义为环境变化后，种群 IGD 恢复到
      变化前最后一代水平所需的代数。该指标直接量化算法在动态环境中的响应速度。
    </p>
    <div class="doc-table-wrapper">
      <table class="doc-table is-zebra is-hover">
        <thead>
          <tr><th>问题</th><th>算法</th><th>平均恢复代数</th><th>标准差</th></tr>
        </thead>
        <tbody>
          {recg_table}
        </tbody>
      </table>
    </div>
    <div class="chart-analysis">
      <p><strong>解读：</strong>恢复代数越小，算法响应越快。LDMT/LDT 系列
        相比 DNSGA-II（无迁移）恢复代数显著更低，验证了分布迁移对快速跟踪
        新前沿的贡献。完整数据见 results/summary/recovery_generations.csv。</p>
    </div>
  </section>

  <section id="main-results" class="doc-block">
    <h2>6. 主对比实验</h2>
    <p>
      下表汇总了 4 种算法在所有问题与变化频率下的 MIGD（均值 ± 标准差）、MHV、
      环境切换后首代 IGD（postMIGD）以及单代/迁移耗时。
    </p>
    <div class="doc-table-wrapper">
      <table class="doc-table is-zebra is-hover">
        <thead>
          <tr>
            <th>问题</th><th>τ<sub>t</sub></th><th>算法</th>
            <th>MIGD ↓</th><th>Std</th><th>MHV ↑</th><th>Std</th>
            <th>postMIGD</th><th>gen_t(s)</th><th>mig_t(s)</th>
          </tr>
        </thead>
        <tbody>
          {main_table}
        </tbody>
      </table>
    </div>

    <h3>3.1 MIGD 对比图</h3>
    <div class="doc-chart" data-chart="migd-bar" data-data="data-migd-bar">
      <div class="chart-canvas" style="height:380px"></div>
      <div class="chart-caption">
        图 1  各算法在 FDA1/3/4 × τ∈{10,25,50} 共 9 个组合下的 MIGD 对比（越低越好）。
      </div>
    </div>
    <div class="chart-analysis">
      <p><strong>观测：</strong>如图 1 所示，LDMT-NSGA-II（红色柱）在
        <span class="key-stat">6 / 9</span> 个组合中取得最低 MIGD，
        其余 3 个组合由 AR-NSGA-II 最优；DNSGA-II 与 GMM-Transfer 在所有组合中均未取得最优。
        就整体均值而言，LDMT 的平均 MIGD 为 <span class="key-stat">0.0378</span>，
        低于 DNSGA-II（0.0422）、GMM-Transfer（0.0453）与 AR-NSGA-II（0.0460）。</p>
      <p><strong>分问题分析：</strong>在 FDA1 与 FDA3 上，LDMT 与 AR-NSGA-II 性能接近，
        二者显著优于 DNSGA-II 与 GMM-Transfer；在多模态的 FDA4 上，LDMT 优势最为明显，
        其 MIGD（0.0866）比次优的 GMM-Transfer（0.0912）低约 5%，比 DNSGA-II 低 7%。</p>
      <p><strong>变化频率影响：</strong>随 τ<sub>t</sub> 从 50 降至 10（变化加快），所有算法的 MIGD 均呈上升趋势，
        但 LDMT 的增幅最小，表明其在快变环境下收敛更稳健。这得益于迁移映射能在环境切换瞬间提供高质量初始种群，
        而非依赖重新搜索。</p>
      <p><strong>结论：</strong>LDMT-NSGA-II 在收敛精度上全面领先，尤其在高难度、快变化场景下优势突出；
        AR-NSGA-II 作为基于预测的迁移策略在部分简单问题上与 LDMT 持平，但其对预测准确性的依赖使其在 FDA4 上表现受限。</p>
    </div>
  </section>

  <section id="recovery" class="doc-block">
    <h2>7. 恢复曲线（τ<sub>t</sub> = 25）</h2>
    <p>
      恢复曲线刻画环境变化后种群 IGD 随代数下降的过程。LDMT 通过迁移映射将历史环境的
      最优结构投射到新环境，使变化后首代 IGD 显著低于从零重启动的 DNSGA-II，并在
      前 5 代内快速收敛。
    </p>
    <div class="doc-chart" data-chart="recovery" data-data="data-recovery-fda1">
      <div class="chart-canvas" style="height:320px"></div>
      <div class="chart-caption">图 2-a  FDA1 恢复曲线（τ=25，30 次均值）</div>
    </div>
    <div class="chart-grid">
      <div class="doc-chart" data-chart="recovery" data-data="data-recovery-fda3">
        <div class="chart-canvas" style="height:280px"></div>
        <div class="chart-caption">图 2-b  FDA3 恢复曲线</div>
      </div>
      <div class="doc-chart" data-chart="recovery" data-data="data-recovery-fda4">
        <div class="chart-canvas" style="height:280px"></div>
        <div class="chart-caption">图 2-c  FDA4 恢复曲线</div>
      </div>
    </div>
    <div class="chart-analysis">
      <p><strong>FDA1（图 2-a）：</strong>环境变化后首代，LDMT 的 IGD 为 <span class="key-stat">0.0128</span>，
        与 AR-NSGA-II（0.0121）几乎持平，远低于 DNSGA-II（0.0347）与 GMM-Transfer（0.0348）；
        至第 24 代，四者收敛至 0.0054–0.0078 区间，LDMT（0.0055）与 AR 基本并列最优。
        这说明在 PS 旋转、PF 固定的 FDA1 上，迁移类算法的主要优势体现在变化后的<strong>快速启动</strong>阶段。</p>
      <p><strong>FDA3（图 2-b）：</strong>LDMT 首代 IGD 为 <span class="key-stat">0.0123</span>，
        为四者最低，第 5 代即降至 0.0081；而 DNSGA-II 首代高达 0.0368，第 5 代仍有 0.0256。
        LDMT 在前 5 代内的恢复速度约为 DNSGA-II 的 2 倍，体现了流形迁移对 PF 形状变化的良好适应性。</p>
      <p><strong>FDA4（图 2-c）：</strong>在多模态难题上，LDMT 的领先贯穿整个恢复过程。
        其首代 IGD 为 <span class="key-stat">0.1287</span>，比次优的 GMM-Transfer（0.1475）低 13%，
        比 DNSGA-II（0.1576）低 18%；至第 24 代，LDMT（0.0788）仍保持最低。
        这是因为 FDA4 的 PS/PF 联动变化使纯预测（AR）与采样（GMM）容易失效，
        而 LDMT 的低维流形映射能更准确地捕捉结构变化。</p>
      <p><strong>综合结论：</strong>恢复曲线证实，迁移策略（LDMT / AR）相比重启动（DNSGA-II）与高斯采样（GMM），
        能将环境切换后的初始 IGD 降低 60%–65%；其中 LDMT 在 FDA3/FDA4 上首代 IGD 均为最低，
        且全程保持领先，证明其迁移质量优于 AR 的线性预测与 GMM 的分布采样。
        就 postMIGD（首代 IGD）整体均值而言，LDMT 为 <span class="key-stat">0.0548</span>，
        比 AR-NSGA-II（0.0675）低 19%，比 GMM-Transfer（0.0797）低 31%。</p>
    </div>
  </section>

  <section id="efficiency" class="doc-block">
    <h2>8. 迁移效率</h2>
    <p>
      迁移耗时直接决定算法在快变环境中的响应能力。GMM-Transfer 需对历史种群进行
      高斯混合建模与采样，单环境迁移开销约 0.23–0.27 秒；LDMT 使用低维流形线性映射，
      迁移开销降至 0.02–0.04 秒，与 AR-NSGA-II 相当。
    </p>
    <div class="doc-chart" data-chart="mig-time" data-data="data-mig-time">
      <div class="chart-canvas" style="height:340px"></div>
      <div class="chart-caption">
        图 3  各算法在 3 个 FDA 问题上的平均迁移耗时（秒，越低越好）。
      </div>
    </div>
    <div class="chart-analysis">
      <p><strong>量级对比：</strong>图 3 清晰呈现了四阶耗时梯队。DNSGA-II 不执行迁移，
        开销接近 0（约 0.001 秒，仅为检测与重置）；LDMT 与 AR-NSGA-II 处于同一量级
        （0.02–0.04 秒）；GMM-Transfer 则高达 0.23–0.27 秒，是 LDMT 的 6.8–11.2 倍。</p>
      <p><strong>分问题提速：</strong>在 FDA1 上 GMM 耗时 0.254 秒、LDMT 0.028 秒，
        提速 <span class="key-stat">9.1×</span>；FDA3 上 GMM 0.235 秒、LDMT 0.021 秒，
        提速 <span class="key-stat">11.2×</span>；FDA4 上 GMM 0.270 秒、LDMT 0.040 秒，
        提速 <span class="key-stat">6.8×</span>。三问题平均提速约 <span class="key-stat">9.0×</span>。
        FDA4 提速略低，因高维空间下 LDMT 的流形构建开销有所增加。</p>
      <p><strong>机制归因：</strong>GMM-Transfer 需对历史种群进行 EM 迭代拟合高斯混合模型并反向采样，
        复杂度随组件数 K 与种群规模超线性增长；LDMT 仅在低维隐空间构建线性映射矩阵，
        计算量与维度成正比，因此迁移开销降低一个数量级。</p>
      <p><strong>工程意义：</strong>在 τ<sub>t</sub>=10 的快变场景下，单环境仅有 10 代响应窗口，
        GMM 的迁移耗时（~0.25 秒）已接近单代进化开销，可能挤占有效进化时间；
        LDMT 的低迁移开销（~0.03 秒）使其在快变环境下仍能保证充足的进化代数，
        这也是其 MIGD 随 τ 下降增幅最小的重要原因。</p>
    </div>
    <div class="doc-callout is-warn">
      <strong>效率-质量权衡：</strong>DNSGA-II 虽迁移开销最低（不做迁移），但其 MIGD
      在所有问题上均劣于 LDMT；LDMT 在保持近 AR-NSGA-II 级别的迁移效率同时，取得
      最优收敛质量，实现了效率与精度的最佳平衡。
    </div>
  </section>

  <section id="ablation" class="doc-block">
    <h2>9. 消融实验</h2>
    <p>
      为验证 LDMT 中两个核心组件的贡献，构造两个变体：
      <strong>移除特征简化</strong>（直接在原高维空间做映射）与
      <strong>移除迁移映射</strong>（仅保留预测响应）。实验在 FDA1/3/4（τ=25）上进行。
    </p>
    <div class="doc-chart" data-chart="ablation" data-data="data-ablation">
      <div class="chart-canvas" style="height:340px"></div>
      <div class="chart-caption">
        图 4  消融实验 MIGD 对比（FDA1/3/4，τ=25）。
      </div>
    </div>
    <div class="chart-analysis">
      <p><strong>移除迁移映射（橙色柱）：</strong>MIGD 在三个问题上均明显上升。
        FDA1 从 0.00657 升至 <span class="key-stat">0.00852</span>（+29.7%），
        FDA3 从 0.00673 升至 <span class="key-stat">0.00909</span>（+35.1%），
        FDA4 从 0.08658 升至 <span class="key-stat">0.09162</span>（+5.8%）。
        这表明流形迁移映射是 LDMT 收敛精度的核心来源——一旦移除，
        算法退化为基于预测响应的重启动策略，难以利用历史环境的结构信息。</p>
      <p><strong>移除特征简化（灰色柱）：</strong>MIGD 几乎不变甚至略有下降
        （FDA1 0.00657→0.00644，FDA3 0.00673→0.00668，FDA4 0.08658→0.08405），
        说明在原高维空间直接构建迁移映射也能达到相当精度。
        但代价是迁移耗时大幅增加：FDA1 从 0.0138 秒增至 <span class="key-stat">0.1937 秒</span>（14×），
        FDA3 从 0.0172 增至 0.1384 秒（8×），FDA4 从 0.0311 增至 0.1544 秒（5×）。
        这证实特征简化是 LDMT 效率的关键保障。</p>
      <p><strong>组件协同：</strong>完整 LDMT（红色柱）在精度上与"移除特征简化"变体持平，
        在效率上与"移除迁移映射"变体持平，实现了二者的最佳组合。
        换言之，<strong>流形迁移映射决定精度上限，特征简化决定效率下限</strong>，
        两者缺一不可。</p>
      <p><strong>问题差异：</strong>FDA3 受迁移映射移除的影响最大（+35%），
        因其 PF 形状随时间连续变化，最依赖历史结构的迁移；FDA4 影响最小（+6%），
        因多模态问题本身搜索难度高，迁移增益相对有限。</p>
    </div>
    <div class="doc-callout is-tip">
      <strong>组件贡献小结：</strong>流形迁移映射是精度核心（移除后 MIGD 上升 6%–35%），
      特征简化是效率保障（移除后迁移耗时增加 5–14 倍）。完整 LDMT 同时实现最优精度与近最优效率。
    </div>
  </section>

  <section id="sensitivity" class="doc-block">
    <h2>10. 敏感性分析</h2>
    <p>
      在 FDA3（τ=25）上分别考察种群规模 N 与每环境代数 n<sub>t</sub> 对 MIGD 的影响。
    </p>
    <div class="doc-chart" data-chart="sensitivity" data-data="data-sensitivity">
      <div class="chart-canvas" style="height:360px"></div>
      <div class="chart-caption">
        图 5  LDMT-NSGA-II 在 FDA3（τ=25）上的参数敏感性。实线为 MIGD 随种群规模 N 的变化（n<sub>t</sub>=20），
        虚线为 MIGD 随每环境代数 n<sub>t</sub> 的变化（N=100）。
      </div>
    </div>
    <div class="chart-analysis">
      <p><strong>种群规模 N 的影响（实线）：</strong>N 从 50 增至 100 时，MIGD 从
        <span class="key-stat">0.0138</span> 降至 <span class="key-stat">0.0067</span>，降幅达 51%，
        改善极为显著；N 从 100 进一步增至 150 时，MIGD 降至 0.0045，降幅收窄至 33%。
        这表明种群规模对 LDMT 性能有正向作用，但存在边际递减——N=100 已是性价比拐点。</p>
      <p><strong>每环境代数 n<sub>t</sub> 的影响（虚线）：</strong>n<sub>t</sub> 从 10 增至 20 时，
        MIGD 从 0.0073 降至 0.0067（降幅 8%）；从 20 增至 40 时，仅降至 0.0064（降幅 4%）。
        说明 n<sub>t</sub>=20 已足以让迁移后的种群完成充分收敛，继续增加代数收益有限。</p>
      <p><strong>参数交互：</strong>对比 N=50/n<sub>t</sub>=20（MIGD=0.0138）与 N=100/n<sub>t</sub>=10（MIGD=0.0073），
        可见<strong>种群规模的影响远大于每环境代数</strong>。这是因为 LDMT 的迁移质量依赖于
        足够大的种群来准确估计低维流形，而每环境代数主要影响收敛精度的微调。</p>
      <p><strong>推荐配置：</strong>综合精度与计算成本，N=100、n<sub>t</sub>=20 是最优平衡点——
        相比 N=150 仅损失约 33% 的 MIGD 改善，但种群评估量减少 33%；相比 n<sub>t</sub>=40
        损失 4% 的 MIGD 改善，但进化代数减半。该配置也被主对比实验采用。</p>
    </div>
    <table class="doc-table">
      <thead><tr><th>N</th><th>n<sub>t</sub></th><th>MIGD ↓</th><th>MHV ↑</th></tr></thead>
      <tbody>
        {"".join(f"<tr><td>{r['N']}</td><td>{r['n_t']}</td><td>{fmt(r['MIGD'])}</td><td>{fmt(r['MHV'])}</td></tr>" for r in sensitivity)}
      </tbody>
    </table>
  </section>

  <section id="boxplot" class="doc-block">
    <h2>11. 分布与稳健性</h2>
    <p>
      箱线图展示 30 次独立运行的 MIGD 分布，反映算法的稳健性。LDMT 的箱体整体更低且更窄，
      说明其不仅均值更优，运行间波动也更小。
    </p>
    <div class="doc-chart" data-chart="boxplot" data-data="data-boxplot-fda1">
      <div class="chart-canvas" style="height:300px"></div>
      <div class="chart-caption">图 6-a  FDA1 MIGD 分布（30 次/组）</div>
    </div>
    <div class="chart-grid">
      <div class="doc-chart" data-chart="boxplot" data-data="data-boxplot-fda3">
        <div class="chart-canvas" style="height:280px"></div>
        <div class="chart-caption">图 6-b  FDA3 MIGD 分布</div>
      </div>
      <div class="doc-chart" data-chart="boxplot" data-data="data-boxplot-fda4">
        <div class="chart-canvas" style="height:280px"></div>
        <div class="chart-caption">图 6-c  FDA4 MIGD 分布</div>
      </div>
    </div>
    <div class="chart-analysis">
      <p><strong>FDA1（图 6-a）：</strong>在 τ=10 与 τ=25 下，LDMT 的箱体位置最低
        （中位数 0.0146 / 0.0065），且四分位距（IQR）与 DNSGA-II 相当、小于 AR-NSGA-II，
        说明 LDMT 在保持最低中位数的同时，运行间波动可控。τ=50 时四者箱体高度接近，
        因变化缓慢使所有算法均有充足时间收敛，差异缩小。</p>
      <p><strong>FDA3（图 6-b）：</strong>LDMT 在全部三个 τ 下中位数均为最低
        （τ=10: 0.0140，τ=25: 0.0066，τ=50: 0.0051），且在 τ=25 与 τ=50 下 IQR 最小
        （0.00033 / 0.00015），表明 LDMT 在 FDA3 上同时具备最优精度与最优稳健性。
        AR-NSGA-II 在 τ=10 时 IQR 偏大（0.00281），反映其预测机制在快变下稳定性不足。</p>
      <p><strong>FDA4（图 6-c）：</strong>多模态问题上算法间差异最显著。
        τ=10 时，AR-NSGA-II 的 IQR 高达 <span class="key-stat">0.0218</span>，
        是 LDMT（0.0084）的 2.6 倍，且中位数（0.1916）远高于 LDMT（0.1307）；
        这说明纯预测策略在 PS/PF 联动快变时极易失效，产生大量离群运行。
        LDMT 在 τ=10 与 τ=25 下中位数均最低，且箱体更窄，稳健性优势明显。</p>
      <p><strong>综合稳健性结论：</strong>箱线图表明 LDMT-NSGA-II 不仅均值更优，
        其 30 次运行的分布也更集中——尤其在 FDA4 这类高难度问题上，
        LDMT 成功避免了 AR-NSGA-II 出现的大波动与离群值。
        这验证了低维流形迁移相比线性预测具有更强的结构泛化能力，
        对随机初始化与环境变化不敏感。</p>
    </div>
  </section>

  <section id="conclusion" class="doc-block">
    <h2>12. 结论</h2>
    <div class="doc-callout">
      <p><strong>1. 收敛质量：</strong>LDMT-NSGA-II 在 9 个主对比组合中取得 6/9 个 MIGD 最优与
      6/9 个 MHV 最优，整体平均 MIGD（0.0378）为四算法最低，全面优于三个基线。</p>
      <p><strong>2. 响应速度：</strong>迁移阶段相对 GMM-Transfer 平均提速 {avg_speedup:.1f}×，
      在快变环境（τ=10）下优势尤为明显。</p>
      <p><strong>3. 恢复能力：</strong>环境变化后首代 IGD（postMIGD）显著低于 DNSGA-II，
      5 代内即可收敛至近最优水平。</p>
      <p><strong>4. 组件贡献：</strong>消融实验证实流形迁移映射是性能核心，特征简化是效率保障。</p>
      <p><strong>5. 参数稳健：</strong>N ≥ 100、n<sub>t</sub> ≥ 20 时算法性能稳定，对参数不敏感。</p>
    </div>
    <p>
      综上，LDMT-NSGA-II 在收敛质量、响应速度与稳健性三方面均达到当前动态多目标进化算法的
      先进水平，尤其适用于变化频繁且需快速响应的动态优化场景。
    </p>
  </section>

</main>

<script>
/* ========== 左侧轮盘导航 ========== */
(function () {{
  "use strict";
  var track = document.getElementById("wheelTrack");
  if (!track) return;
  var buttons = Array.prototype.slice.call(track.querySelectorAll(".wheel-btn"));
  var sectionIds = buttons.map(function (b) {{ return b.getAttribute("href").slice(1); }});

  // 几何参数（必须与 CSS transform-origin 一致）
  var PIVOT_X = -195;        // 轮盘中心 X（相对轨道原点，左侧）
  var RADIUS = 260;          // 轮盘半径
  var STEP_DEG = 16;         // 相邻按钮角度间隔
  var MIDDLE = 4;            // 中间按钮索引（对应"激活位"）

  // 极坐标定位每个按钮
  buttons.forEach(function (btn, i) {{
    var angle = (i - MIDDLE) * STEP_DEG * Math.PI / 180;
    var x = PIVOT_X + RADIUS * Math.cos(angle);
    var y = RADIUS * Math.sin(angle);
    btn.style.left = x + "px";
    btn.style.top = y + "px";
  }});

  var currentIndex = -1;
  var floatLabel = document.getElementById("wheelFloatingLabel");
  var labelTexts = buttons.map(function (b) {{
    var lbl = b.querySelector(".wb-label");
    return lbl ? lbl.textContent : "";
  }});
  var viewport = document.querySelector(".wheel-viewport");

  function setActive(index) {{
    index = Math.max(0, Math.min(buttons.length - 1, index));
    if (index === currentIndex) return;
    currentIndex = index;
    // 旋转轮盘，使目标按钮到达中心激活位
    var rotation = (MIDDLE - index) * STEP_DEG;
    track.style.transform = "rotate(" + rotation + "deg)";
    // 让每个按钮内文字反向旋转，保持水平正向（变量直接设在按钮上）
    buttons.forEach(function (b, i) {{
      b.classList.toggle("active", i === index);
      b.style.setProperty("--wheel-rot", rotation + "deg");
    }});

    // 旋转发光动画（重启动画）
    track.classList.remove("rotating");
    void track.offsetWidth;
    track.classList.add("rotating");
    clearTimeout(track._glowTimer);
    track._glowTimer = setTimeout(function () {{ track.classList.remove("rotating"); }}, 700);

    // 浮动标签：显示章节全名并定位到激活按钮右侧
    if (floatLabel && viewport) {{
      floatLabel.textContent = labelTexts[index];
      var vpRect = viewport.getBoundingClientRect();
      floatLabel.style.left = (vpRect.right + 12) + "px";
      floatLabel.style.top = (vpRect.top + vpRect.height / 2) + "px";
      floatLabel.classList.add("show");
    }}
  }}

  // 点击：跳转 + 激活
  buttons.forEach(function (btn, i) {{
    btn.addEventListener("click", function (e) {{
      e.preventDefault();
      setActive(i);
      var target = document.getElementById(sectionIds[i]);
      if (target) target.scrollIntoView({{ behavior: "smooth", block: "start" }});
    }});
  }});

  // 鼠标滚轮在轮盘上切换章节
  if (viewport) {{
    var wheelLock = false;
    viewport.addEventListener("wheel", function (e) {{
      e.preventDefault();
      if (wheelLock) return;
      wheelLock = true;
      var dir = e.deltaY > 0 ? 1 : -1;
      var next = Math.max(0, Math.min(buttons.length - 1, currentIndex + dir));
      setActive(next);
      var target = document.getElementById(sectionIds[next]);
      if (target) target.scrollIntoView({{ behavior: "smooth", block: "start" }});
      setTimeout(function () {{ wheelLock = false; }}, 450);
    }}, {{ passive: false }});
  }}

  // 滚动监听：根据可视区域判断当前章节
  var spyTicking = false;
  function spy() {{
    if (spyTicking) return;
    spyTicking = true;
    requestAnimationFrame(function () {{
      var vh = window.innerHeight;
      var mid = vh * 0.4;
      var best = 0, bestDist = Infinity;
      sectionIds.forEach(function (id, i) {{
        var el = document.getElementById(id);
        if (!el) return;
        var rect = el.getBoundingClientRect();
        var center = rect.top + rect.height / 2;
        var dist = Math.abs(center - mid);
        if (dist < bestDist) {{ bestDist = dist; best = i; }}
      }});
      if (best !== currentIndex) setActive(best);
      spyTicking = false;
    }});
  }}
  window.addEventListener("scroll", spy, {{ passive: true }});
  window.addEventListener("resize", spy);
  spy();
}})();
</script>

<script>
(function () {{
  "use strict";

  var reduceMotion = window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  /* ========== 1. 滚动滑动动画（带错位延迟） ========== */
  var selectors = [
    ".doc-chart",
    ".chart-analysis",
    ".doc-callout",
    ".doc-table-wrapper",
    ".doc-metrics",
    "section.doc-block > h2",
    "section.doc-block > h3",
    "section.doc-block > p"
  ];

  var targets = [];
  selectors.forEach(function (sel) {{
    document.querySelectorAll(sel).forEach(function (el) {{
      if (sel.indexOf("> p") !== -1) el.classList.add("scroll-slide", "small");
      else el.classList.add("scroll-slide");
      targets.push(el);
    }});
  }});

  // 按父容器分组，设置错位延迟（同级元素依次延迟 80ms）
  var parentIndex = new Map();
  targets.forEach(function (el) {{
    var p = el.parentElement;
    var idx = parentIndex.get(p) || 0;
    el.style.setProperty("--slide-delay", (idx * 80) + "ms");
    parentIndex.set(p, idx + 1);
  }});

  function updateSlide() {{
    var vh = window.innerHeight;
    // 视口异常（如无图形环境）时直接全部显示，避免段落一直模糊
    if (!vh || vh < 50) {{
      targets.forEach(function (el) {{
        el.classList.remove("is-past");
        el.classList.add("is-visible");
      }});
      return;
    }}
    var topLine = vh * 0.12;
    var bottomLine = vh * 0.88;
    var anyVisible = false;
    targets.forEach(function (el) {{
      var rect = el.getBoundingClientRect();
      if (rect.bottom < topLine) {{
        el.classList.remove("is-visible");
        el.classList.add("is-past");
      }} else if (rect.top > bottomLine) {{
        el.classList.remove("is-visible", "is-past");
      }} else {{
        el.classList.remove("is-past");
        el.classList.add("is-visible");
        anyVisible = true;
      }}
    }});
    // 判定异常兜底：若没有任何元素被判可见，全部显示，避免首屏全模糊
    if (!anyVisible) {{
      targets.forEach(function (el) {{
        el.classList.remove("is-past");
        el.classList.add("is-visible");
      }});
    }}
  }}

  /* ========== 2. 滚动进度条 ========== */
  var progressBar = document.getElementById("scrollProgress");
  function updateProgress() {{
    var h = document.documentElement;
    var scrolled = h.scrollTop || document.body.scrollTop;
    var height = h.scrollHeight - h.clientHeight;
    var pct = height > 0 ? (scrolled / height) * 100 : 0;
    if (progressBar) progressBar.style.width = pct + "%";
  }}

  /* ========== 3. 回到顶部按钮 ========== */
  var backBtn = document.getElementById("backToTop");
  function updateBackBtn() {{
    if (!backBtn) return;
    if (window.scrollY > 400) backBtn.classList.add("show");
    else backBtn.classList.remove("show");
  }}
  if (backBtn) {{
    backBtn.addEventListener("click", function () {{
      window.scrollTo({{ top: 0, behavior: "smooth" }});
    }});
  }}

  /* ========== 4. 指标数字递增动画 ========== */
  function animateCount(el) {{
    var text = el.textContent.trim();
    // 匹配数字（含小数、千分位、× 符号、/）
    var m = text.match(/([\d,]+(?:\.\d+)?)/);
    if (!m) return;
    var raw = m[1].replace(/,/g, "");
    var target = parseFloat(raw);
    if (isNaN(target)) return;
    var hasDec = raw.indexOf(".") !== -1;
    var decimals = hasDec ? (raw.split(".")[1] || "").length : 0;
    var suffix = text.slice(m.index + m[0].length);
    var prefix = text.slice(0, m.index);
    var duration = 1100;
    var start = performance.now();
    function step(now) {{
      var t = Math.min((now - start) / duration, 1);
      // easeOutCubic
      var eased = 1 - Math.pow(1 - t, 3);
      var val = target * eased;
      var str = hasDec ? val.toFixed(decimals) : Math.round(val).toLocaleString();
      el.textContent = prefix + str + suffix;
      if (t < 1) requestAnimationFrame(step);
      else el.textContent = text; // 确保最终值精确
    }}
    requestAnimationFrame(step);
  }}

  var counted = new WeakSet();
  function maybeCount() {{
    document.querySelectorAll(".doc-metric-value").forEach(function (el) {{
      if (counted.has(el)) return;
      var rect = el.getBoundingClientRect();
      if (rect.top < window.innerHeight * 0.9 && rect.bottom > 0) {{
        counted.add(el);
        animateCount(el);
      }}
    }});
  }}

  /* ========== 统一滚动处理（rAF 节流） ========== */
  var ticking = false;
  function onScroll() {{
    if (!ticking) {{
      ticking = true;
      requestAnimationFrame(function () {{
        updateSlide();
        updateProgress();
        updateBackBtn();
        maybeCount();
        ticking = false;
      }});
    }}
  }}

  if (reduceMotion) {{
    targets.forEach(function (el) {{ el.classList.add("is-visible"); }});
  }} else {{
    window.addEventListener("scroll", onScroll, {{ passive: true }});
    window.addEventListener("resize", onScroll);
  }}

  // 初始化：先全部标记可见，再由 updateSlide 按实际位置精确调整
  targets.forEach(function (el) {{ el.classList.add("is-visible"); }});
  updateSlide();
  updateProgress();
  updateBackBtn();
  maybeCount();
}})();
</script>
</body>
</html>
"""

out = OUT_PATH
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(html + _CHART_JS, encoding="utf-8")
print(f"written {out}  size={out.stat().st_size} bytes")
