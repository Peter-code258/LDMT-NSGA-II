"""LDT-NSGA-II 一键流水线（v2.4）。

按项目标准顺序自动执行全部流程：

    1. pytest -q            63 个单元测试（质量门，失败即中止）
    2. analyze.py           统计汇总 + 统计检验 + 12 张 PNG 图表
    3. build_report.py      重新生成单文件自包含 HTML 实验报告
    4. 自动用默认浏览器打开  ldmt-experiment-report/index.html

用法（在项目根目录）：

    python main.py                      # 测试 → 分析 → 报告 → 打开网页
    python main.py --skip-tests         # 跳过 pytest（仅分析+报告）
    python main.py --with-experiments   # 先重跑四组实验（约 3780 次运行，
                                        #   断点续跑，耗时以小时计）
    python main.py --no-open            # 不自动打开浏览器

说明：v2.4 的实验数据（v2.3 K1 修复后重跑版）已完整落盘于 results/，
日常更新文档/报告无需重跑实验，故实验组默认不执行。
"""
import argparse
import subprocess
import sys
import time
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPORT = ROOT / "ldmt-experiment-report" / "index.html"

EXPERIMENT_GROUPS = [
    ("主比较 comparison", [sys.executable, "-m", "experiments.run_comparison"]),
    ("DMOP 补充 supplementary", [sys.executable, "-m", "experiments.run_supplementary"]),
    ("消融 ablation", [sys.executable, "-m", "experiments.run_ablation"]),
    ("敏感性 sensitivity", [sys.executable, "-m", "experiments.run_sensitivity"]),
]


def run_step(title, cmd):
    """Run one pipeline step; abort the whole pipeline on failure."""
    print(f"\n{'=' * 64}\n[main] {title}\n  $ {' '.join(str(c) for c in cmd)}\n{'-' * 64}")
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT)
    dt = time.perf_counter() - t0
    if proc.returncode != 0:
        print(f"\n[main] ✗ 步骤失败（exit {proc.returncode}）：{title}")
        print("[main] 流水线中止；后续步骤未执行。请先修复上面的错误。")
        sys.exit(proc.returncode)
    print(f"[main] ✓ {title} 完成（{dt:.1f} s）")


def main():
    ap = argparse.ArgumentParser(
        description="LDT-NSGA-II 一键流水线：测试 → analyze → build_report → 打开报告")
    ap.add_argument("--skip-tests", action="store_true", help="跳过 pytest 质量门")
    ap.add_argument("--with-experiments", action="store_true",
                    help="先重跑四组实验（约 3780 次运行，断点续跑）")
    ap.add_argument("--no-open", action="store_true", help="生成后不打开浏览器")
    args = ap.parse_args()

    print("[main] LDT-NSGA-II v2.4 一键流水线")
    print(f"[main] 项目根目录：{ROOT}")

    if args.with_experiments:
        for title, cmd in EXPERIMENT_GROUPS:
            run_step(f"实验组：{title}", cmd)
    else:
        print("[main] 跳过四组实验（--with-experiments 可开启；v2.4 数据已完整落盘）")

    if not args.skip_tests:
        run_step("pytest 单元测试（63 个）", [sys.executable, "-m", "pytest", "tests", "-q"])

    run_step("analyze.py（统计汇总 + 图表）", [sys.executable, "analyze.py"])
    run_step("build_report.py（HTML 实验报告）", [sys.executable, "build_report.py"])

    if not REPORT.exists():
        print(f"\n[main] ✗ 报告文件不存在：{REPORT}")
        sys.exit(1)
    size_mb = REPORT.stat().st_size / 1024 / 1024
    print(f"\n[main] 报告已生成：{REPORT}（{size_mb:.2f} MB，单文件自包含）")

    if not args.no_open:
        print("[main] 正在用默认浏览器打开报告 …")
        webbrowser.open(REPORT.as_uri())

    print("\n[main] 全部流程完成 ✓")


if __name__ == "__main__":
    main()
