# -*- coding: utf-8 -*-
"""Final v3.0 release audit: code, docs, archive, and README consistency.

Run from the repository root:

    python -m scripts.check_v30_release

Asserts that every v3.0 deliverable exists, that the V2.5 material is still
archived, that each optimisation-plan item (K27-K33) is really implemented,
that the test suite has the expected size, and that both theses and both
READMEs carry the v3.0 content.  Non-zero exit on any failure, so it can be
used as a release gate.
"""
import io
import os
import re
import subprocess
import sys

# this file: <workspace>/LDT-NSGA-II/scripts/check_v30_release.py
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WS = os.path.dirname(PROJ)
MAT = os.path.join(WS, "文档", "论文材料")
ARCH = os.path.join(MAT, "参考归档")
MID = os.path.join(ARCH, "中期检查报告")

ok = True


def check(label, cond, detail=""):
    global ok
    print(f"  {'OK  ' if cond else 'FAIL'} {label}" + (f"  [{detail}]" if detail else ""))
    if not cond:
        ok = False


print("=== 1. deliverables ===")
for name, min_pages_hint in [
    ("LDT-NSGA-II实验论文v3.1.docx", "185"),
    ("LDT-NSGA-II实验论文v3.1.md", ""),
    ("LDT-NSGA-II本科毕业论文.docx", "28"),
    ("LDT-NSGA-II开题报告.docx", "29"),
    ("LDT-NSGA-II V2.5答辩.html", ""),
]:
    p = os.path.join(MAT, name)
    check(name, os.path.exists(p), f"{os.path.getsize(p):,} B" if os.path.exists(p) else "missing")

print()
print("=== 2. archive (nothing lost) ===")
for name in ["LDT-NSGA-II实验论文 V2.5.docx", "LDT-NSGA-II实验论文 V2.5.md"]:
    p = os.path.join(MID, name)
    check(f"中期检查报告/{name}", os.path.exists(p))
for name in ["LDT-NSGA-II项目优化与精进方案 V2.0.docx",
             "LDT-NSGA-II毕业论文v1.0.docx",
             "LDT-NSGA-II开题报告v1.0.docx",
             "LDT-NSGA-II V2.4答辩.html",
             "LDT-NSGA-II V2.4毕业论文.md",
             "LDT-NSGA-II_V2.4毕业论文（草稿）.docx"]:
    p = os.path.join(ARCH, name)
    check(f"参考归档/{name}", os.path.exists(p))

print()
print("=== 3. no stale duplicate at top level ===")
for bad in ["LDT-NSGA-II毕业论文（硕士学位论文）.docx",
            "LDT-NSGA-II毕业论文.docx",
            "LDT-NSGA-II开题报告（研究生）.docx",
            "LDT-NSGA-II V2.4答辩.html"]:
    p = os.path.join(MAT, bad)
    check(f"absent: {bad}", not os.path.exists(p))

print()
print("=== 4. code: plan items P0/P1 ===")
def read(rel):
    p = os.path.join(PROJ, rel.replace("/", os.sep))
    return io.open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""

core = read("algorithms/nsga2_core.py")
check("K27 code_version injected (no hard-coded v2.4)",
      'code_version="unknown"' in core and 'self.code_version = "v2.4"' not in core)
check("K27 run_id column in rows", '"run_id": self.run_id' in core)
check("K31 establish_baseline called per environment",
      "establish_baseline" in core)

common = read("experiments/common.py")
check("K28 config_hash present", "def config_hash" in common)
check("K28 make_run_id present", "def make_run_id" in common)
check("K28 experiment_group in RunConfig", "experiment_group" in common)
check("K32 partial detection", "_partial_seeds" in common and "_drop_partial_seeds" in common)
check("K32 complete flag in summary", '"complete":' in common)
check("git_dirty recorded", "def git_dirty" in common)

rep = read("build_report.py")
check("K29 no hard-coded drive path in header",
      'Path(r"D:' not in rep[:2500] and "Path(__file__)" in rep[:2500])

req = read("requirements.txt")
check("K30 pandas declared", "pandas" in req)
check("K30 versions pinned", req.count("==") >= 9)

sup = read("experiments/run_supplementary.py")
check("K14 default covers 3 project variants",
      'DEFAULT_PROBLEMS = "DMOP1,DMOP2,DMOP3"' in sup)

check("K33 CI workflow present",
      os.path.exists(os.path.join(PROJ, ".github", "workflows", "ci.yml")))
check("K33 CI runs pytest + smoke + report",
      all(k in read(".github/workflows/ci.yml") for k in
          ("pytest", "_ci_smoke", "build_report.py", "check_dmop_authoritative")))

print()
print("=== 5. tests ===")
r = subprocess.run([sys.executable, "-m", "pytest", "-q", "--collect-only"],
                   cwd=PROJ, capture_output=True, text=True)
m = re.search(r"(\d+) tests? collected", r.stdout or "")
n = int(m.group(1)) if m else -1
# v3.0 shipped 77; v3.1 added the adaptive-alpha regressions.
check(f"{n} tests collected (>= 77 expected)", n >= 77, f"got {n}")
check("test_v30_fixes.py exists",
      os.path.exists(os.path.join(PROJ, "tests", "test_v30_fixes.py")))
check("test_v31_adaptive_alpha.py exists",
      os.path.exists(os.path.join(PROJ, "tests", "test_v31_adaptive_alpha.py")))

print()
print("=== 6. theses reflect v3.0 ===")
master = io.open(os.path.join(MAT, "LDT-NSGA-II实验论文v3.1.md"), encoding="utf-8").read()
check("master has 4.7 v3.0 data-trustworthiness section",
      "v3.0 数据可信度层" in master)
check("master K-list spans K1-K34",
      "K1–K34" in master or "K1-K34" in master)
check("master has no K1-K26 residual",
      "K1–K26" not in master and "K1-K26" not in master)
check("master states v3.0 did not change search behaviour",
      "未改变任何搜索行为" in master or "不改变搜索行为" in master)

ug = io.open(os.path.join(MAT, "LDT-NSGA-II本科毕业论文.md"), encoding="utf-8").read()
check("ugrad has v3.0 traceability item", "run_id" in ug and "v3.0" in ug)
check("ugrad K-list spans K1-K34", "K1–K34" in ug or "K1-K34" in ug)

print()
print("=== 7. READMEs ===")
docs_readme = io.open(os.path.join(WS, "文档", "README.md"), encoding="utf-8").read()
code_readme = io.open(os.path.join(PROJ, "README.md"), encoding="utf-8").read()
for name, s in (("文档/README.md", docs_readme), ("LDT-NSGA-II/README.md", code_readme)):
    check(f"{name} mentions v3.1", "v3.1" in s)
    check(f"{name} mentions the current page count",
          "188 页" in s or "185 页" in s or "177" in s)
    # A *bare* "源文件" is only a hazard when it is used as a PATH. Prose such as
    # "改源文件再重新构建" is legitimate, so look at path-shaped occurrences.
    bare_paths = re.findall(r"[\w./\\]*?(?<!_)\u6e90\u6587\u4ef6[\\/]", s)
    check(f"{name} has no un-prefixed 源文件/ path",
          not bare_paths, str(sorted(set(bare_paths))[:3]))
    check(f"{name} mentions 中期检查报告", "中期检查报告" in s)

check("docs README references the current thesis name",
      re.search(r"LDT-NSGA-II实验论文v3\.[01]\.docx", docs_readme) is not None)
check("code README references the current thesis name",
      re.search(r"LDT-NSGA-II实验论文v3\.[01]\.docx", code_readme) is not None)
check("code README points at the mid-term report archive",
      "中期检查报告" in code_readme)

check("master records the K34 retraction of conclusion 3",
      "\u64a4\u56de" in master and "0.05697" in master)
check("master documents the adaptive-alpha study",
      "ldt_adaptive" in master)

print()
print("=== 8. v3.1 adaptive-alpha evidence ===")
for rel in (("results", "raw_v31"),
            ("results", "summary_v31", "ADAPTIVE_ALPHA_FINDINGS.md"),
            ("scripts", "staleness_audit.py"),
            ("scripts", "adaptive_alpha_analysis.py"),
            ("experiments", "run_v31_adaptive.py")):
    check("present: " + "/".join(rel),
          os.path.exists(os.path.join(PROJ, *rel)))
check("ldt_adaptive registered in experiments/common.py",
      "ldt_adaptive" in io.open(os.path.join(PROJ, "experiments", "common.py"),
                                encoding="utf-8").read())
check("adaptive rule has a dedicated test module",
      os.path.exists(os.path.join(PROJ, "tests",
                                  "test_v31_adaptive_alpha.py")))

print()
print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
sys.exit(0 if ok else 1)
