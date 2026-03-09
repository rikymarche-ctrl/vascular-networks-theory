#!/usr/bin/env python3
"""
reproduce.py — Cross-platform reproduction script (alternative to `make`).

Usage:
    python reproduce.py           # both papers
    python reproduce.py paper1    # Paper 1 only
    python reproduce.py paper2    # Paper 2 only

Requirements: Python 3.8+, numpy, scipy, matplotlib, pdflatex (optional).
Install Python deps: pip install -r requirements.txt
"""

import os
import sys
import subprocess

ROOT = os.path.dirname(os.path.abspath(__file__))

PAPER1_TITLE = "Beyond Murray's Law"
PAPER2_TITLE = "A Unified Variational Principle for Branching Transport Networks"


def run(cmd, cwd=None):
    print(f"  > {cmd}")
    result = subprocess.run(cmd, cwd=cwd, shell=True)
    if result.returncode != 0:
        print(f"\nError: command failed with exit code {result.returncode}")
        sys.exit(result.returncode)


def has_pdflatex():
    return subprocess.run(
        "pdflatex --version", shell=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    ).returncode == 0


def compute(script_rel):
    run(f'python "{os.path.join(ROOT, script_rel)}"')


def compile_latex(paper_dir, jobname, bib=False):
    mdir = os.path.join(ROOT, paper_dir, "manuscript")
    if not os.path.isdir(mdir):
        print(f"  Skipping LaTeX: {mdir} not found")
        return
    flags = "-interaction=nonstopmode"
    jflag = f'-jobname "{jobname}"'
    run(f'pdflatex {flags} {jflag} main.tex', cwd=mdir)
    if bib:
        run(f'bibtex "{jobname}"', cwd=mdir)
        run(f'pdflatex {flags} {jflag} main.tex', cwd=mdir)
    run(f'pdflatex {flags} {jflag} main.tex', cwd=mdir)


def paper1():
    print("\n=== Paper 1: Beyond Murray's Law ===")
    compute("paper1-murray/compute_paper_murray.py")
    if has_pdflatex():
        compile_latex("paper1-murray", jobname=PAPER1_TITLE, bib=True)
    else:
        print("  pdflatex not found — skipping LaTeX compilation")


def paper2():
    print("\n=== Paper 2: Unified Variational Principle ===")
    compute("paper2-variational/compute_paper_variational.py")
    if has_pdflatex():
        compile_latex("paper2-variational", jobname=PAPER2_TITLE, bib=True)
    else:
        print("  pdflatex not found — skipping LaTeX compilation")


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "all"
    if target == "paper1":
        paper1()
    elif target == "paper2":
        paper2()
    elif target == "all":
        paper1()
        paper2()
    else:
        print(f"Unknown target '{target}'. Use: paper1, paper2, or omit for all.")
        sys.exit(1)
    print("\nDone.")
