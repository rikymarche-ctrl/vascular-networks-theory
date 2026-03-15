"""
run_all.py -- Master script for the branching-networks monorepo.

Regenerates all dynamic variables and figures for all papers by running
each compute script in sequence. Called by build.ps1 -Target compute (or all).

Usage:
    python run_all.py
"""

import subprocess
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def run(script_name):
    """Run a script in SCRIPT_DIR, inherit stdout/stderr, raise on failure."""
    script_path = os.path.join(SCRIPT_DIR, script_name)
    print(f"\n{'=' * 60}")
    print(f"Running: {script_name}")
    print('=' * 60)
    result = subprocess.run([sys.executable, script_path], cwd=SCRIPT_DIR)
    if result.returncode != 0:
        print(f"\n[ERROR] {script_name} exited with code {result.returncode}", file=sys.stderr)
        sys.exit(result.returncode)


if __name__ == '__main__':
    run('compute_paper_murray.py')
    run('compute_paper_variational.py')
    run('compute_supplemental.py')
    print(f"\n{'=' * 60}")
    print("[OK] All computations complete.")
    print('=' * 60)
