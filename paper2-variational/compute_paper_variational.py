"""
compute_paper_variational.py — Entry point for Paper 2.

Delegates to shared/scripts/compute_paper_variational.py, which contains
the full two-level minimax model. This wrapper exists so that `make paper2`
and `python reproduce.py paper2` work from the repository root without
requiring manual PYTHONPATH configuration.
"""

import subprocess
import sys
import os

if __name__ == "__main__":
    script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "shared", "scripts", "compute_paper_variational.py"
    )
    result = subprocess.run([sys.executable, os.path.normpath(script)])
    sys.exit(result.returncode)
