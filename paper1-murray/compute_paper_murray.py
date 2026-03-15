"""
compute_paper_murray.py — Entry point for Paper 1.

Delegates to shared/scripts/compute_paper_murray.py, which contains
the full generalized cost function model. This wrapper exists so that
`make paper1` and `python reproduce.py paper1` work from the repository
root without requiring manual PYTHONPATH configuration.
"""

import subprocess
import sys
import os

if __name__ == "__main__":
    script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "shared", "scripts", "compute_paper_murray.py"
    )
    result = subprocess.run([sys.executable, os.path.normpath(script)])
    sys.exit(result.returncode)
