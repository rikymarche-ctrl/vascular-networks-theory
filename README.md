# Vascular Networks Theory

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![arXiv Paper 1](https://img.shields.io/badge/arXiv-2603.13687-b31b1b.svg)](https://arxiv.org/abs/2603.13687)
[![arXiv Paper 2](https://img.shields.io/badge/arXiv-2603.14691-b31b1b.svg)](https://arxiv.org/abs/2603.14691)
[![arXiv Paper 3](https://img.shields.io/badge/arXiv-Submitted-b31b1b.svg)]()
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

A unified theoretical framework for biological branching transport networks.
Starting from measured tissue properties alone — with **zero fitted parameters** — the framework derives the empirically observed branching exponent α ≈ 2.7 from first principles, resolving a long-standing discrepancy with Murray's classical cubic law.

---

## Papers

| # | Title | Links | Key result |
|---|---|---|---|
| **1** | [*Beyond Murray's Law: Non-Universal Branching Exponents from Vessel-Wall Metabolic Costs*](paper1-murray/output/) | [arXiv:2603.13687](https://arxiv.org/abs/2603.13687) | Wall sub-linearity (*p* < 1) rigorously breaks Murray's universality; α\* ∈ [2.90, 2.94] |
| **2** | [*A Unified Variational Principle for Branching Transport Networks*](paper2-variational/output/) | [arXiv:2603.14691](https://arxiv.org/abs/2603.14691) | Minimax game between wave and transport costs yields **α\* = 2.72** with zero free parameters |
| **3** | [*The Dynamic Origin of Kleiber's Law*](paper3-kleiber/output/) | [arXiv:Submitted]() | Parameter-free exact prediction of macroscopic metabolic scaling laws across five phyla |

---

## Quickstart

```powershell
# 1. Clone and install dependencies
git clone https://github.com/rikymarche-ctrl/vascular-networks-theory.git
cd vascular-networks-theory
pip install -r requirements.txt

# 2. Reproduce everything (Python + LaTeX)
.\build.ps1

# Or one paper at a time:
.\build.ps1 -Target paper1
.\build.ps1 -Target paper2
.\build.ps1 -Target paper3
```

> **Note:** The `build.ps1` script acts as a master dispatcher, delegating to the modular PowerShell scripts inside each paper's directory. LaTeX compilation requires a standard TeX distribution (TeX Live, MacTeX, or MiKTeX) with `pdflatex` and `bibtex`.

---

## Repository Structure

```
vascular-networks-theory/
├── build.ps1                           # Master build pipeline
├── requirements.txt
├── CITATION.cff
├── LICENSE
├── shared/
│   └── scripts/
│       └── params.py                   # Single source of truth for all physical constants
├── paper1-murray/
│   ├── build.ps1                       # Local build script
│   ├── scripts/
│   │   └── compute.py                  # Generates dynamic variables
│   ├── manuscript/
│   │   ├── main.tex                    # LaTeX source
│   │   ├── references.bib
│   │   ├── dynamic_variables.tex       # Auto-generated — do not edit by hand
│   │   └── figures/
│   └── output/                         # Compiled PDFs
└── paper2-variational/
    ├── build.ps1                       # Local build script
    ├── scripts/
    │   └── compute.py                  # Generates dynamic variables
    ├── manuscript/
    ├── supplements/
    └── output/                         # Compiled PDFs

└── paper3-kleiber/
    ├── build.ps1                       # Local build script
    ├── scripts/
    │   └── compute.py                  # Generates dynamic variables
    ├── manuscript/
    └── output/                         # Compiled PDFs
```

---

## Reproducibility Pipeline

Every number that appears in the manuscripts is generated automatically:

1. **`shared/scripts/params.py`** defines all physical constants, each traced to its published source.
2. **`compute_paper_*.py`** runs deterministic numerical computations and writes `dynamic_variables.tex`.
3. **LaTeX** reads every quantity via `\input{dynamic_variables.tex}` — no value is typed by hand.

Modifying a single parameter in `params.py` and re-running the pipeline updates all tables, figures, and inline values atomically.

---

## Citation

If this work is useful to you, please cite the relevant paper:

```bibtex
@article{marchesi2026murray,
  author  = {Marchesi, Riccardo},
  title   = {Beyond {Murray's} Law: Non-Universal Branching Exponents
             from Vessel-Wall Metabolic Costs},
  journal = {arXiv preprint arXiv:2603.13687},
  year    = {2026},
  url     = {https://arxiv.org/abs/2603.13687}
}

@article{marchesi2026variational,
  author  = {Marchesi, Riccardo},
  title   = {A Unified Variational Principle for Branching Transport Networks},
  journal = {arXiv preprint arXiv:2603.14691},
  year    = {2026},
  url     = {https://arxiv.org/abs/2603.14691}
}

@article{marchesi2026kleiber,
  author  = {Marchesi, Riccardo},
  title   = {The Dynamic Origin of Kleiber's Law},
  journal = {arXiv preprint arXiv:Submitted},
  year    = {2026},
  url     = {}
}
```

---

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
You are free to share and adapt it for any purpose, provided appropriate credit is given.
