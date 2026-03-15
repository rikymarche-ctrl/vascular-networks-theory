# Vascular Networks Theory

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI Paper 1](https://zenodo.org/badge/DOI/10.5281/zenodo.18975093.svg)](https://zenodo.org/records/18975093)
[![DOI Paper 2](https://zenodo.org/badge/DOI/10.5281/zenodo.19027393.svg)](https://zenodo.org/records/19027393)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

A unified theoretical framework for biological branching transport networks.
Starting from measured tissue properties alone — with **zero fitted parameters** — the framework derives the empirically observed branching exponent α ≈ 2.7 from first principles, resolving a long-standing discrepancy with Murray's classical cubic law.

---

## Papers

| # | Title | Key result |
|---|---|---|
| **1** | [*Beyond Murray's Law: Non-Universal Branching Exponents from Vessel-Wall Metabolic Costs*](paper1-murray/manuscript/) | Wall sub-linearity (*p* < 1) rigorously breaks Murray's universality; α\* ∈ [2.90, 2.94] |
| **2** | [*A Unified Variational Principle for Branching Transport Networks*](paper2-variational/manuscript/) | Minimax game between wave and transport costs yields **α\* = 2.72** with zero free parameters |

---

## Quickstart

```bash
# 1. Clone and install dependencies
git clone https://github.com/rikymarche-ctrl/vascular-networks-theory.git
cd vascular-networks-theory
pip install -r requirements.txt

# 2. Reproduce everything (Python + LaTeX)
make all

# Or one paper at a time:
make paper1
make paper2
```

If `make` is not available (plain Windows), use the cross-platform script:

```bash
python reproduce.py           # both papers
python reproduce.py paper1    # Paper 1 only
python reproduce.py paper2    # Paper 2 only
```

> **Note:** LaTeX compilation requires a standard TeX distribution
> (TeX Live, MacTeX, or MiKTeX) with `pdflatex` and `bibtex`.
> Running only the Python step generates all numerical outputs and figures
> even without LaTeX.

---

## Repository Structure

```
vascular-networks-theory/
├── shared/
│   └── scripts/
│       ├── params.py                   # Single source of truth for all physical constants
│       ├── compute_paper_murray.py     # Paper 1 numerical computations
│       ├── compute_paper_variational.py# Paper 2 numerical computations
│       ├── compute_supplemental.py     # Supplemental computations
│       └── run_all.py                  # Master script (runs all three above)
├── paper1-murray/
│   ├── compute_paper_murray.py         # Entry point → delegates to shared/scripts/
│   ├── manuscript/
│   │   ├── main.tex                    # LaTeX source
│   │   ├── references.bib
│   │   ├── dynamic_variables.tex       # Auto-generated — do not edit by hand
│   │   └── figures/
│   └── supplements/
├── paper2-variational/
│   ├── compute_paper_variational.py    # Entry point → delegates to shared/scripts/
│   ├── manuscript/
│   └── supplements/
├── reproduce.py                        # Cross-platform build script
├── Makefile
├── requirements.txt
├── CITATION.cff
└── LICENSE
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
  year    = {2026},
  doi     = {10.5281/zenodo.18975093},
  url     = {https://zenodo.org/records/18975093}
}

@article{marchesi2026variational,
  author  = {Marchesi, Riccardo},
  title   = {A Unified Variational Principle for Branching Transport Networks},
  year    = {2026},
  doi     = {10.5281/zenodo.19027393},
  url     = {https://zenodo.org/records/19027393}
}
```

---

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
You are free to share and adapt it for any purpose, provided appropriate credit is given.
