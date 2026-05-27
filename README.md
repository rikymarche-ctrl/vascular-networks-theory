# Branching Transport Networks — Research Monorepo

A unified theoretical framework for the optimal geometry of biological
branching networks. The series derives branching exponents from first
principles — no free parameters — across vascular, pulmonary, neural,
plant, and invertebrate systems.

**Author:** Riccardo Marchesi — University of Pavia

---

## Papers

### Paper I — *Beyond Murray's Law*
> `paper1-murray/` · [arXiv:2603.13687](https://arxiv.org/abs/2603.13687)

Extends Murray's cubic law by incorporating vessel-wall tissue metabolic
cost. Derives α*(Q) ∈ [2.885, 3.000] via seven rigorous theorems
(existence, uniqueness, classification, angle bounds, N=2 selection).

### Paper II — *A Unified Variational Principle for Branching Transport Networks*
> `paper2-variational/` · [arXiv:2603.14691](https://arxiv.org/abs/2603.14691)

Network-level two-level Lagrangian combining transport cost and wave
impedance. Derives α* = 2.720 with zero free parameters; η* = 0.833 as
the minimax duty cycle. Includes supplemental material.

### Paper III — *The Dynamic Origin of Kleiber's Law*
> `paper3-kleiber/` · [arXiv:2604.10476](https://arxiv.org/abs/2604.10476)

Proves that α_t = (n+m)/2 has genuine predictive power across nine
biological systems (vascular, pulmonary, neural, plant, invertebrate).
Shows that Kleiber's 3/4 scaling emerges from the wave-impedance attractor
α_w = 2, not from Murray's law. Derives a general allometric equation of
state β(α, d) = dα/(2d+α).

### Paper IV — *The Incommensurability Principle in Biological Transport*
> `paper4-incommensurability/` · [arXiv:2605.03219](https://arxiv.org/abs/2605.03219)

Proves that universal vascular branching exponents cannot emerge from local
optimization: any junction-level coupling of incommensurable costs requires
scale-dependent fine-tuning varying by O(10²–10³) across the hierarchy.
Establishes Topological Rigidity (optimal exponent depends only on dimensionless
structural parameters) and Architectural Invariance (η* is an exact allometric
class invariant), with a dual-threshold Womersley framework.

---

## Key Results

| Paper | Core result | Status |
|-------|-------------|--------|
| I | α* ∈ [2.885, 3.000] from vessel-wall cost | Published |
| II | α* = 2.720, η* = 0.833, zero free params | Published |
| III | β(α,d) = dα/(2d+α); Kleiber from α_w = 2 | Published |
| IV | Incommensurability → unique minimax structure | Published |

---

## Repository Structure

```text
branching-papers/
├── build.ps1                       # Master build pipeline (delegates to per-paper scripts)
├── tools/
│   ├── format.py                   # Unified LaTeX formatter (sections + 80-char wrap)
│   └── template_manuscript.tex     # Standard preamble template for new papers
├── shared/
│   └── scripts/
│       └── params.py               # Physical constants (single source of truth)
│
├── paper[N]-[id]/                  # One directory per paper — identical template
│   ├── build.ps1                   # Local build script
│   ├── manuscript/
│   │   ├── main.tex                # LaTeX source
│   │   ├── references.bib
│   │   ├── dynamic_variables.tex   # Auto-generated numbers — do not edit by hand
│   │   └── figures/
│   ├── scripts/
│   │   └── compute.py              # Deterministic numerical pipeline
│   ├── output/                     # Compiled PDFs and figure assets
│   └── [supplements/ | data/ | …]  # Optional — paper-specific appendices or datasets
│
└── Output_PDFs/                    # Mirror of all compiled PDFs across papers
```

---

## Build

**All papers** (Python → LaTeX → PDFs in `Output_PDFs/`):
```powershell
.\build.ps1
```

**Single paper:**
```powershell
.\build.ps1 -Target paper1
.\build.ps1 -Target paper3
```

**Individual paper** (from its own folder):
```powershell
cd paper3-kleiber
.\build.ps1
```

**Reformat LaTeX source** (section headers + 80-char wrap):
```powershell
python tools/format.py          # all papers
python tools/format.py 3 4      # papers 3 and 4
```

**Requirements:** Python ≥ 3.8 with `numpy`, `scipy`, `matplotlib`;
MiKTeX or TeX Live.

---

## License

Research & Code: CC BY 4.0
[github.com/rikymarche-ctrl/vascular-networks-theory](https://github.com/rikymarche-ctrl/vascular-networks-theory)
