# The Dynamic Origin of Kleiber's Law
### and the Generalized Metabolic Scaling Theorem

**Riccardo Marchesi** — University of Pavia

---

## What this paper does

Kleiber's 3/4 metabolic scaling law has long been attributed to fractal geometry and viscous dissipation minimization. This paper inverts that narrative: **Kleiber's law is a signature of pulsatile wave physics, not steady-state geometry**.

The core result: by coupling local branching optimization to global allometry, the paper derives the exact generalized metabolic exponent β = dα/(2d+α), which strictly maps local transport microphysics to organismal-scale energetics. Wave-impedance matching in the proximal vasculature uniquely enforces α = 2, yielding β = 3/4 in three dimensions — and this bound is **dynamically protected**: no static viscous optimization can reproduce it.

The paper further:

- Derives the critical body mass M* for the wave-to-viscous transition, explaining the empirical shift to β ≈ 0.9 in small mammals with no free parameters
- Proves that the classical West–Brown–Enquist (WBE) derivation is structurally divergent under its own assumptions
- Validates the framework across **nine biological systems spanning five phyla** (vertebrate vasculature, insect tracheae, plant xylem, sponge canals)
- Organizes diverse networks into discrete **allometric universality classes** via the exact formula β(n, m, d)

---

## Key results at a glance

| Quantity | Value | Source |
|---|---|---|
| Generalized exponent β(n, m, d) | dα/(2d+α) | Theorem 1 |
| Kleiber fixed point β | 3/4 | Corollary A (α=2, d=3) |
| Allometric spectrum | β ∈ [3/4, 1] | Theorem 3 (Universal Bounds) |
| Minimax exponent α* | ≈ 2.77 | Proposition (physical bisection) |
| Minimax β | ≈ 0.947 | Corollary C |
| WBE inconsistency | proximal-dominance limit diverges | Theorem 7 |
| Biological systems validated | 9 (5 phyla) | Fig. 1 |

---

## Repository structure

```
paper3-kleiber/
├── manuscript/
│   ├── main.tex              ← main LaTeX source
│   ├── dynamic_variables.tex ← auto-generated numerical macros (do not edit manually)
│   ├── references.bib        ← bibliography
│   └── figures/
│       ├── fig1_spectrum.pdf  ← allometric universality classes β(n, m, d=3)
│       ├── fig2_transition.pdf← β(M) wave-to-viscous transition curve
│       └── fig3_minimax.pdf   ← minimax cost curves, α* saddle point
├── scripts/
│   └── compute.py            ← generates all numerical results and figures
├── build.ps1                 ← PowerShell build script (runs compute.py + LaTeX)
└── output/
    └── The Dynamic Origin of Kleibers Law.pdf
```

---

## Reproducing the results

### 1. Generate numerical variables and figures

```bash
cd scripts
python compute.py
```

Writes `manuscript/dynamic_variables.tex` and all three figures. No free parameters — all values derive from independent literature.

### 2. Compile the paper

```bash
cd manuscript
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Or use the PowerShell script from the `paper3-kleiber` directory:

```powershell
.\build.ps1
```

### Dependencies

- Python ≥ 3.8 with `numpy`, `scipy`, `matplotlib`
- LaTeX distribution with `amsmath`, `mathpazo`, `booktabs`, `hyperref`, `authblk`

---

## Theorem structure

```
Theorem 0    Optimal branching exponent α_t = (n + m) / 2
Corollary 1  Murray's law: α_t = 3  (n=4, m=2)
Corollary 2  Wall-dominated limit: α_t = (5+p)/2  (m = 1+p)
Theorem 1    Generalized Metabolic Scaling: β = dα / (2d + α)
Corollary A  Dynamic origin of Kleiber's law: β = 3/4  (α=2, d=3)
Proposition  Minimax equilibrium exponent α* ≈ 2.77 (physical bisection)
Corollary C  β(α*, d=3) ≈ 0.947
Theorem 2    Dimensional universality: β_d = d/(d+1) for any d
Theorem 3    Universal allometric bounds: β ∈ [3/4, 1]  (m ∈ [1,2], d=3)
Theorem 7    WBE geometric inconsistency
```

---

## License

Research & Code: CC BY 4.0
All manuscripts, scripts, numerical data, and figure-generation code: [github.com/rikymarche-ctrl/vascular-networks-theory](https://github.com/rikymarche-ctrl/vascular-networks-theory)
