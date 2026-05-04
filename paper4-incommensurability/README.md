# The Incommensurability Principle in Biological Transport
### A No-Go Theorem, Metabolic Gauge Invariance, and the Architectural Origin of the Minimax

**Riccardo Marchesi** — University of Pavia

---

## What this paper does

This paper establishes the rigorous theoretical foundation for the network-level Lagrangian framework introduced in Paper II, proving that its structure is not a modelling convenience but a **mathematical necessity** imposed by the physical incommensurability of biological cost functions.

Three theorems are proved:

1. **No-Go Theorem**: any Lagrangian combining an extensive transport-metabolic penalty (watts) with a dimensionless wave-reflection penalty at a single junction requires a coupling parameter $\mu$ varying by $10^2$--$10^3$ across the vascular hierarchy — making a universal branching exponent $\alpha^*$ impossible under local optimization.

2. **Metabolic Gauge Invariance**: the unique dimensionless network-level transport penalty consistent with scale invariance, positivity, and the linear thermodynamics of entropy production is the fractional excess cost $\mathcal{C}_\mathrm{transport}^\mathrm{net} = (\Phi_\mathrm{net} - \Phi_\mathrm{opt})/\Phi_\mathrm{opt}$. All alternatives, including logarithmic measures, violate thermodynamic linearity.

3. **Architectural Invariance**: the minimax duty cycle $\eta^*$ is an exact invariant of the network's allometric class, strictly orthogonal to absolute metabolic scales — resolving the ontogenetic paradox of stable vascular geometry throughout growth.

A corollary recovers Papers I and II as the degenerate boundary cases $\eta \to 0$ and $\eta \to 1$ of the unified principle.

---

## Key results at a glance

| Quantity | Value | Source |
|---|---|---|
| $\mu$ variation across vascular hierarchy | $10^2$--$10^3$ | Theorem 1 (No-Go) |
| Unique admissible functional | $(\Phi_\mathrm{net} - \Phi_\mathrm{opt})/\Phi_\mathrm{opt}$ | Theorem 2 (Gauge) |
| Minimax duty cycle $\eta^*$ | $\approx 0.833$ | Theorem 3 |
| $\eta^*$ dependence on body mass | none (exact invariant) | Corollary |
| Papers I, II recovered as | $\eta \to 0$, $\eta \to 1$ limits | Corollary |

---

## Repository structure

```
paper4-incommensurability/
├── manuscript/
│   ├── main.tex              ← main LaTeX source
│   ├── dynamic_variables.tex ← auto-generated numerical macros (do not edit manually)
│   ├── references.bib        ← BibTeX bibliography
│   └── figures/
│       └── fig_phase_diagram.tex ← TikZ phase diagram figure
├── scripts/
│   └── compute.py            ← generates all numerical results
├── build.ps1                 ← PowerShell build script (runs compute.py + LaTeX)
└── output/
    └── The Incommensurability Principle in Biological Transport.pdf
```

---

## Reproducing the results

### 1. Generate numerical variables

```bash
cd scripts
python compute.py
```

Writes `manuscript/dynamic_variables.tex`. All values derived from first principles — no free parameters.

### 2. Compile the paper

```bash
cd manuscript
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Or use the PowerShell script from the `paper4-incommensurability` directory:

```powershell
.\build.ps1
```

### Dependencies

- Python ≥ 3.8 with `numpy`, `scipy`
- LaTeX distribution with `amsmath`, `mathpazo`, `booktabs`, `hyperref`, `authblk`, `tikz`

---

## Theorem structure

```
Theorem 1    No-Go Theorem: local incommensurable optimization has no solution
Theorem 2    Metabolic Gauge Invariance: unique admissible network functional
Theorem 3    Architectural Invariance: η* is an exact allometric class invariant
Corollary    Papers I and II recovered as η→0, η→1 degenerate limits
```

---

## License

Research & Code: CC BY 4.0
All manuscripts, scripts, numerical data, and figure-generation code: [github.com/rikymarche-ctrl/vascular-networks-theory](https://github.com/rikymarche-ctrl/vascular-networks-theory)
