# The Incommensurability Principle in Biological Transport

### A No-Go Proposition, Metabolic Gauge Invariance, and the Architectural Origin of the Minimax

**Riccardo Marchesi** — University of Pavia

**Status:** Major revision in progress for *Transport Phenomena* (De Gruyter Brill). Internal audit + revision: 2026-06-20.

---

## What this paper does

This paper establishes the theoretical foundation for the network-level Lagrangian framework introduced in Paper II, arguing that its structure is not a modelling convenience but a structural consequence of the physical incommensurability of biological cost functions.

The framework successfully predicts:

- **Critical Womersley number:** $\mathrm{Wo}_c = \sqrt{3} \approx 1.732$ (0.5% from exact 1.740)
- **Allometric transition:** $M^* \approx 0.8\,$g separates viscous/wave regimes
- **Branching exponent:** $\alpha^* \approx 2.72$ from minimax balance (vs empirical 2.70±0.20)
- **Ontogenetic trajectory:** Testable prediction of $\alpha \approx 3.0 \to 2.7$ transition in embryonic development

Three core results are proved (one proposition, two theorems):

1. **No-Go Proposition** (Biologically Implausible Local Optimization): Any Lagrangian combining an extensive transport-metabolic penalty (watts) with a dimensionless wave-reflection penalty at a single junction requires a coupling parameter $\mu$ varying by $10^2$--$10^3$ across the vascular hierarchy — making a universal branching exponent $\alpha^*$ biologically implausible under local optimization (symmetric rules).

   - **Information-accounting corollary:** regulating a generation-dependent profile requires only the generation index ($\sim \log_2 G$ bits); the obstruction is that this index is *non-local* — it cannot be reconstructed from purely local hemodynamic ratios — not a Shannon channel-capacity bound.
2. **Metabolic Gauge Invariance**: The unique dimensionless network-level transport penalty consistent with scale invariance, positivity, and the linear thermodynamics of entropy production is the fractional excess cost $\mathcal{C}_\mathrm{transport}^\mathrm{net} = (\Phi_\mathrm{net} - \Phi_\mathrm{opt})/\Phi_\mathrm{opt}$. All alternatives, including logarithmic measures, violate thermodynamic linearity.

   - **Extensivity caveat:** Linearity follows from extensivity OR Onsager near-equilibrium regime OR physiological near-optimality.
   - **Gauge terminology:** Global metabolic scaling symmetry, distinct from field theory gauge symmetries (U(1), SU(2)).
3. **Architectural Invariance**: The minimax saddle weight $\eta^*$ is an exact invariant of the network's allometric class, strictly orthogonal to absolute metabolic scales — resolving the ontogenetic paradox of stable vascular geometry throughout growth.

   - **Cross-Class Transition:** Body mass crossing $M^* \approx 0.8\,$g triggers viscous→wave transition.

A corollary recovers Papers I and II as the degenerate boundary cases $\eta \to 0$ and $\eta \to 1$ of the unified principle.

---

## Key results at a glance

| Quantity                           | Value                                                         | Source                               |
| ---------------------------------- | ------------------------------------------------------------- | ------------------------------------ |
| Critical Womersley number          | $\mathrm{Wo}_c = \sqrt{3} \approx 1.732$                    | Kinematic Matching Criterion         |
| Allometric transition mass         | $M^* \approx 0.8\,$g                                        | Theorem 3 (Cross-Class Transition)   |
| Symmetric minimax prediction       | $\alpha^* \approx 2.629$ (Womersley-rigorous)               | Section 7 (Model Hierarchy)          |
| Large-mammal empirical             | $\alpha^* \approx 2.72$                                     | Kassab et al. 1993                   |
| $\mu$ variation (local coupling) | $10^2$--$10^3$                                            | Proposition 1 (No-Go)                |
| Information-accounting argument    | generation index, $\sim \log_2 G$ bits (non-local)         | Corollary (No-Go)                    |
| Admissible functional under stated assumptions       | $(\Phi_\mathrm{net} - \Phi_\mathrm{opt})/\Phi_\mathrm{opt}$ | Theorem 2 (Gauge)                    |
| Minimax saddle weight$\eta^*$    | $\approx 0.814$                                             | Theorem 3                            |
| $\eta^*$ dependence on body mass | none (exact invariant)                                        | Theorem 3 (Architectural Invariance) |
| Ontogenetic prediction             | $\alpha \approx 3.0 \to 2.7$ (E10→Adult)                   | Section 8 (Conclusion)               |
| Papers I, II recovered as          | $\eta \to 0$, $\eta \to 1$ limits                         | Corollary                            |

---

## New theoretical content (post peer review)

### Renormalization Group Interpretation

The minimax attractor $\alpha^*$ functions as a **renormalization group fixed point** of hierarchical scale transformations. The branching exponent is a critical exponent; metabolic parameters are irrelevant operators. Finite-size scaling prediction: $\alpha^*(G) = \alpha^*_\infty + c/G + \mathcal{O}(G^{-2})$.

### Model Hierarchy and the Residual Gap

Three-tier physical hierarchy explains the gap between theory and experiment:

1. **Static transport optimization** (Paper I): $\alpha_t \approx 2.92$
2. **Symmetric minimax** (this work): $\alpha^* \approx 2.629$ (Womersley-rigorous)
3. **Large-mammal empirical** (with heterogeneities): $\alpha^* \approx 2.72$

The Womersley correction ($F_{10}$) accounts for $\Delta\alpha \approx 0.02$; remaining shift ($\Delta\alpha \approx 0.07$) reflects scale-dependent morphometric heterogeneities (asymmetry, taper, flow partitioning).

### Critical Falsifiable Test: Ontogenetic Phase Transition

**Prediction:** Embryonic vasculature should exhibit Murray's viscous scaling ($\alpha \approx 3.0$) at early stages ($M < M^* \approx 0.8\,$g), then transition to wave-influenced minimax ($\alpha^* \approx 2.7$) as body mass crosses the Womersley threshold.

**Falsification criterion:** If embryonic vasculature exhibits $\alpha \approx 2.7$ from the earliest stages of angiogenesis (when wave reflections are physically absent), the theory is falsified.

**Experimental protocol:** Time-resolved micro-CT of staged mouse/zebrafish embryos measuring $\alpha(t)$ vs developmental time.

### The Hummingbird Test

**Heart rate overrides body mass:**

- Mouse (20g, 600bpm): $\mathrm{Wo} = 2.5 \to \alpha^* \approx 2.60$ (viscous)
- Hummingbird (4g, 1000bpm): $\mathrm{Wo} = 1.76 \to \alpha^* \approx 2.6$ (transition)
- Human (70kg, 70bpm): $\mathrm{Wo} = 17.9 \to \alpha^* \approx 2.72$ (minimax)

### Retinal Paradox as Direct Experimental Proof

The simultaneous divergence of two attractors in retinal vasculature:

- **Diameters:** $\alpha \approx 2.0$ (2D wave-matching, $\mathrm{Wo}_c = \sqrt{6}$)
- **Bifurcation angles:** $\theta \approx 75°$ (3D Murray equilibrium)

This dual-attractor state is **impossible under local optimization** (would converge to single attractor) — direct experimental support for the No-Go Theorem.

---

## Repository structure

```
paper4-incommensurability/
├── manuscript/
│   ├── main.tex              ← main LaTeX source (uses \input for modular sections)
│   ├── sections/             ← modular section files (13 files)
│   │   ├── Section-01-Introduction.tex
│   │   ├── Section-02-Scaling-Conflict.tex
│   │   ├── Section-03-Gauge-Invariance.tex
│   │   ├── Section-04-Architectural-Invariance.tex
│   │   ├── Section-05-Single-Mechanism-Limits.tex
│   │   ├── Section-06-Architectural-Transition.tex
│   │   ├── Section-07-Discussion.tex
│   │   ├── Section-08-Conclusion.tex
│   │   ├── Section-09-Retinal-Paradox.tex
│   │   ├── Section-10-Data-Availability.tex
│   │   ├── Appendix-A-Fragility.tex
│   │   ├── Appendix-B-Transfer-Matrix.tex
│   │   └── Appendix-C-Sensitivity.tex
│   ├── dynamic_variables.tex ← auto-generated numerical macros (do not edit manually)
│   ├── references.bib        ← BibTeX bibliography
│   └── figures/
│       └── fig_phase_diagram.tex ← TikZ phase diagram figure
├── supplements/
│   ├── supplemental.tex      ← supplemental material (Appendices A-F + empirical data)
│   ├── references_supp.bib   ← supplemental bibliography
│   └── figures/
│       └── womersley_verification.png
├── scripts/
│   ├── main/
│   │   └── compute.py        ← generates all numerical results (with dynamic pruning)
│   ├── verification/
│   │   ├── generate_womersley_figure.py
│   │   ├── verify_heterogeneity.py
│   │   └── verify_shielding.py
│   └── utility/
│       ├── wrap_project_tex.py
│       └── check_orphaned_macros.py
├── build.ps1                 ← PowerShell build script (runs compute.py + LaTeX)
└── output/
    ├── The Incommensurability Principle Clean.pdf
    └── Supplemental Material Clean.pdf
```

---

## Reproducing the results

### 1. Generate numerical variables

```bash
cd scripts/main
python compute.py
```

Writes `manuscript/dynamic_variables.tex` with automatic dynamic pruning (removes orphaned macros). All target-exponent predictions use independently constrained inputs; no parameter is fitted to the target morphometric exponent.

### 2. Compile the paper

**Recommended:** Use the automated build script from the `paper4-incommensurability` directory:

```powershell
.\build.ps1
```

This script automatically:

1. Runs `compute.py` to generate numerical variables
2. Generates verification figures
3. Compiles main manuscript (3 pdfLaTeX passes + BibTeX)
4. Compiles supplemental material (3 pdfLaTeX passes + BibTeX)
5. Copies PDFs to `output/` and `../Output_PDFs/`

**Manual compilation** (if needed):

```bash
cd manuscript
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

**Build detail:** The manuscript uses a modular structure with 13 separate section files in `manuscript/sections/`. The build process handles these automatically via `\input{}` commands.

### Dependencies

- Python ≥ 3.8 with `numpy`, `scipy`
- LaTeX distribution with `amsmath`, `mathpazo`, `booktabs`, `hyperref`, `authblk`, `tikz`

---

## Modular manuscript structure

The manuscript follows a modular organization pattern for easier editing and maintenance:

- **Main file** (`manuscript/main.tex`): Contains only document preamble, `\input{}` commands, and bibliography
- **Section files** (`manuscript/sections/*.tex`): Each section/appendix in a separate file
- **Benefits**:
  - Edit individual sections without touching the full document
  - Better version control (smaller diffs)
  - Easier collaboration and review
  - Clear logical organization

**Editing workflow:**

1. Edit the relevant section file in `manuscript/sections/`
2. Run `.\build.ps1` to recompile
3. Changes are automatically integrated into the final PDF

**Supplemental material:**

- Organized as Appendices A-F (properly numbered)
- Empirical validation sections follow appendices
- Bibliography at the end (correct LaTeX structure)

---

## Theoretical structure

```
Theorem 0         Kinematic Matching Criterion (Woc = √3)
                  Epistemological status: physical matching criterion, with exact-Bessel validation
                  Cellular mechanobiological feedback loop to be developed in Paper V

Proposition 1     No-Go: Local incommensurable optimization is biologically implausible (symmetric rules)
  Corollary 1.1   Information-accounting: generation index (~log₂G bits) is non-locally inaccessible

Theorem 2         Metabolic Gauge Invariance: Admissible network functional under stated assumptions
  Remark 2.1      Extensivity caveat: linearity conditional on extensivity/Onsager/near-optimality
  Remark 2.2      Gauge terminology: global scaling symmetry vs field theory gauges

Theorem 3         Architectural Invariance: η* is an allometric class invariant
  Theorem 3.2     Cross-Class Transition: M* ≈ 0.8g separates viscous/wave regimes
  Corollary 3.1   Topological Rigidity: α* independent of metabolic parameters

Corollary         Papers I and II recovered as η→0, η→1 degenerate limits
```


---

## License

Research & Code: CC BY 4.0
All manuscripts, scripts, numerical data, and figure-generation code: [github.com/rikymarche-ctrl/vascular-networks-theory](https://github.com/rikymarche-ctrl/vascular-networks-theory)
