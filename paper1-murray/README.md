# Beyond Murray's Law
### Non-Universal Branching Exponents from Vessel-Wall Metabolic Costs

**Riccardo Marchesi** — University of Pavia  
Preprint: [arXiv](https://arxiv.org/abs/2603.13687)

---

## What this paper does

Murray's cubic law (α = 3) has been the standard prediction for vascular branching geometry since 1926. This paper shows that the universality of α = 3 is not a property of biological networks — it is an artifact of using an incomplete cost function.

The core result: adding the vessel-wall metabolic cost $h(r) = c_0 r^p$ ($p \approx 0.77$, measured histologically) to Murray's two-term Lagrangian renders the cost function **inhomogeneous with incommensurate scaling exponents**. By Cauchy's functional equation, this inhomogeneity is both necessary and sufficient to destroy universality. Murray's law is thereby identified as a singular degeneracy of the cost-function family, not a general biological principle.

The paper then derives, without fitting any morphometric data:

- a general classification of any three-term additive cost function showing that non-universality is the generic state for distinct maintenance exponents (Corollary 7)
- strict bounds **(5+p)/2 < α*(Q) < 3** for all flows Q (Theorem 4)
- a static thermodynamic attractor at **α* ∈ [2.90, 2.94]** for porcine coronary arteries
- tight bounds on bifurcation angles **74.9° < 2θ* < 80.2°**, confirmed by 3D morphometry (Theorem 10)
- a proof that purely static optimization **cannot** select binary branching (N=2) — the metabolic cost ratio $K_2/K_1 \approx 0.19$ falls three orders of magnitude short of the threshold $\tau(p) \approx 10^3$ (Corollary 14)

This last result — the **Static Insufficiency Corollary** — is the
deepest finding: it constitutes a mathematical proof that binary
branching cannot be explained by steady-flow metabolic optimisation
alone. The most physically motivated candidate for the missing
constraint is pulsatile wave-reflection (impedance mismatch) at
high-degree nodes — a mechanism whose necessity is mandated by the
static analysis, though its formal treatment lies at the network
level beyond the scope of this work.

---

## Key results at a glance

| Quantity | Value | Source |
|---|---|---|
| Wall exponent $p$ | 0.77 | Histology (Rhodin 1967, Kassab 1993) |
| Theoretical lower bound $(5+p)/2$ | 2.885 | Theorem 4 |
| Static attractor $\alpha^*$ | 2.90 – 2.94 | Numerical (Table 2) |
| Empirical value $\alpha_{\exp}$ | 2.70 ± 0.20 | Kassab 1993 |
| Bifurcation angle bounds | 74.9° – 80.2° | Theorem 10 |
| $K_2/K_1$ (coronary, actual) | 0.19 | Corollary 14 |
| $\tau(p)$ threshold for N=2 | $10^3$ | Proposition 13 |

The gap between $\alpha^* \approx 2.90$ and $\alpha_{\exp} \approx 2.70$ is not a model failure. It is a diagnostic signature: any exponent significantly below 2.90 mathematically signals the competitive action of pulsatile wave dynamics.

---

## Repository structure

```
paper1-murray/
├── manuscript/
│   ├── main.tex              ← main LaTeX source
│   ├── dynamic_variables.tex ← auto-generated numerical macros (do not edit manually)
│   ├── references.bib        ← bibliography
│   └── figures/
│       └── fig_alpha_scale.pdf
├── scripts/
│   └── compute.py            ← generates all numerical results and fig_alpha_scale.pdf
├── build.ps1                 ← PowerShell build script (runs compute.py + LaTeX)
└── output/
    └── Beyond Murray's Law.pdf
```

---

## Reproducing the results

### 1. Generate numerical variables and figure

```bash
cd scripts
python compute.py
```

This writes `manuscript/dynamic_variables.tex` and `manuscript/figures/fig_alpha_scale.pdf`. All parameters are drawn from independent literature — no morphometric fitting.

### 2. Compile the paper

```bash
cd manuscript
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Or use the PowerShell script (Windows) from the `paper1-murray` directory:

```powershell
.\build.ps1
```

### Dependencies

- Python ≥ 3.8 with `numpy`, `scipy`, `matplotlib`
- LaTeX distribution with `amsmath`, `mathpazo`, `booktabs`, `hyperref`, `authblk`

---

## Physical parameters (no free parameters)

| Symbol | Value | Unit | Source |
|---|---|---|---|
| μ (viscosity) | 3.5 | mPa·s | Caro 1978 |
| b (blood metabolism) | 1500 | W/m³ | Murray 1926, Taber 1998 |
| $m_w$ (wall metabolism) | 5 – 35 | kW/m³ | Paul 1980 |
| $c_0$ (wall thickness prefactor) | 0.041 | m^(1−p) | Rhodin 1967 |
| $p$ (wall thickness exponent) | 0.77 | — | Rhodin 1967, Kassab 1993 |
| $Q_0$ (coronary flow) | 1.3 | mL/s | Kassab 1993 |
| $r_0$ (coronary radius) | 1.5 | mm | Kassab 1993 |

No parameter is fitted to branching exponent data.

---

## Theorem structure

```
Theorem 1    Existence and uniqueness of optimal radius r*(Q)
Lemma 2      Power law ⟺ Murray's branching law (Cauchy)
Theorem 3    Single-term classification: α = (4+γ)/2
Theorem 4    Strict bounds: (5+p)/2 < α*(Q) < 3
Corollary 5  Bounds hold for all flow asymmetries f ∈ (0,1)
Corollary 6  Murray's law is the unique universal member of this family
Corollary 7  General three-term incommensurability (for distinct m, k)
Proposition 8 Physical determination of Bennett's parameter m = 1+p
Lemma 9      Optimal cost Φ*(r) in closed form
Theorem 10   Bifurcation angles and bounds (Angle-Reflection tradeoff)
Theorem 11   Topological bounding: N* finite and small
Corollary 12 Shallow gradient explains occasional trifurcations
Proposition 13 N=2 is unique iff K₂/K₁ > τ(p)
Corollary 14 [Static Insufficiency] K₂/K₁ ≈ 0.19 ≪ τ(p) ≈ 10³:
             static optimization cannot select N=2
Remark       Capillary anastomoses predicted by κ → 1 crossover
```

---

## License

Research & Code: CC BY 4.0  
All manuscripts, scripts, numerical data, and figure-generation code: [github.com/rikymarche-ctrl/vascular-networks-theory](https://github.com/rikymarche-ctrl/vascular-networks-theory)
