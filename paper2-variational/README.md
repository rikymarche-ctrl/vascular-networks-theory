# A Unified Variational Principle for Branching Transport Networks
### Wave Impedance, Viscous Flow, and Tissue Metabolism

**Riccardo Marchesi** — University of Pavia

---

## What this paper does

This paper presents a network-level Lagrangian framework that unifies pulsatile wave-reflection penalties with steady transport-metabolic costs. By casting the morphological optimization as a minimax game between network architecture and environmental duty cycle, the model predicts the empirically observed cardiovascular branching exponent ($\alpha \approx 2.7$) without fitted parameters.

The core result: the apparent conflict between impedance matching ($\alpha \approx 2.1$) and minimum dissipation ($\alpha \approx 3$) is rigorously resolved by formulating the optimization at the network level and incorporating the biological cost of the conduit wall. The arterial tree sits precisely at the Pareto-optimal frontier where the marginal return on signal integrity equals the marginal cost of transport maintenance.

The paper demonstrates:
- Deterministic prediction of $\alpha^* = 2.72$ for porcine coronary arteries.
- A proof that binary branching ($N=2$) is the unique dynamic stiffness maximizer.
- Derivation of the minimax duty cycle $\eta^* \approx 0.833$ as a structural invariant rather than a variable parameter.

---

## Key results at a glance

| Quantity | Value | Source |
|---|---|---|
| Predicted optimal exponent $\alpha^*$ | 2.72 | Minimax game |
| Empirical value $\alpha_{\exp}$ | $2.70 \pm 0.20$ | Kassab 1993 |
| Minimax duty cycle $\eta^*$ | 0.833 | Network equilibrium |
| Peak transport exponent $\alpha_t$ | 2.90 | Single-vessel static optimum |
| Network wave cost $\mathcal{C}_\mathrm{wave}^\mathrm{net}$ at $\alpha^*$ | 6.3% | Predicted marginal penalty |

---

## Repository structure

```text
paper2-variational/
├── manuscript/
│   ├── main.tex              ← main LaTeX source
│   ├── dynamic_variables.tex ← auto-generated numerical macros
│   ├── references.bib        ← bibliography
│   └── figures/
├── supplements/
│   ├── supplemental.tex      ← supplemental materials and proofs
│   ├── dynamic_variables_supplemental.tex
│   └── figures/
├── scripts/
│   ├── compute.py            ← generates variables for main text
│   └── compute_supp.py       ← generates variables for supplements
├── build.ps1                 ← PowerShell build script
└── output/
    ├── A Unified Variational Principle for Branching Transport Networks.pdf
    └── Supplemental Material - Unified Variational Principle.pdf
```

---

## Reproducing the results

### 1. Generate numerical variables and figures

```bash
cd scripts
python compute.py
python compute_supp.py
```

This writes the necessary `dynamic_variables.tex` and `.pdf` figures in their respective directories.

### 2. Compile the paper

Or, simply use the unified PowerShell script (Windows) from the `paper2-variational` directory:

```powershell
.\build.ps1
```

### Dependencies

- Python ≥ 3.8 with `numpy`, `scipy`, `matplotlib`
- LaTeX distribution with standard packages (`amsmath`, `hyperref`, etc.)

---

## License

Research & Code: CC BY 4.0
All manuscripts, scripts, numerical data, and figure-generation code: [github.com/rikymarche-ctrl/vascular-networks-theory](https://github.com/rikymarche-ctrl/vascular-networks-theory)
