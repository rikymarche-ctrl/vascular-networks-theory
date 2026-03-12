# A Unified Variational Principle for Branching Transport Networks: Wave Impedance, Viscous Flow, and Tissue Metabolism

This paper presents a network-level Lagrangian framework that unifies pulsatile wave-reflection penalties with steady transport-metabolic costs. By casting the morphological optimization as a minimax game between network architecture and environmental duty cycle, the model predicts the empirically observed cardiovascular branching exponent ($\alpha \approx 2.7$) without fitted parameters.

**Key results:**
- Deterministic prediction: $\alpha^* = 2.72$ for porcine coronary arteries.
- Proof of binary branching ($N=2$) as the unique dynamic stiffness maximizer.
- Minimax duty cycle $\eta^* \approx 0.83$ as a derived structural property.

## Reproduce

From the repository root:

```bash
# Python only (generates figures + dynamic_variables.tex)
python paper2-variational/compute_paper_variational.py

# Python + LaTeX (full PDF)
make paper2

# or
python reproduce.py paper2
```

## Outputs

| File | Description |
|---|---|
| `manuscript/dynamic_variables.tex` | Numerical constants and minimax eigenvalues |
| `manuscript/figures/fig1_kappa.pdf` | Emergence of the network stiffness ratio |
| `manuscript/figures/fig2_minimax.pdf` | Convergence to the minimax saddle point |
| `manuscript/A Unified Variational Principle for Branching Transport Networks.pdf` | Compiled manuscript |
| `supplements/Supplemental Material - Unified Variational Principle.pdf` | Coherent analysis and UQ tables |
