# Beyond Murray's Law: Non-Universal Branching Exponents from Vessel-Wall Metabolic Costs

Murray's law predicts a universal cubic branching exponent (α = 3), yet empirical data
consistently show α ≈ 2.7–2.9. This paper demonstrates mathematically that incorporating
the sub-linear metabolic cost of the structural vessel wall (*h* ∝ *r*^0.77) into the
optimization objective strictly breaks this universality, providing a parameter-free
prediction that matches mammalian cardiovascular morphometry.

**Key result:** α\* ∈ [2.90, 2.94] with zero fitted parameters.

## Reproduce

From the repository root:

```bash
# Python only (generates figures + dynamic_variables.tex)
python paper1-murray/compute_paper_murray.py

# Python + LaTeX (full PDF)
make paper1

# or
python reproduce.py paper1
```

## Outputs

| File | Description |
|---|---|
| `manuscript/dynamic_variables.tex` | All numerical results as LaTeX macros |
| `manuscript/figures/fig_alpha_scale.pdf` | Branching exponent vs flow scale |
| `manuscript/Beyond Murray's Law.pdf` | Compiled manuscript |
