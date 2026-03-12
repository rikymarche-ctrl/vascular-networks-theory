"""
params.py — Single source of truth for all physical constants.

Every numerical value used across the compute scripts originates here.
No other file should define physical constants.

References:
  [1] Murray 1926, Proc. Natl. Acad. Sci.
  [2] Kassab 1993, Am. J. Physiol. 265, H350-H365
  [3] Rhodin 1967, J. Ultrastruct. Res. 18, 181-223
  [4] Nichols 2011, McDonald's Blood Flow in Arteries, 6th ed.
  [5] Wolinsky & Glagov 1967, Circ. Res. 20, 99-111
  [6] Caro 1978, Mechanics of the Circulation
  [7] Paul 1980, Handbook of Physiology, vol. 2, 201-235
  [8] Taber 1998, Biophys. J. 74, 109-114
"""

import numpy as np

# ── Fluid properties ─────────────────────────────────────────
MU       = 3.5e-3       # dynamic viscosity [Pa·s]                [6]
RHO      = 1060.0       # blood density [kg/m³]                   [4]

# ── Metabolic costs ──────────────────────────────────────────
b        = 1500.0        # blood volume metabolic cost [W/m³]     [1,8]
MW_LOW   = 5e3           # wall tissue cost, low  [W/m³]          [7]
MW_MID   = 20e3          # wall tissue cost, mid  [W/m³]          [7]
MW_HIGH  = 35e3          # wall tissue cost, high [W/m³]          [7]

# ── Wall thickness allometry ─────────────────────────────────
c0       = 0.041         # wall-thickness prefactor [m^{1-p}]     [3,4,5]
p        = 0.77          # wall-thickness exponent (dimensionless) [3,4,5]

# ── Derived cost coefficients ────────────────────────────────
# Cost per unit length = Viscous(r,Q) + Blood_volume(r) + Wall_tissue(r)
# Blood:  π r² × b           → coefficient B_blood = π b
# Wall:   2π r h(r) × m_w    → coefficient B_wall  = 2π m_w c0
#         (2πr is circumference, h = c0 r^p is wall thickness)
B_blood  = np.pi * b                              # [W/m³ × m²] = [W/m]
B_wall   = lambda mw: 2 * np.pi * mw * c0         # [W/m^{2+p}]

# ── Coronary artery morphometry ──────────────────────────────
r0_coronary  = 1.5e-3    # proximal LAD radius [m]                [2]
Q0_coronary  = 1.3e-6    # resting coronary flow [m³/s]           [4]
ell0_coronary = 15e-3    # proximal segment length [m]            [2]
G_coronary   = 11        # tree depth (Kassab orders 10→0)        [2]
beta_coronary = 2.0**(-1.0 / 2.90) # daughter/parent length ratio (~0.787)
# Note: Scripts use beta=1.0 (uniform lengths) as a baseline simplifying assumption.
# Sensitivity analysis confirms alpha* ≈ 2.7 is highly robust to variations in beta.

# ── Impedance physics ────────────────────────────────────────
# Z ∝ r^{-(5-p)/2} for h ∝ r^p
alpha_w  = (5 - p) / 2   # = 2.115 for p = 0.77

# ── Information-theoretic ────────────────────────────────────
K_fidelity = 50           # phenomenological SNR constant

# ── Compact cost-function A(Q) ───────────────────────────────
# Φ(r,Q,ℓ) = A(Q)/r⁴ ℓ + B_blood r² ℓ + B_wall r^{1+p} ℓ
A_of_Q = lambda Q: 8 * MU * Q**2 / np.pi

# ── Wave mechanics (transfer-matrix, Supplemental S1) ────────
E_elastic          = 0.4e6   # Pa — coronary elastic modulus       [Nichols 2011]
f0_cardiac         = 1.2     # Hz — fundamental cardiac frequency
ell_factor         = 10      # ℓ_g = ell_factor × r_g (self-similar morphometry)

# ── Power budget (Supplemental S3) ───────────────────────────
CPO_watts          = 1.0     # W  — cardiac power output at rest
pulsatile_fraction = 0.15    # pulsatile fraction of CPO
coronary_fraction  = 0.05    # coronary fraction of cardiac output

# ── Pulmonary artery parameters ──────────────────────────────
p_pulmonary = 0.60        # thinner walls                        [Huang 1996]

# ── Convenience ──────────────────────────────────────────────
def wall_alpha_limit(p_val=p):
    """Single-term limit: α = (5+p)/2"""
    return (5 + p_val) / 2

def murray_alpha():
    """Murray's cubic law."""
    return 3.0

def classification_alpha(gamma):
    """Theorem 3: α = (4+γ)/2 for two-term cost ~ r^γ."""
    return (4 + gamma) / 2
