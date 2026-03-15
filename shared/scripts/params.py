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
MU       = 3.5e-3        # dynamic viscosity [Pa·s]                 [6]
RHO      = 1060.0        # blood density [kg/m³]                    [4]

# ── Metabolic costs ──────────────────────────────────────────
b        = 1500.0        # blood volume metabolic cost [W/m³]      [1,8]
MW_LOW   = 5e3           # wall tissue cost, low  [W/m³]           [7]
MW_MID   = 20e3          # wall tissue cost, mid  [W/m³]           [7]
MW_HIGH  = 35e3          # wall tissue cost, high [W/m³]           [7]

# ── Wall thickness allometry ─────────────────────────────────
c0       = 0.041         # wall-thickness prefactor [m^{1-p}]      [3,4,5]
p        = 0.77          # wall-thickness exponent (dimensionless) [3,4,5]

# ── Derived cost coefficients ────────────────────────────────
B_blood  = np.pi * b     # [W/m³ × m²] = [W/m]

def B_wall(mw):
    """
    Calculate the wall-tissue metabolic cost coefficient per unit length.
    Coefficient B_wall is defined such that Wall_cost = B_wall * r^{1+p}.
    B_wall = 2π * m_w * c0 (units: W/m^{2+p}).
    """
    return 2 * np.pi * mw * c0

# ── Coronary artery morphometry ──────────────────────────────
r0_coronary  = 1.5e-3              # proximal LAD radius [m]                 [2]
Q0_coronary  = 1.3e-6              # resting coronary flow [m³/s]            [4]
Q0_peak      = 4.0 * Q0_coronary   # peak systolic coronary flow [m³/s] (~5.2e-6)
ell0_coronary = 15e-3              # proximal segment length [m]             [2]
G_coronary   = 11                  # tree depth (Kassab orders 10→0)         [2]
beta_coronary = 2.0**(-1.0 / 2.90) # daughter/parent length ratio (~0.787)
# Note: Scripts use beta=1.0 (uniform lengths) as a baseline simplifying assumption.
# Sensitivity analysis confirms alpha* ≈ 2.7 is highly robust to variations in beta.

# ── Impedance physics ────────────────────────────────────────
# Z ∝ r^{-(5-p)/2} for h ∝ r^p
alpha_w  = (5 - p) / 2   # = 2.115 for p = 0.77
# Empirical elastic stiffening (distal stiffening, Kassab/Nichols)
ke       = 0.23           # Z ∝ r^{-(5-p+ke)/2} when E ∝ r^{-ke}

# ── Compact cost-function A(Q) ───────────────────────────────
# The total metabolic cost per unit length is defined as:
# Φ(r,Q,ℓ) = ℓ ⋅ [ A(Q)r⁻⁴ + B_blood r² + B_wall r^{1+p} ]
def A_of_Q(Q):
    """
    Calculate the viscous dissipation coefficient A(Q).
    A(Q) = 8 * μ * Q² / π (units: W·m³).
    """
    return 8 * MU * Q**2 / np.pi

# ── Wave mechanics (transfer-matrix, Supplemental S1) ────────
E_elastic          = 0.4e6   # Pa — coronary elastic modulus       [Nichols 2011]
f0_cardiac         = 1.2     # Hz — fundamental cardiac frequency
ell_factor         = 10      # ℓ_g = ell_factor × r_g (self-similar morphometry)

# ── Power budget (Supplemental S2) ───────────────────────────
CPO_watts          = 1.0     # W  — cardiac power output at rest
pulsatile_fraction = 0.15    # pulsatile fraction of CPO
coronary_fraction  = 0.05    # coronary fraction of cardiac output

# ── Cross-system validation (literature parameters) ──────────
p_pulmonary = 0.60                          # thinner walls [Huang 1996]
alpha_w_pulmonary = (5 - p_pulmonary) / 2   # ≈ 2.200
G_pulmonary = 15                            # tree depth
alpha_exp_pulmonary = 2.75                  # morphometric central value
alpha_exp_err_pulmonary = 0.15              # 1σ uncertainty
alpha_exp_pulmonary_low = 2.70              # range lower bound
alpha_exp_pulmonary_high = 2.80             # range upper bound

alpha_w_aortic_low  = 2.10                  # p ∈ [0.50, 0.80] lower
alpha_w_aortic_high = 2.25                  # p ∈ [0.50, 0.80] upper
G_aortic = 15
alpha_exp_aortic = 2.5                      # morphometric central value
alpha_exp_err_aortic = 0.3                  # 1σ uncertainty

alpha_w_neural = 1.50                       # electrotonic, Rall 1959
G_neural = 8                                # indicative tree depth
alpha_star_neural_indicative = 2.3          # indicative (Cuntz MST estimate, footnote ‡)
alpha_exp_neural = 2.4                      # Cuntz 2010 (MST metric)
sigma_neural_indicative = 0.3               # indicative, no published σ

G_airways = 23                              # Weibel 1963 bronchial tree
alpha_exp_airways_low  = 2.8                # experimental range
alpha_exp_airways_high = 3.0

# ── Minimax band parameter ranges ────────────────────────────
p_low  = 0.72             # physiological range for p
p_high = 0.82
G_low  = 9                # tree depth range
G_high = 13
alpha_t_low_band  = 2.80  # transport optimum range
alpha_t_high_band = 3.00

# ── Sensitivity table physical baselines ─────────────────────
b_Wm3           = int(b)                     # 1500 W/m³
MW_MID_kWm3     = int(MW_MID / 1000)         # 20 kW/m³
MU_mPas         = MU * 1000                  # 3.5 mPa·s
Q_sensitivity_mL = Q0_peak * 1e6             # mL/s (peak pulsatile flow, = 4×resting)
ell0_mm         = int(ell0_coronary * 1000)  # 15 mm

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
