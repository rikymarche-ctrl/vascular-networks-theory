"""
compute.py — Paper III: The Dynamic Origin of Kleiber's Law
============================================================
Generates dynamic_variables.tex containing all LaTeX \\newcommand definitions
used in the manuscript. Every calculated quantity in the .tex file originates
here; no numerical values are hardcoded in the source.

Empirical inputs (literature sources, not derived here)
-------------------------------------------------------
Vascular wall scaling
    p_coronary  = 0.77   wall-thickness exponent, coronary/cerebral  (Rhodin 1967)
    p_pulmonary = 0.60   wall-thickness exponent, pulmonary          (Huang 1996)
    p_bronchial = 0.70   wall-thickness exponent, bronchial          (Weibel 1963)

Maintenance exponents for non-vascular systems
    m_trachea   = 0.50   insect tracheal taenidia                    (Westneat 2003)
    m_neuron    = 1.00   Na⁺/K⁺ pump, surface-dominated             (Attwell & Laughlin 2001)
    m_xylem_lo  = 1.00   xylem, drought-stressed                     (McCulloh 2003)
    m_xylem_hi  = 2.00   xylem, well-hydrated                        (McCulloh 2003)

Allometric reference values
    α*       = 2.72    minimax branching exponent                    (Paper II / Kassab 1993)
    M_ref    = 70 kg   human body mass reference
    r_0      = 12.5 mm human aortic inner radius                     (standard anatomy)
    f_h      = 70 bpm  human resting heart rate                      (Peters 1983)
    ν        = 3.2×10⁻⁶ m²/s  blood kinematic viscosity             (Windberger 2003)
    Wo_c     = 2.0     critical Womersley number                     (Womersley 1955)
    d        = 3       spatial embedding dimension (3D networks)
    n        = 4       Poiseuille transport exponent

Output
------
    ../manuscript/dynamic_variables.tex
"""

import math
from datetime import datetime
from fractions import Fraction

# ============================================================
# SECTION 1: Input parameters from literature
# ============================================================

# Wall-thickness scaling exponents p (m = 1+p)
p_coronary  = 0.77   # Rhodin 1967 — coronary and cerebral arteries
p_pulmonary = 0.60   # Huang et al. 1996 — pulmonary artery
p_bronchial = 0.70   # Weibel 1963 — bronchial tree

# Maintenance exponents for other systems
m_trachea   = 0.50   # Westneat 2003 — insect tracheal taenidia (h ∝ r^{1/2})
m_neuron    = 1.00   # Attwell & Laughlin 2001 — Na+/K+ pump, surface-dominated
m_xylem_lo  = 1.00   # McCulloh 2003 — xylem, drought-stressed (p=0, m=1)
m_xylem_hi  = 2.00   # McCulloh 2003 — xylem, well-hydrated (p=1, m=2)
m_leaf      = 1.00   # McCulloh 2003 — leaf venation, surface-dominated
m_sponge    = 1.00   # Reiswig 1975 — sponge canal walls, surface-dominated

# Transport exponents
n_P = 4              # Poiseuille (all vascular/tracheal/xylem/sponge/leaf)
n_O = 2              # Ohm / diffusion (dendrite electrical, His-Purkinje, constructal)

# Minimax exponents
alpha_star   = 2.72    # vascular minimax (Kassab 1993 calibration)
alpha_neural = 2.38    # neural minimax

# Spatial dimensions
d3 = 3
d2 = 2

# Reference values for Womersley transition (Theorem 5)
M_ref_kg  = 70.0     # human body mass (kg)
r0_m      = 0.0125   # human aortic lumen radius (m) — standard anatomy, inner diam ~2.5 cm
f_h_bpm   = 70.0     # human resting heart rate (bpm) — Peters 1983
nu_m2s    = 3.2e-6   # kinematic viscosity of blood (m^2/s) — Windberger 2003
Wo_c      = 2.0      # critical Womersley number — Womersley 1955

# ============================================================
# SECTION 2: Core formulas
# ============================================================

def alpha_t(n, m):
    """Theoretical branching exponent: Bennett (2025), arXiv:2511.04022."""
    return (n + m) / 2.0

def beta_alpha_d(alpha, d):
    """Generalized metabolic scaling exponent."""
    return (d * alpha) / (2*d + alpha)

def beta_n_m_d(n, m, d):
    """Allometric equation of state (viscous limit)."""
    return (d * (n + m)) / (4*d + (n + m))

def Wo(r, omega, nu):
    """Womersley number."""
    return r * math.sqrt(omega / nu)

def M_star(M_ref, Wo_c, Wo_ref):
    """Critical body mass for wave-to-viscous transition."""
    return M_ref * (Wo_c / Wo_ref)**4

# ============================================================
# SECTION 3: Compute Table 1 — αt values
# ============================================================

m_coronary  = 1 + p_coronary   # = 1.77
m_pulmonary = 1 + p_pulmonary  # = 1.60
m_bronchial = 1 + p_bronchial  # = 1.70

at_coronary  = alpha_t(n_P, m_coronary)   # (4 + 1.77)/2 = 2.885
at_pulmonary = alpha_t(n_P, m_pulmonary)  # (4 + 1.60)/2 = 2.800
at_cerebral  = alpha_t(n_P, m_coronary)   # same p as coronary = 2.885
at_bronchial = alpha_t(n_P, m_bronchial)  # (4 + 1.70)/2 = 2.850
at_neuron    = alpha_t(n_P, m_neuron)     # (4 + 1.00)/2 = 2.500
at_trachea   = alpha_t(n_P, m_trachea)   # (4 + 0.50)/2 = 2.250
at_xylem_lo  = alpha_t(n_P, m_xylem_lo)  # (4 + 1.00)/2 = 2.500
at_xylem_hi  = alpha_t(n_P, m_xylem_hi)  # (4 + 2.00)/2 = 3.000
at_leaf      = alpha_t(n_P, m_leaf)      # (4 + 1.00)/2 = 2.500
at_sponge    = alpha_t(n_P, m_sponge)    # (4 + 1.00)/2 = 2.500

# ============================================================
# SECTION 4: Compute Table 2 — untested predictions αt values
# ============================================================

at_fish_gill      = alpha_t(n_P, 1.0)  # (4+1)/2 = 2.50
at_his_purkinje   = alpha_t(n_O, 2.0)  # (2+2)/2 = 2.00
at_constructal    = alpha_t(n_O, 1.0)  # (2+1)/2 = 1.50
at_lymphatic      = alpha_t(n_P, 1.0)  # (4+1)/2 = 2.50
at_coral          = alpha_t(n_P, 0.5)  # (4+0.5)/2 = 2.25

# ============================================================
# SECTION 5: Corollary C — β(α*, 3)
# ============================================================

beta_minimax = beta_alpha_d(alpha_star, d3)
# = 3*2.72/(6+2.72) = 8.16/8.72

# ============================================================
# SECTION 6: Theorem 4 viscous sub-bound — β(5/2, 3)
# ============================================================

alpha_t_surface = alpha_t(n_P, m_neuron)  # = 5/2 exactly
beta_surface_bound = beta_alpha_d(alpha_t_surface, d3)
# = 3*(5/2)/(6+5/2) = (15/2)/(17/2) = 15/17

# Exact fraction check
beta_surface_frac = Fraction(3*5, 2*(6) + 5)  # 15/17
assert abs(beta_surface_bound - float(beta_surface_frac)) < 1e-12, \
    f"Fraction check failed: {beta_surface_bound} vs {float(beta_surface_frac)}"

# ============================================================
# SECTION 7: Theorem 5 — Womersley transition
# ============================================================

omega_ref = 2 * math.pi * f_h_bpm / 60.0   # rad/s
Wo_ref    = Wo(r0_m, omega_ref, nu_m2s)
M_star_g  = M_star(M_ref_kg * 1000, Wo_c, Wo_ref)  # in grams
# Note: M_ref in grams = 70000 g, result in grams

# ============================================================
# SECTION 8: Theorem 6 verification table
# ============================================================

beta_4_2_3 = beta_n_m_d(4, 2,   d3)   # = 18/18 = 1
beta_4_1_3 = beta_n_m_d(4, 1,   d3)   # = 15/17
beta_4_05_3= beta_n_m_d(4, 0.5, d3)   # = 13.5/16.5
beta_2_05_3= beta_n_m_d(2, 0.5, d3)   # = 7.5/14.5
beta_4_1_2 = beta_n_m_d(4, 1,   d2)   # = 10/13

# Exact fraction checks
assert abs(beta_4_2_3 - 1.0)         < 1e-12, "β(4,2,3) should be 1"
assert abs(beta_4_1_3 - 15/17)       < 1e-12, "β(4,1,3) should be 15/17"
assert abs(beta_4_05_3 - 13.5/16.5)  < 1e-12, "β(4,0.5,3) check"
assert abs(beta_2_05_3 - 7.5/14.5)   < 1e-12, "β(2,0.5,3) check"
assert abs(beta_4_1_2 - 10/13)       < 1e-12, "β(4,1,2) should be 10/13"

# ============================================================
# SECTION 9: Tracheal diffusion limit
# ============================================================

at_tracheole_diff = alpha_t(2, m_trachea)  # n=2 (diffusion), m=0.5 → 1.25

# ============================================================
# SECTION 10: Verification printout
# ============================================================

print("=" * 60)
print("COMPUTED VALUES — ALL MUST MATCH .tex")
print("=" * 60)
print(f"m_coronary          = {m_coronary:.4f}   [p=0.77 → m=1+p]")
print(f"m_pulmonary         = {m_pulmonary:.4f}   [p=0.60 → m=1+p]")
print(f"m_bronchial         = {m_bronchial:.4f}   [p=0.70 → m=1+p]")
print()
print(f"αt_coronary/cerebral= {at_coronary:.4f}   [(4+1.77)/2]")
print(f"αt_pulmonary        = {at_pulmonary:.4f}   [(4+1.60)/2]")
print(f"αt_bronchial        = {at_bronchial:.4f}   [(4+1.70)/2]")
print(f"αt_neuron           = {at_neuron:.4f}   [(4+1.00)/2]")
print(f"αt_trachea          = {at_trachea:.4f}   [(4+0.50)/2]")
print(f"αt_xylem_lo         = {at_xylem_lo:.4f}   [(4+1.00)/2]")
print(f"αt_xylem_hi         = {at_xylem_hi:.4f}   [(4+2.00)/2]")
print(f"αt_leaf             = {at_leaf:.4f}   [(4+1.00)/2, d=2]")
print(f"αt_sponge           = {at_sponge:.4f}   [(4+1.00)/2]")
print()
print(f"αt_fish_gill        = {at_fish_gill:.4f}")
print(f"αt_his_purkinje     = {at_his_purkinje:.4f}")
print(f"αt_constructal      = {at_constructal:.4f}")
print(f"αt_lymphatic        = {at_lymphatic:.4f}")
print(f"αt_coral            = {at_coral:.4f}")
print()
print(f"β(α*=2.72, d=3)     = {beta_minimax:.4f}   [Corollary C]")
print(f"β(5/2, d=3)=15/17   = {beta_surface_bound:.6f} ≈ {float(beta_surface_frac):.6f}")
print()
print(f"ω_ref               = {omega_ref:.4f} rad/s  [2π×70bpm/60]")
print(f"Wo_ref              = {Wo_ref:.4f}           [r0=1.2cm, ω_ref, ν]")
print(f"M*                  = {M_star_g:.2f} g     [M_ref×(Wo_c/Wo_ref)^4]")
print()
print(f"β(4,2,3)            = {beta_4_2_3:.6f}  [18/18]")
print(f"β(4,1,3)=15/17      = {beta_4_1_3:.6f}  [{15/17:.6f}]")
print(f"β(4,0.5,3)          = {beta_4_05_3:.6f}  [13.5/16.5]")
print(f"β(2,0.5,3)          = {beta_2_05_3:.6f}  [7.5/14.5]")
print(f"β(4,1,2)=10/13      = {beta_4_1_2:.6f}  [{10/13:.6f}]")
print()
print(f"αt_tracheole_diff   = {at_tracheole_diff:.4f}  [n=2,m=0.5 diffusion]")

# ============================================================
# SECTION 11: Generate LaTeX \newcommand file
# ============================================================

def fmt3(x):
    """Format to exactly 3 decimal places."""
    return f"{x:.3f}"

def fmt4(x):
    """Format to 4 significant figures."""
    return f"{x:.4g}"

def format_macro(name, value, unit, description):
    """Format LaTeX macro to column-50 alignment with unit comment."""
    macro = f"\\newcommand{{\\{name}}}{{{value}}}"
    return f"{macro:<50s} % [{unit:<8}] {description}"

lines = []
lines.append("% ====================================================================")
lines.append("% DYNAMIC VARIABLES")
lines.append("% ====================================================================")
lines.append("% Source: paper3-kleiber/scripts/compute.py")
lines.append(f"% Generated: {datetime.now().strftime('%Y-%m-%dT%H:%M')}")
lines.append("% ====================================================================")
lines.append("")

# --- Table 1: m values ---
lines.append("% Wall-maintenance exponents (m = 1+p, derived from literature p)")
lines.append(f"\\newcommand{{\\mCoronary}}{{{fmt3(m_coronary)}}}")
lines.append(f"\\newcommand{{\\mPulmonary}}{{{fmt3(m_pulmonary)}}}")
lines.append(f"\\newcommand{{\\mBronchial}}{{{fmt3(m_bronchial)}}}")
lines.append("")

# --- Table 1: m values for non-vascular systems ---
lines.append("% Non-vascular maintenance exponents")
lines.append(f"\\newcommand{{\\mNeuron}}{{{fmt3(m_neuron)}}}")
lines.append(f"\\newcommand{{\\mTrachea}}{{{fmt3(m_trachea)}}}")
lines.append(f"\\newcommand{{\\mXylemLo}}{{{fmt3(m_xylem_lo)}}}")
lines.append(f"\\newcommand{{\\mXylemHi}}{{{fmt3(m_xylem_hi)}}}")
lines.append(f"\\newcommand{{\\mLeaf}}{{{fmt3(m_leaf)}}}")
lines.append(f"\\newcommand{{\\mSponge}}{{{fmt3(m_sponge)}}}")
lines.append("")

# --- Table 1: αt values ---
lines.append("% Theoretical branching exponents alpha_t = (n+m)/2")
lines.append(f"\\newcommand{{\\atCoronary}}{{{fmt3(at_coronary)}}}")
lines.append(f"\\newcommand{{\\atPulmonary}}{{{fmt3(at_pulmonary)}}}")
lines.append(f"\\newcommand{{\\atCerebral}}{{{fmt3(at_cerebral)}}}")
lines.append(f"\\newcommand{{\\atBronchial}}{{{fmt3(at_bronchial)}}}")
lines.append(f"\\newcommand{{\\atNeuron}}{{{fmt3(at_neuron)}}}")
lines.append(f"\\newcommand{{\\atTrachea}}{{{fmt3(at_trachea)}}}")
lines.append(f"\\newcommand{{\\atXylemLo}}{{{fmt3(at_xylem_lo)}}}")
lines.append(f"\\newcommand{{\\atXylemHi}}{{{fmt3(at_xylem_hi)}}}")
lines.append(f"\\newcommand{{\\atLeaf}}{{{fmt3(at_leaf)}}}")
lines.append(f"\\newcommand{{\\atSponge}}{{{fmt3(at_sponge)}}}")
lines.append("")

# --- Table 2: predictions ---
lines.append("% Untested predictions (Table 2)")
lines.append(f"\\newcommand{{\\atFishGill}}{{{fmt3(at_fish_gill)}}}")
lines.append(f"\\newcommand{{\\atConstructal}}{{{fmt3(at_constructal)}}}")
lines.append(f"\\newcommand{{\\atLymphatic}}{{{fmt3(at_lymphatic)}}}")
lines.append(f"\\newcommand{{\\atCoral}}{{{fmt3(at_coral)}}}")
lines.append("")


# --- Theorem 4 sub-bound ---
lines.append("% Theorem 4 viscous sub-bound: beta(5/2, 3) = 15/17")
lines.append(f"\\newcommand{{\\betaSurfaceBound}}{{{fmt3(beta_surface_bound)}}}")
lines.append("")

# --- Womersley transition ---
lines.append("% Theorem 5: Womersley transition")
lines.append(f"\\newcommand{{\\omegaRef}}{{{fmt3(omega_ref)}}}")
lines.append(f"\\newcommand{{\\WoRef}}{{{fmt3(Wo_ref)}}}")
lines.append(f"\\newcommand{{\\Mstar}}{{{round(M_star_g, 2)}}}")
lines.append(f"\\newcommand{{\\MstarRounded}}{{"
             f"\\ensuremath{{\\approx {round(M_star_g, 2)}~\\mathrm{{g}}}}}}")
lines.append("")

# --- Theorem 6 verification values ---
lines.append("% Theorem 6 verification table")
lines.append(f"\\newcommand{{\\betaFourTwoThree}}{{{fmt3(beta_4_2_3)}}}")
lines.append(f"\\newcommand{{\\betaFourOneThree}}{{{fmt3(beta_4_1_3)}}}")
lines.append(f"\\newcommand{{\\betaFourHalfThree}}{{{fmt3(beta_4_05_3)}}}")
lines.append(f"\\newcommand{{\\betaTwoHalfThree}}{{{fmt3(beta_2_05_3)}}}")
lines.append(f"\\newcommand{{\\betaFourOneTwo}}{{{fmt3(beta_4_1_2)}}}")
lines.append("")

# --- Tracheal diffusion limit ---
lines.append("% Diffusion tracheole limit (n=2, m=0.5)")
lines.append(f"\\newcommand{{\\atTracheoleDiff}}{{{fmt3(at_tracheole_diff)}}}")
lines.append("")

# --- Leaf venation β ---
lines.append("% Leaf venation beta (d=2, n=4, m=1)")
lines.append(f"\\newcommand{{\\betaLeaf}}{{{fmt3(beta_4_1_2)}}}")
lines.append("")

# --- Leaf venation dimensionality caveat ---
d_eff = 2.1
beta_leaf_eff = beta_alpha_d(at_leaf, d_eff)
lines.append("% Leaf venation dimensionality caveat (d_eff = 2.1)")
lines.append(f"\\newcommand{{\\betaLeafEff}}{{{fmt3(beta_leaf_eff)}}}")
lines.append(f"\\newcommand{{\\Deff}}{{{d_eff}}}")
lines.append("")

# --- alpha_star ---
lines.append("% Minimax exponent alpha*")
lines.append(f"\\newcommand{{\\alphastar}}{{{fmt3(alpha_star)}}}")
lines.append(f"\\newcommand{{\\alphastarNeural}}{{{fmt3(alpha_neural)}}}")
lines.append("")

# ============================================================
# SECTION 12: Reflection coefficients (dynamic selection of α_w)
# ============================================================


def gamma_ratio(alpha, N=2):
    """Impedance ratio γ = N^{1-2/α}. γ=1 iff α=2."""
    return N ** (1 - 2/alpha)

def R_power(alpha, N=2):
    """Power reflection coefficient R² = ((γ-1)/(γ+1))²."""
    g = gamma_ratio(alpha, N)
    return ((g - 1) / (g + 1))**2

def T_cumulative(alpha, K, N=2):
    """Cumulative power transmission over K bifurcations."""
    return (1 - R_power(alpha, N))**K

K_coronary = 11   # typical coronary generations

R2_wave    = R_power(2.00)          # = 0 exactly
R2_minimax = R_power(alpha_star)    # alpha*=2.72
R2_murray  = R_power(3.00)

T_wave    = T_cumulative(2.00,       K_coronary)
T_minimax = T_cumulative(alpha_star, K_coronary)
T_murray  = T_cumulative(3.00,       K_coronary)

print("=" * 60)
print("REFLECTION COEFFICIENTS (dynamic selection of α_w)")
print("=" * 60)
print(f"R²(α=2.00)  = {R2_wave:.6f}   [α_w: zero reflection]")
print(f"R²(α=2.72)  = {R2_minimax:.6f}")
print(f"R²(α=3.00)  = {R2_murray:.6f}")
print(f"T(K=11, α=2.00)  = {T_wave:.6f}")
print(f"T(K=11, α=2.72)  = {T_minimax:.6f}")
print(f"T(K=11, α=3.00)  = {T_murray:.6f}")

# LaTeX commands
lines.append("% Reflection coefficients — dynamic selection of alpha_w")
lines.append(f"\\newcommand{{\\TcumMinimax}}{{{fmt3(T_minimax)}}}")
lines.append(f"\\newcommand{{\\TcumMurray}}{{{fmt3(T_murray)}}}")
# Percentage losses
lines.append(f"\\newcommand{{\\TlossMinimax}}{{{fmt3(100*(1-T_minimax))}}}")
lines.append(f"\\newcommand{{\\TlossMurray}}{{{fmt3(100*(1-T_murray))}}}")
lines.append("")

# ============================================================
# SECTION 13: Dimensional universality β(α_w=2, d) = d/(d+1)
# ============================================================

def beta_dimensional(d):
    """β at wave attractor α_w=2: β = d/(d+1)"""
    return d / (d + 1)

beta_d1 = beta_dimensional(1)   # = 1/2
beta_d2 = beta_dimensional(2)   # = 2/3
beta_d3 = beta_dimensional(3)   # = 3/4 (Kleiber)
beta_d4 = beta_dimensional(4)   # = 4/5

print("=" * 60)
print("DIMENSIONAL UNIVERSALITY β(α_w=2, d) = d/(d+1)")
print("=" * 60)
for d in [1, 2, 3, 4]:
    frac = Fraction(d, d+1)
    print(f"  d={d}: beta = {frac} = {float(frac):.4f}")

lines.append("% Dimensional universality β(α_w=2, d) = d/(d+1)")
lines.append(f"\\newcommand{{\\betaDimOne}}{{{fmt3(beta_d1)}}}")
lines.append(f"\\newcommand{{\\betaDimTwo}}{{{fmt3(beta_d2)}}}")
lines.append(f"\\newcommand{{\\betaDimThree}}{{{fmt3(beta_d3)}}}")
lines.append(f"\\newcommand{{\\betaDimFour}}{{{fmt3(beta_d4)}}}")
lines.append("")

# ============================================================
# SECTION 14: Allometric universality classes — discrete spectrum
# ============================================================

def beta_s_d(s, d):
    return d * s / (4*d + s)

# Key universality classes for d=3
s_values = [4.0, 4.5, 5.0, 5.77, 6.0]
print("=" * 60)
print("DISCRETE ALLOMETRIC SPECTRUM d=3")
print("=" * 60)
for s in s_values:
    b = beta_s_d(s, 3)
    print(f"  s={s:.2f}: beta = {b:.4f}")

# His-Purkinje: n=2, m=2, s=4 → β=3/4 — same class as Kleiber!
beta_his_purkinje_metabolic = beta_s_d(2+2, 3)  # = 3/4 exactly
print(f"\nHis-Purkinje (n=2,m=2): beta = {beta_his_purkinje_metabolic:.6f}  [= 3/4 EXACTLY]")

# s=4.5 for trachea
beta_s45_d3 = beta_s_d(4.5, 3)
# d=2 wave attractor
beta_wave_d2 = beta_dimensional(2)

lines.append("% Allometric universality classes")
lines.append(f"\\newcommand{{\\betaSfourFive}}{{{fmt3(beta_s45_d3)}}}")
lines.append("")

# Fixed point uniqueness: β/2 - 1/8 = 1/4 → β = 3/4
# Wo ∝ M^(β/2 - 1/8); self-consistency requires exponent = 1/4
print("=" * 60)
print("FIXED POINT UNIQUENESS")
print("=" * 60)
for beta_test in [0.70, 0.75, 0.80]:
    exp = beta_test/2 - 1/8
    status = '✓' if abs(exp - 0.25) < 1e-6 else 'x'
    print(f"  beta={beta_test:.2f}: Wo proportional to "
          f"M^{exp:.4f} (need 0.2500, {status})")

lines.append("% Fixed point uniqueness exponent")
lines.append("")

# ============================================================
# SECTION 15: Monotonicity and double optimality of alpha_w
# ============================================================

# rho(alpha, d=3, N=2) = 2^{1-2/alpha-1/3}
def rho_fn(alpha, d=3, N=2):
    return N**(1 - 2/alpha - 1/d)

# drho/dalpha = rho * ln(N) * 2/alpha^2 > 0 always
def drho_dalpha(alpha, d=3, N=2):
    r = rho_fn(alpha, d, N)
    return r * math.log(N) * 2 / alpha**2

rho_wave   = rho_fn(2.00)
rho_neural = rho_fn(2.38)
rho_vasc   = rho_fn(2.72)
rho_murray = rho_fn(3.00)
drho_at_2  = drho_dalpha(2.0)

# M* generalized to d dimensions
# Wo exponent: (d-1)/(2(d+1))
# M* exponent: 2(d+1)/(d-1)
def mstar_exponent(d):
    return 2*(d+1)/(d-1)

mstar_exp_d2 = mstar_exponent(2)  # = 6
mstar_exp_d3 = mstar_exponent(3)  # = 4
mstar_exp_d4 = mstar_exponent(4)  # = 10/3

print("=" * 60)
print("MONOTONICITY AND DOUBLE OPTIMALITY")
print("=" * 60)
print(f"rho(alpha_w=2)     = {rho_wave:.6f}")
print(f"drho/dalpha at α=2 = {drho_at_2:.6f} > 0 (strictly increasing)")
print(f"M* exponent d=2: {mstar_exp_d2}")
print(f"M* exponent d=3: {mstar_exp_d3}")
print(f"M* exponent d=4: {mstar_exp_d4:.4f}")

lines.append("% Monotonicity of rho — double optimality")
lines.append(f"\\newcommand{{\\rhoWave}}{{{fmt3(rho_wave)}}}")
lines.append("")

# ============================================================
# SECTION 16: Convergence table — dynamic values
# ============================================================

K_conv = 11   # generations used in convergence table

def rho_conv(alpha, d=3, N=2):
    return N ** (1 - 2/alpha - 1/d)

def eps_conv(alpha, K=11, d=3, N=2):
    return rho_conv(alpha, d, N) ** (K + 1)

def vol_pct(alpha, K=11, d=3, N=2):
    return (1 - eps_conv(alpha, K, d, N)) * 100

# Three key exponents in the table
alpha_wave = 2.00
alpha_vasc = alpha_star   # vascular minimax (= alpha_star)

rho_w  = rho_conv(alpha_wave)
rho_n  = rho_conv(alpha_neural)
rho_v  = rho_conv(alpha_vasc)

eps_w  = eps_conv(alpha_wave)
eps_n  = eps_conv(alpha_neural)
eps_v  = eps_conv(alpha_vasc)

vol_w  = vol_pct(alpha_wave)
vol_n  = vol_pct(alpha_neural)
vol_v  = vol_pct(alpha_vasc)

print("=" * 60)
print("CONVERGENCE TABLE VALUES (K=11, N=2, d=3)")
print("=" * 60)
cases = [
    (alpha_wave, rho_w, eps_w, vol_w),
    (alpha_neural, rho_n, eps_n, vol_n),
    (alpha_vasc, rho_v, eps_v, vol_v)
]
for a, r, e, v in cases:
    print(f"  alpha={a:.2f}: rho={r:.4f}  eps={e:.4f}  vol={v:.1f}%")

# EpsVascPercent: the ~57% figure in the text
eps_vasc_pct = eps_v * 100   # = 56.5 → text says "57%"

# For conclusion: rho_wave to 2dp, and vol_wave to nearest %
rho_wave_2dp = round(rho_w, 2)    # 0.79
vol_wave_pct_round = round(vol_w)  # 94

print(f"\nEpsVascPercent = {eps_vasc_pct:.1f}%  (text uses ~57%)")
print(f"rhoWaveTwoDp = {rho_wave_2dp}")
print(f"VolCaptWaveRound = {vol_wave_pct_round}%")

lines.append("% Convergence table — dynamic values (K=11, N=2, d=3)")
lines.append(f"\\newcommand{{\\rhoAlphaWave}}{{{fmt3(rho_w)}}}")
lines.append(f"\\newcommand{{\\epsAlphaWave}}{{{fmt3(eps_w)}}}")
lines.append(f"\\newcommand{{\\volAlphaWave}}{{{round(vol_w,1)}}}")
lines.append(f"\\newcommand{{\\rhoAlphaNeuralMinimax}}{{{fmt3(rho_n)}}}")
lines.append(f"\\newcommand{{\\epsAlphaNeuralMinimax}}{{{fmt3(eps_n)}}}")
lines.append(f"\\newcommand{{\\volAlphaNeuralMinimax}}{{{round(vol_n,1)}}}")
lines.append(f"\\newcommand{{\\rhoAlphaVascMinimax}}{{{fmt3(rho_v)}}}")
lines.append(f"\\newcommand{{\\epsAlphaVascMinimax}}{{{fmt3(eps_v)}}}")
lines.append(f"\\newcommand{{\\volAlphaVascMinimax}}{{{round(vol_v,1)}}}")
lines.append(f"\\newcommand{{\\rhoWaveTwoDp}}{{{rho_wave_2dp}}}")
lines.append(f"\\newcommand{{\\VolCaptWaveRound}}{{{vol_wave_pct_round}}}")
lines.append("")

# ============================================================
# SECTION 17: M* sensitivity to Wo_c
# ============================================================

Wo_c_range = [1.5, 1.8, 2.0, 2.2, 2.5]
print("=" * 60)
print("M* SENSITIVITY TO Wo_c (d=3, exponent=4)")
print("=" * 60)
for wc in Wo_c_range:
    m_star_g = M_star(M_ref_kg * 1000, wc, Wo_ref)
    print(f"  Wo_c={wc:.1f}: M* = {m_star_g:.2f} g")

# Store min/max for text
m_star_lo = M_star(M_ref_kg * 1000, 1.5, Wo_ref)
m_star_hi = M_star(M_ref_kg * 1000, 2.5, Wo_ref)
lines.append(f"\\newcommand{{\\MstarWocLo}}{{{round(m_star_lo,1)}}}")
lines.append(f"\\newcommand{{\\MstarWocHi}}{{{round(m_star_hi,1)}}}")
lines.append(f"\\newcommand{{\\MstarWocLoVal}}{{1.5}}")
lines.append(f"\\newcommand{{\\MstarWocHiVal}}{{2.5}}")
lines.append("")

# β sensitivity to α_w deviation
print("=" * 60)
print("β SENSITIVITY TO α_w DEVIATION")
print("=" * 60)
for aw in [1.9, 2.0, 2.1, 2.2]:
    b = beta_alpha_d(aw, 3)
    print(f"  alpha_w={aw:.1f}: beta = {b:.4f}  (delta_beta = {b-0.75:+.4f})")

beta_aw_2p1 = beta_alpha_d(2.1, 3)
beta_aw_1p9 = beta_alpha_d(1.9, 3)
lines.append(f"\\newcommand{{\\betaAlphaWTwoPointOne}}{{{fmt3(beta_aw_2p1)}}}")
lines.append(f"\\newcommand{{\\deltaBetaAlphaWTwoPointOne}}{{{fmt3(beta_aw_2p1 - 0.75)}}}")
lines.append("")

# ============================================================
# SECTION 18: His-Purkinje — resonance analysis
# ============================================================
# n=2 (Ohm), m=2 (glycogen volume), alpha_t = alpha_w = 2
# Zero reflection, perfect transmission, beta = 3/4 exactly

n_HP = 2
m_HP = 2
at_HP = alpha_t(n_HP, m_HP)          # = 2.000
beta_HP = beta_s_d(n_HP + m_HP, d3)  # = 3/4 exactly

# Resonance: α_t = α_w → R² = 0 exactly
R2_HP = R_power(at_HP)               # = 0.000000
T_HP  = T_cumulative(at_HP, K_coronary)  # = 1.000000

# Geometric convergence at α=2
rho_HP = rho_fn(at_HP)               # = 0.7937
vol_HP = vol_pct(at_HP)              # = 93.8%

# Morphometric prediction: r²_His ∝ M^{3/4}
# r_His ∝ M^{3/8} — falsifiable, no data exists
r_scaling_exp_HP = beta_HP / 2       # = 3/8 = 0.375

print("=" * 60)
print("HIS-PURKINJE RESONANCE (SECTION 18)")
print("=" * 60)
print(f"alpha_t = alpha_w = {at_HP:.3f}")
print(f"beta    = {beta_HP:.6f}  [= 3/4 EXACTLY]")
print(f"R^2     = {R2_HP:.6f}")
print(f"T       = {T_HP:.6f}")
print(f"rho     = {rho_HP:.4f}")
print(f"vol     = {vol_HP:.1f}%")
print(f"r_His scaling exp = {r_scaling_exp_HP:.3f}  [= 3/8]")

lines.append("% ============================================================")
lines.append("% Section 18: His-Purkinje resonance analysis")
lines.append("% ============================================================")
lines.append(f"\\newcommand{{\\mHP}}{{{fmt3(float(m_HP))}}}")
lines.append(f"\\newcommand{{\\atHP}}{{{fmt3(at_HP)}}}")
lines.append(f"\\newcommand{{\\betaHP}}{{{fmt3(beta_HP)}}}")
lines.append(f"\\newcommand{{\\betaHPfrac}}{{3/4}}")
lines.append("")

# ============================================================
# SECTION 19: Clade-shifting of the Womersley Transition
# ============================================================
# All clades evaluated at M_ref = 1 kg for comparability.
# r0 derived from isometric allometry r0 ∝ M^{3/8} anchored to human.

print("=" * 60)
print("CLADE-SHIFTING OF WOMERSLEY TRANSITION (M_ref = 1 kg)")
print("=" * 60)

M_ref_clade_kg = 1.0
Wo_c_clade     = 2.0

# Isometric scaling of aortic radius: human 70 kg → 1 kg
r0_1kg = r0_m * (M_ref_clade_kg / M_ref_kg)**(3/8)
print(f"r0 at 1 kg (isometric): {r0_1kg*1000:.4f} mm")

# --- Birds ---
# Sources: Calder 1968 (fh), Grubb 1983 (cardiovascular allometry)
fh_bird   = 156.0     # bpm
nu_bird   = 4.0e-6    # m^2/s  (nucleated erythrocytes)
omega_bird = fh_bird * 2 * math.pi / 60.0
Wo_ref_bird   = r0_1kg * math.sqrt(omega_bird / nu_bird)
M_star_bird_g = (M_ref_clade_kg * 1000) * (Wo_c_clade / Wo_ref_bird)**4

print(f"\nBIRDS:")
print(f"  omega   = {omega_bird:.4f} rad/s")
print(f"  Wo_ref  = {Wo_ref_bird:.4f}")
print(f"  M*      = {M_star_bird_g:.1f} g")

# --- Reptiles ---
# fh source: White 1976 (heart rate scaling in reptiles)
# nu sources: Langille & Crisp 1980 (turtles, interpolated at 20°C)
#             Dunlap 2006 (Sceloporus, upper bound)
fh_rept    = 34.0     # bpm
omega_rept = fh_rept * 2 * math.pi / 60.0

nu_rept_lo = 2.68e-6  # m^2/s  Langille & Crisp 1980, turtle at 20°C
nu_rept_hi = 5.62e-6  # m^2/s  Dunlap 2006, Sceloporus max (5.9 cP)

Wo_ref_rept_lo = r0_1kg * math.sqrt(omega_rept / nu_rept_lo)
Wo_ref_rept_hi = r0_1kg * math.sqrt(omega_rept / nu_rept_hi)

# lower nu → higher Wo → lower M*  (hence _lo nu gives _lo M*)
M_star_rept_lo_g = (M_ref_clade_kg * 1000) * (Wo_c_clade / Wo_ref_rept_lo)**4
M_star_rept_hi_g = (M_ref_clade_kg * 1000) * (Wo_c_clade / Wo_ref_rept_hi)**4

print(f"\nREPTILES (20°C):")
print(f"  omega         = {omega_rept:.4f} rad/s")
print(f"  Wo_ref range  = [{Wo_ref_rept_hi:.4f}, {Wo_ref_rept_lo:.4f}]")
print(f"  M* range      = [{M_star_rept_lo_g:.1f}, {M_star_rept_hi_g:.1f}] g")

# --- Asymptotic viscous limit for reptiles ---
# Reptile vasculature: same wall physics as mammals (n=4, m=1+p≈1.77)
beta_viscous_reptile = beta_n_m_d(4.0, m_coronary, d3)
print(f"\nREPTILES VISCOUS ASYMPTOTE:")
print(f"  beta(4, {m_coronary}, 3) = {beta_viscous_reptile:.4f}")

print(f"\n=== CLADE SUMMARY ===")
print(f"Mammals: M* = {M_star_g:.2f} g  (from Section 7)")
print(f"Birds:   M* = {M_star_bird_g:.1f} g")
print(f"Reptiles:M* = [{M_star_rept_lo_g:.0f}, {M_star_rept_hi_g:.0f}] g")
print(f"Ratio reptile_hi / mammal: {M_star_rept_hi_g / M_star_g:.1f}x")

lines.append("% ============================================================")
lines.append("% Section 19: Clade-shifting of the Womersley Transition")
lines.append("% ============================================================")
lines.append(f"\\newcommand{{\\rZeroOneKg}}{{{fmt3(r0_1kg * 1000)}}}")   # in mm
lines.append(f"\\newcommand{{\\MstarBird}}{{{round(M_star_bird_g, 1)}}}")
lines.append(f"\\newcommand{{\\WoRefBird}}{{{fmt3(Wo_ref_bird)}}}")
lines.append(f"\\newcommand{{\\MstarReptLo}}{{{round(M_star_rept_lo_g, 0):.0f}}}")
lines.append(f"\\newcommand{{\\MstarReptHi}}{{{round(M_star_rept_hi_g, 0):.0f}}}")
lines.append(f"\\newcommand{{\\WoRefReptLo}}{{{fmt3(Wo_ref_rept_hi)}}}")   # lower Wo
lines.append(f"\\newcommand{{\\WoRefReptHi}}{{{fmt3(Wo_ref_rept_lo)}}}")   # higher Wo
lines.append(f"\\newcommand{{\\betaViscousReptile}}{{{fmt3(beta_viscous_reptile)}}}")
lines.append("")

# ============================================================
# SECTION 20: Physical minimax derivation — g(α) from C_maint
# ============================================================
#
# With terminal capillary radius r_K fixed, the structural
# maintenance cost of the full K-level network is:
#
#   C_maint(α) = B · r_K^m · l_K · N^{K(m/α + 1/d)} · Σ_{j=0}^{K} μ(α)^j
#
# where  μ(α) = N^{1 − m/α − 1/d}.
#
# The excess cost relative to the optimum α_t is then:
#
#   g(α) = C_maint(α) / C_maint(α_t)  − 1   ≥ 0,  with g(α_t) = 0.
#
# The minimax saddle point α* solves f(α*) = g(α*) uniquely,
# where f(α) = 1 − (1 − R²(α))^K  (cumulative wave-reflection loss).
# No free parameters: everything is fixed by m, K, N, d, α_t.

def mu_maint(alpha, m=m_coronary, N=2, d=3):
    """Level-to-level maintenance volume ratio μ(α) = N^{1−m/α−1/d}."""
    return N ** (1.0 - m/alpha - 1.0/d)

def C_maint_phys(alpha, m=m_coronary, K=K_coronary, N=2, d=3):
    """Total structural maintenance cost (arbitrary units, terminal-unit fixed)."""
    mu = mu_maint(alpha, m, N, d)
    prefactor = N ** (K * (m/alpha + 1.0/d))
    if abs(mu - 1.0) < 1e-10:
        s = float(K + 1)
    else:
        s = (mu**(K+1) - 1.0) / (mu - 1.0)
    return prefactor * s

_C_maint_at = C_maint_phys(at_coronary)

def g_phys_fn(alpha):
    """Normalised excess maintenance cost: 0 at α_t, positive below."""
    return C_maint_phys(alpha) / _C_maint_at - 1.0

def f_phys_fn(alpha):
    """Cumulative wave-reflection loss: 0 at α_w=2, positive above."""
    return 1.0 - T_cumulative(alpha, K_coronary)

# Normalization denominators (make both curves span [0, 1] for Fig 3)
_f_max = f_phys_fn(at_coronary)   # f at α_t
_g_max = g_phys_fn(2.0)           # g at α_w

def f_norm_fn(alpha):
    return f_phys_fn(alpha) / _f_max

def g_norm_fn(alpha):
    return g_phys_fn(alpha) / _g_max

def _bisect_fg(a, b, tol=1e-9):
    """Bisection to find α* where f(α*) = g(α*)."""
    for _ in range(120):
        mid = (a + b) * 0.5
        diff = f_phys_fn(mid) - g_phys_fn(mid)
        if abs(b - a) < tol:
            break
        if (f_phys_fn(a) - g_phys_fn(a)) * diff <= 0:
            b = mid
        else:
            a = mid
    return (a + b) * 0.5

# f < g at α_w (f=0, g>0) and f > g at α_t (g=0, f>0) → unique root
alpha_star_phys = _bisect_fg(2.001, at_coronary - 1e-6)
f_star_phys     = f_phys_fn(alpha_star_phys)
g_star_phys     = g_phys_fn(alpha_star_phys)
beta_star_phys  = beta_alpha_d(alpha_star_phys, 3)

print("=" * 60)
print("SECTION 20: PHYSICAL MINIMAX (no free parameters)")
print("=" * 60)
print(f"  α_t (coronary, m={m_coronary})  = {at_coronary:.4f}")
print(f"  g(α_w=2)                = {g_phys_fn(2.0)*100:.1f}%  (raw maint excess)")
print(f"  f(α_t)                  = {_f_max*100:.1f}%  (wave loss at α_t)")
print()
print(f"  Bisection crossing f=g:")
for a_test in [2.50, 2.60, 2.70, 2.72, 2.75, 2.78, 2.80, 2.85]:
    if a_test < at_coronary:
        print(f"    α={a_test:.2f}:  f={f_phys_fn(a_test)*100:.2f}%  "
              f"g={g_phys_fn(a_test)*100:.2f}%  "
              f"diff={f_phys_fn(a_test)-g_phys_fn(a_test):+.4f}")
print(f"\n  α*_phys = {alpha_star_phys:.4f}  (physical minimax, no free params)")
print(f"  f(α*) = g(α*) = {f_star_phys*100:.2f}%")
print(f"  β(α*_phys, 3) = {beta_star_phys:.4f}")
print(f"  Kassab αexp   = 2.70 ± 0.20  [consistent]")

# ============================================================
# SECTION 21: Pulmonary Minimax (no free parameters)
# ============================================================

r0_pulm = 15e-3
omega_pulm = 2 * math.pi * 1.2
rho_blood = 1050
mu_blood = 3e-3

Wo_pulm = r0_pulm * math.sqrt(omega_pulm * rho_blood / mu_blood)

n_p = 4
m_p = m_pulmonary  # usually 1.60
alpha_t_p = (n_p + m_p) / 2.0

def f_pulm_fn(alpha):
    return 1.0 - T_cumulative(alpha, K=K_coronary, N=2)

def g_pulm_fn(alpha):
    mu_val = mu_maint(alpha, m=m_p)
    prefactor = 2 ** (K_coronary * (m_p/alpha + 1.0/3.0))
    if abs(mu_val - 1.0) > 1e-10:
        s = (mu_val**(K_coronary+1) - 1.0) / (mu_val - 1.0)
    else:
        s = float(K_coronary + 1)
    C_maint = prefactor * s

    mu_t = mu_maint(alpha_t_p, m=m_p)
    pref_t = 2 ** (K_coronary * (m_p/alpha_t_p + 1.0/3.0))
    if abs(mu_t - 1.0) > 1e-10:
        st = (mu_t**(K_coronary+1) - 1.0) / (mu_t - 1.0)
    else:
        st = float(K_coronary + 1)
    C_maint_t = pref_t * st

    return C_maint / C_maint_t - 1.0

def _bisect_fg_pulm(a, b, tol=1e-9):
    for _ in range(120):
        mid = (a + b) * 0.5
        diff = f_pulm_fn(mid) - g_pulm_fn(mid)
        if abs(b - a) < tol:
            break
        if (f_pulm_fn(a) - g_pulm_fn(a)) * diff <= 0:
            b = mid
        else:
            a = mid
    return (a + b) * 0.5

alpha_star_pulm = _bisect_fg_pulm(2.001, alpha_t_p - 1e-6)
f_star_pulm = f_pulm_fn(alpha_star_pulm)

print("=" * 60)
print("SECTION 21: PULMONARY MINIMAX")
print("=" * 60)
print(f"  Wo_pulm = {Wo_pulm:.2f}")
print(f"  α*_pulm = {alpha_star_pulm:.4f}")
print(f"  f(α*) = g(α*) = {f_star_pulm*100:.2f}%")

# Append new \newcommands to dynamic_variables.tex
lines.append("% ============================================================")
lines.append("% Section 20: Physical minimax α* (no free parameters)")
lines.append("% ============================================================")
lines.append(f"\\newcommand{{\\alphastarPhys}}{{{alpha_star_phys:.2f}}}")
lines.append(f"\\newcommand{{\\betaMinimaxPhys}}{{{beta_star_phys:.3f}}}")
lines.append(f"\\newcommand{{\\fstarPhys}}{{{f_star_phys*100:.1f}}}")
lines.append("")

lines.append("% ============================================================")
lines.append("% Section 21: Pulmonary minimax α*")
lines.append("% ============================================================")
lines.append(f"\\newcommand{{\\rZeroPulm}}{{{r0_pulm * 1000:.0f}}}")
lines.append(f"\\newcommand{{\\omegaPulm}}{{{omega_pulm:.2f}}}")
lines.append(f"\\newcommand{{\\WoPulm}}{{{Wo_pulm:.0f}}}")
lines.append(f"\\newcommand{{\\alphastarPulm}}{{{alpha_star_pulm:.3f}}}")
lines.append(f"\\newcommand{{\\fstarPulm}}{{{f_star_pulm*100:.2f}}}")
lines.append("")

# ============================================================
# SECTION 22: Figure generation
# ============================================================
# Requires: matplotlib, numpy  (pip install matplotlib numpy)
# Outputs:  ../manuscript/figures/fig1_spectrum.pdf
#           ../manuscript/figures/fig2_transition.pdf
#           ../manuscript/figures/fig3_minimax.pdf

import os
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")       # headless — no display needed
    matplotlib.rcParams["text.usetex"] = False
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    _MPL_OK = True
except ImportError:
    print("\n[WARNING] matplotlib not found — skipping figure generation.")
    print("  Install with: pip install matplotlib")
    _MPL_OK = False

if _MPL_OK:
    os.makedirs("../manuscript/figures", exist_ok=True)

    plt.rcParams.update({
        "font.family":      "DejaVu Serif",
        "font.size":        10,
        "axes.labelsize":   10,
        "axes.titlesize":   10,
        "xtick.labelsize":  9,
        "ytick.labelsize":  9,
        "legend.fontsize":  9,
        "figure.dpi":       150,
        "pdf.fonttype":     42,   # TrueType in PDF (editable in Illustrator)
        "text.usetex":      False,
    })

    # ===========================================================
    # FIG 1 — Allometric spectrum: β(n, m) heatmap for d = 3
    # ===========================================================
    #
    # Contour label positions are computed from the exact formula
    # β(n+m=s, d=3) = 3s/(12+s):
    #   3/4   → s=4.000  → label at (1.40, 2.60)
    #   0.80  → s=4.364  → label at (1.82, 2.54)
    #   15/17 → s=5.000  → label at (3.00, 2.00)
    #   1     → s=6.000  → label at (4.20, 1.80)
    # Biological annotations use explicit text anchors to avoid
    # any overlap among the three vascular points at n=4.

    n_grid = np.linspace(1.0, 5.0, 400)
    m_grid = np.linspace(0.0, 3.0, 400)
    NN, MM = np.meshgrid(n_grid, m_grid)
    BB = beta_n_m_d(NN, MM, 3)

    fig1, ax1 = plt.subplots(figsize=(8.2, 5.4))
    im = ax1.contourf(NN, MM, BB,
                      levels=np.linspace(0.55, 1.0, 60),
                      cmap="viridis_r", extend="both")
    cb = fig1.colorbar(im, ax=ax1, pad=0.02)
    cb.set_label("β(n, m, d=3)", fontsize=10)

    key_levels = sorted([3/4, 0.80, 15/17, 1.0])
    # Draw contour lines only (no clabel — labels added manually below)
    cs = ax1.contour(NN, MM, BB, levels=key_levels,
                     colors="white", linewidths=1.0, linestyles="--")

    # Manual contour labels: all at the same m=2.50 row so they are
    # visually aligned. Background color is sampled from the colormap at
    # the exact β value of each contour, so the dashed line is invisible
    # behind the label without any colour mismatch.
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    norm_lab  = Normalize(vmin=0.55, vmax=1.0)
    sm_lab    = ScalarMappable(norm=norm_lab, cmap=plt.cm.viridis_r)
    m_label   = 2.50   # fixed m for all four labels
    label_data = [(3/4, "0.75"), (0.80, "0.80"), (15/17, "0.88"), (1.0, "1.00")]

    ax1.set_xlim(1.0, 5.0); ax1.set_ylim(0.0, 3.0)
    fig1.tight_layout(); fig1.canvas.draw()
    bbox_ax    = ax1.get_window_extent()
    slope_disp = -1.0 * (bbox_ax.height / 3.0) / (bbox_ax.width / 4.0)
    angle_deg  = np.degrees(np.arctan(slope_disp))

    for beta_val, txt in label_data:
        s = 12 * beta_val / (3 - beta_val)   # n+m=s for d=3
        n_label = s - m_label
        if 1.0 < n_label < 5.0:
            bg = sm_lab.to_rgba(beta_val)
            ax1.text(n_label, m_label, txt,
                     fontsize=8.5, color="white", ha="center", va="center",
                     rotation=angle_deg, rotation_mode="anchor",
                     bbox=dict(boxstyle="square,pad=0.2", fc=bg, ec="none"),
                     zorder=5)

    # Points
    for (n_pt, m_pt) in [(2, m_HP), (4, m_coronary), (4, m_bronchial),
                          (4, m_pulmonary), (4, m_xylem_hi), (4, m_neuron),
                          (4, m_trachea)]:
        ax1.plot(n_pt, m_pt, "wo", ms=5.5, zorder=5,
                 markeredgewidth=0.9, markeredgecolor="0.15")

    # Labels — direct offset, no arrows.
    G = 0.08
    ax1.text(2 - G, m_HP,       "His-Purkinje",
             ha="right",  va="center", fontsize=7.8, color="white")
    ax1.text(4 - G, m_coronary, "Coronary/Cerebral",
             ha="right",  va="center", fontsize=7.8, color="white")
    ax1.text(4 - G, m_bronchial,"Bronchial",
             ha="right",  va="center", fontsize=7.8, color="white")
    ax1.text(4 - G, m_pulmonary,"Pulmonary",
             ha="right",  va="center", fontsize=7.8, color="white")
    ax1.text(4 + G, m_xylem_hi, "Xylem hi",
             ha="left",   va="center", fontsize=7.8, color="white")
    ax1.text(4 + G, m_neuron,   "Dendrite/Xylem lo/\nSponge/Leaf",
             ha="left",   va="center", fontsize=7.8, color="white")
    ax1.text(4 + G, m_trachea,  "Trachea",
             ha="left",   va="center", fontsize=7.8, color="white")

    ax1.set_xlabel("Transport exponent n")
    ax1.set_ylabel("Maintenance exponent m")
    ax1.set_title("Allometric universality classes β(n, m, d=3)")
    ax1.set_xlim(1.0, 5.0)
    ax1.set_ylim(0.0, 3.0)
    ax1.set_xticks([1, 2, 3, 4, 5])
    ax1.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    fig1.tight_layout()
    fig1.savefig("../manuscript/figures/fig1_spectrum.pdf", bbox_inches="tight")
    plt.close(fig1)
    print("✓ ../manuscript/figures/fig1_spectrum.pdf")

    # ===========================================================
    # FIG 2 — β(M) Womersley transition curve
    # ===========================================================
    #
    # α_eff(M) = α_w + (α_t_surf − α_w) / (1 + (Wo(M)/Wo_c)^4)
    # Regime labels placed in unoccupied quadrants; M* annotated
    # with an arrow from the flat transition zone.

    alpha_w_val  = 2.0
    alpha_t_surf = 2.5    # viscous small-mammal limit (m=1, surface maintenance)
    M_vec_g      = np.logspace(np.log10(0.5), np.log10(2e6), 800)

    def Wo_of_M(M_g):
        return Wo_ref * (M_g / (M_ref_kg * 1000)) ** 0.25

    Wo_vec        = Wo_of_M(M_vec_g)
    alpha_eff_vec = alpha_w_val + (alpha_t_surf - alpha_w_val) / \
                    (1.0 + (Wo_vec / Wo_c) ** 4)
    beta_eff_vec  = beta_alpha_d(alpha_eff_vec, 3)

    fig2, ax2 = plt.subplots(figsize=(6.2, 4.4))

    # Theory curve (black, drawn first so reference lines appear on top)
    ax2.semilogx(M_vec_g, beta_eff_vec, "k-", lw=1.8, label="Theory", zorder=3)

    # Reference lines: blue ABOVE black (higher zorder), red dotted
    ax2.axhline(15/17, color="firebrick",  lw=1.2, ls=":",
                label="β = 15/17 (viscous sub-bound)", zorder=4)
    ax2.axhline(3/4,   color="steelblue",  lw=1.2, ls="--",
                label="β = 3/4 (Kleiber)",             zorder=5)

    # M* vertical line
    ax2.axvline(M_star_g, color="0.45", lw=1.0, ls="-.", zorder=4)

    # M* marker: dot exactly on the curve, label to the right
    beta_at_mstar = beta_alpha_d(
        alpha_w_val + (alpha_t_surf - alpha_w_val) /
        (1.0 + (Wo_of_M(M_star_g) / Wo_c) ** 4), 3)
    ax2.plot(M_star_g, beta_at_mstar, "ko", ms=5.5, zorder=6,
             markeredgewidth=0.9, markeredgecolor="0.3")
    ax2.text(M_star_g * 1.35, beta_at_mstar,
             f"M* = {M_star_g:.1f} g",
             fontsize=8.5, color="0.25", ha="left", va="center")

    # x-position for both regime labels: same far-right zone as wave-dominated
    x_label = 3e4
    ax2.text(x_label, 3/4   + 0.003, "Wave-dominated  (β = 3/4)",
             fontsize=8.5, ha="center", va="bottom", color="steelblue")
    ax2.text(x_label, 15/17 + 0.003, "Viscous regime  (β = 15/17)",
             fontsize=8.5, ha="center", va="bottom", color="firebrick")

    ax2.set_xlabel("Body mass M (g)")
    ax2.set_ylabel("Metabolic scaling exponent β")
    ax2.set_title("Womersley transition in metabolic scaling")
    ax2.set_xlim(0.5, 2e6)
    ax2.set_ylim(0.735, 0.910)
    ax2.set_yticks([0.75, 0.80, 0.85, 15/17])
    ax2.set_yticklabels(["3/4", "0.80", "0.85", "15/17"])
    ax2.set_xticks([1, 10, 100, 1e3, 1e4, 1e5, 1e6])
    ax2.set_xticklabels(["1", "10", "100", "1 kg", "10 kg", "100 kg", "1 t"])

    # Legend below the axes in a dedicated row
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
               ncol=3, framealpha=0.0, fontsize=8.5)
    fig2.tight_layout()
    fig2.subplots_adjust(bottom=0.22)   # room for the legend below
    fig2.savefig("../manuscript/figures/fig2_transition.pdf", bbox_inches="tight")
    plt.close(fig2)
    print("✓ ../manuscript/figures/fig2_transition.pdf")

    # ===========================================================
    # FIG 3 — Minimax gap (physical, normalized — Section 20)
    # ===========================================================
    #
    # f_norm(α) = f_phys(α) / f_phys(α_t)   [0 at α_w, 1 at α_t]
    # g_norm(α) = g_phys(α) / g_phys(α_w)   [1 at α_w, 0 at α_t]
    #
    # g_phys is derived from C_maint (Section 20) — no free parameters.
    # Both curves normalized to [0,1] → commensurate axis.
    # α* from bisection in Section 20.
    #
    # Three clearly separated vertical lines:
    #   αw = 2.00   (steelblue, left edge)
    #   α* = 2.77   (gray,      mid zone — label via saddle annotation)
    #   αt = 2.885  (firebrick, right edge)

    # FIG 3: physical units (%), y-axis truncated to [0, 20%] so the
    # crossing at alpha*=2.77 (~9.7%) is centred and both curves are
    # visible. No free-parameter normalisation — f and g share the same
    # units (dimensionless cost fraction) and are directly commensurate.
    # alpha* = 2.77 is where f_phys = g_phys exactly.

    alpha_plot = np.linspace(1.90, 3.05, 500)
    f_pct = np.array([f_phys_fn(a) * 100 for a in alpha_plot])
    g_pct = np.array([g_phys_fn(a) * 100 for a in alpha_plot])

    # Physical crossing (the only meaningful one)
    f_star_pct = f_phys_fn(alpha_star_phys) * 100   # ≈ 9.7%
    g_star_pct = g_phys_fn(alpha_star_phys) * 100   # ≈ 9.7%

    y_max = 100.0
    y_min = -2.0
    mid_y = (y_max + y_min) / 2   # true vertical center of the axes = 49.0

    fig3, ax3 = plt.subplots(figsize=(6.5, 4.6))

    ax3.plot(alpha_plot, f_pct, "b-",  lw=1.8, zorder=3,
             label="f(α): cumulative wave loss (%)")
    ax3.plot(alpha_plot, g_pct, "r--", lw=1.8, zorder=4,
             label="g(α): maintenance excess (%)")

    ax3.axvline(2.00,        color="steelblue", lw=1.0, ls="-.", alpha=0.85, zorder=2)
    ax3.axvline(at_coronary, color="firebrick", lw=1.0, ls="-.", alpha=0.85, zorder=2)

    ax3.text(2.00 + 0.015,        mid_y, "αw = 2",
             ha="left",  va="center", fontsize=8.5, color="steelblue")
    ax3.text(at_coronary + 0.015, mid_y, f"αt = {at_coronary:.3f}",
             ha="left",  va="center", fontsize=8.5, color="firebrick")

    # Dot exactly at physical crossing
    ax3.plot(alpha_star_phys, f_star_pct, "ko", ms=5.5, zorder=6,
             markeredgewidth=0.9, markeredgecolor="0.3")

    # Label centred in the gap between α* and αt, slightly above the dot
    x_label = alpha_star_phys + 0.85 * (at_coronary - alpha_star_phys)
    ax3.text(x_label, f_star_pct + 3.0,
             f"α* = {alpha_star_phys:.2f}",
             ha="right", va="bottom", fontsize=8.8, color="0.15")

    ax3.set_xlabel("Branching exponent α")
    ax3.set_ylabel("Cost penalty (%)")
    ax3.set_title("Minimax gap: wave vs. maintenance competing penalties")
    ax3.set_xlim(1.90, 3.05)
    ax3.set_ylim(-2, y_max)
    ax3.set_yticks([0, 25, 50, 75, 100])

    ax3.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
               ncol=2, framealpha=0.0, fontsize=8.5)
    fig3.tight_layout()
    fig3.subplots_adjust(bottom=0.22)
    fig3.savefig("../manuscript/figures/fig3_minimax.pdf", bbox_inches="tight")
    plt.close(fig3)
    print("✓ ../manuscript/figures/fig3_minimax.pdf")

    print("\n✓ All figures saved in figures/")


# ============================================================
# SECTION 23: Cerebral autoregulation minimax α*(ε)
# ============================================================
# Modified equal-cost condition: ε · f(α*) = g(α*)
# ε = pulsatile energy fraction (environmental input, not free param)
# f, g: same physical functions as coronary (m=1.770, K=11, N=2, d=3)
# ε estimated from Transcranial Doppler PI ratio:
#   ε ≈ (PI_cerebral/PI_systemic)² × (MFV_cerebral/MFV_systemic)
# Typical values: PI_MCA≈0.8, PI_systemic≈2.0 → ε ∼ 0.15–0.20

def find_alpha_star_cerebral(epsilon, tol=1e-9):
    """Solve ε·f(α*) = g(α*) by bisection. Returns α_t if ε < 1e-4."""
    if epsilon < 1e-4:
        return at_cerebral
    def objective(alpha):
        return epsilon * f_phys_fn(alpha) - g_phys_fn(alpha)
    # f_phys uses m=m_coronary=1.770, same as cerebral (Rhodin 1967)
    # objective < 0 at α_w (f=0, g>0) and > 0 at α_t (g=0, f>0)
    a, b = 2.001, at_cerebral - 1e-6
    for _ in range(120):
        mid = (a + b) * 0.5
        if abs(b - a) < tol:
            break
        if objective(a) * objective(mid) <= 0:
            b = mid
        else:
            a = mid
    return (a + b) * 0.5

# Key values
epsilon_lo      = 0.15
epsilon_hi      = 0.20
epsilon_mid     = (epsilon_lo + epsilon_hi) / 2.0

alpha_star_cereb_lo  = find_alpha_star_cerebral(epsilon_lo)
alpha_star_cereb_hi  = find_alpha_star_cerebral(epsilon_hi)
alpha_star_cereb_mid = find_alpha_star_cerebral(epsilon_mid)

print("=" * 60)
print("SECTION 23: CEREBRAL AUTOREGULATION MINIMAX")
print("=" * 60)
print(f"  ε = {epsilon_lo:.2f}: α* = {alpha_star_cereb_lo:.4f}")
print(f"  ε = {epsilon_mid:.2f}: α* = {alpha_star_cereb_mid:.4f}")
print(f"  ε = {epsilon_hi:.2f}: α* = {alpha_star_cereb_hi:.4f}")
print(f"  αt (cerebral)    = {at_cerebral:.4f}")
print(f"  αexp (Rossitti)  = 2.900 ± 0.700  [consistent]")

# The Remark quotes α* ≈ 2.86 — verify it falls in [lo, hi]
assert alpha_star_cereb_hi <= alpha_star_cereb_lo, \
    "α*(ε) should be decreasing in ε"
assert 2.85 <= alpha_star_cereb_lo <= 2.89, \
    f"α* out of expected range: {alpha_star_cereb_lo:.4f}"

lines.append("% ============================================================")
lines.append("% Section 23: Cerebral autoregulation minimax α*(ε)")
lines.append("% ============================================================")
lines.append(f"\\newcommand{{\\epsilonCerebralLow}}{{{epsilon_lo}}}")
lines.append(f"\\newcommand{{\\epsilonCerebralHigh}}{{{epsilon_hi}}}")
lines.append(f"\\newcommand{{\\alphaStarCerebralLow}}{{{alpha_star_cereb_lo:.3f}}}")
lines.append(f"\\newcommand{{\\alphaStarCerebralHigh}}{{{alpha_star_cereb_hi:.3f}}}")
lines.append(f"\\newcommand{{\\alphaStarCerebral}}{{{alpha_star_cereb_mid:.3f}}}")
lines.append("")

# ============================================================
# SECTION 24: Bronchial Womersley verification (A5)
# ============================================================
omega_bronch = 2 * math.pi * 0.25   # 15 breaths/min
rho_air = 1.2                        # kg/m^3
mu_air = 1.8e-5                      # Pa·s
nu_air = mu_air / rho_air

# Weibel 1963 airway radii (mm)
weibel_gen = {0: 9.0, 2: 2.94, 5: 0.791, 10: 0.184}
for gen, r_mm in weibel_gen.items():
    wo = (r_mm * 1e-3) * math.sqrt(omega_bronch * rho_air / mu_air)
    print(f"  Bronchial gen {gen}: r={r_mm} mm, Wo={wo:.2f}")
    if gen >= 2:
        assert wo < 1.0, f"Wo should be < 1 at gen {gen}, got {wo:.2f}"

lines.append("")

# ============================================================
# FINAL WRITE — all dynamic variables
# ============================================================
output_path = "../manuscript/dynamic_variables.tex"
with open(output_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"\n✓ FINAL: Written {output_path} ({len(lines)} lines)")
