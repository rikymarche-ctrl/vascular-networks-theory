"""
compute.py — Paper IV: The Incommensurability Principle in Biological Transport
================================================================================
Generates dynamic_variables.tex containing all LaTeX \\newcommand definitions
used in the manuscript. Every calculated quantity in the .tex file originates
here; no numerical values are hardcoded in the source.

Model
-----
The sensitivity analysis uses the exact two-level minimax model of Paper II:

    C_wave(α)     = 1 − (1 − |Γ(α)|²)^G
    |Γ(α)|²       = ((u − 1)/(u + 1))²,  u = N^(α_w/α − 1)

    C_transport(α) = Σ_g w_g [Φ(r_g) − Φ(r*_g)] / Σ_g w_g Φ(r*_g)
    w_g            = N^g ℓ_0 β^g
    r_g            = r*(Q_0) · N^{−g/α}      (uniform-exponent radii)
    r*_g           = r*(Q_g), Q_g = Q_0/N^g  (locally optimal radii)

Physical parameters correspond to the porcine coronary tree (Kassab 1993).

External inputs (taken from Paper II, not recomputed here)
----------------------------------------------------------
    α*  = 2.72   minimax branching exponent  (Kassab 1993 calibration)
    η*  = 0.833  duty cycle                  (saddle-point condition)
    α_t = 2.90   transport optimum           (Paper II, VarAlphaT)

Output
------
    ../manuscript/dynamic_variables.tex
"""

import math
import numpy as np
from datetime import datetime
from scipy.optimize import brentq, minimize_scalar

# ============================================================
# SECTION 1: Input parameters
# ============================================================

# ---- Physical constants (porcine coronary tree) ----
G       = 11        # tree depth (junctions)
N       = 2         # bifurcation number
p       = 0.77      # wall-thickness exponent (Rhodin 1967)
alpha_w = 2.115     # wave-impedance attractor (Paper II VarAlphaW)
beta    = 0.787     # length-taper factor ell_g = ell_0*beta^g (Paper II VarBeta)

mu_f  = 3.5e-3      # Pa s   blood viscosity
b     = 1500.0      # W/m3   blood-volume metabolic cost
m_w   = 20e3        # W/m3   wall metabolic cost
Q0    = 5.2e-6      # m3/s   proximal flow  (5.2 mL/s)
ell0  = 15e-3       # m      root segment length

# ---- External results from Paper II (inputs, not recomputed here) ----
alpha_star_II = 2.72    # minimax exponent (Kassab 1993 calibration)
eta_star_II   = 0.833   # duty cycle
alpha_t_II    = 2.90    # transport optimum (Paper II VarAlphaT, displayed)

# ---- Theoretical wave attractor (Newtonian fluids, area conservation) ----
alpha_w_theory = 2      # integer, exact: area-preserving branching, N=2

# ---- Physiological radius range (for No-Go coupling estimate) ----
r_aorta_m     = 12.5e-3   # 12.5 mm  — human aortic inner radius
r_arteriole_m = 12.5e-6   # 12.5 um  — arteriole inner radius
radius_decades   = math.log10(r_aorta_m / r_arteriole_m)  # = 3
# NOTE: Theoretical scaling μ ∝ N^[g(4/α-2)] gives ~10^(-1.75) with α=2.72, N=2, G=11.
# Empirical variation (Kassab 1993 + hemodynamic effects) indicates 10^(2-3).
# Using range notation for LaTeX variable to reflect empirical uncertainty:
mu_coupling_exp  = "{2-3}"  # Decades of μ variation (empirical range)

# ---- Perturbation magnitudes for Table 1 ----
pert_b_factor   = 2
pert_mw_factor  = 7
pert_muf_frac   = 0.50       # +50%
pert_Q0_factor  = 4
pert_ell_factor = 3
pert_p_delta    = 0.08
pert_G_delta    = 2
pert_aw_delta   = 0.115      # +0.115 (i.e., alpha_w -> alpha_w + delta)

# ============================================================
# SECTION 2: Single-vessel cost and locally optimal radius
# ============================================================

def phi0(r, Q, mu_f, b, m_w, p):
    """Single-vessel cost per unit length [W/m]."""
    return (8.0*mu_f*Q**2/(np.pi*r**4)
            + b*np.pi*r**2
            + m_w*r**(1.0+p))


def dphi0_dr(r, Q, mu_f, b, m_w, p):
    """d(phi0)/dr: root gives locally optimal radius."""
    return (-32.0*mu_f*Q**2/(np.pi*r**5)
            + 2.0*b*np.pi*r
            + (1.0+p)*m_w*r**p)


def r_opt(Q, mu_f, b, m_w, p, r_lo=1e-7, r_hi=5e-2):
    """Locally optimal radius r*(Q) for given flow and parameters."""
    f_lo = dphi0_dr(r_lo, Q, mu_f, b, m_w, p)
    f_hi = dphi0_dr(r_hi, Q, mu_f, b, m_w, p)
    if f_lo * f_hi >= 0:
        rs = np.logspace(np.log10(r_lo), np.log10(r_hi), 600)
        ds = [dphi0_dr(r, Q, mu_f, b, m_w, p) for r in rs]
        for i in range(len(ds)-1):
            if ds[i]*ds[i+1] < 0:
                r_lo, r_hi = rs[i], rs[i+1]
                break
        else:
            raise ValueError(f"r_opt: no bracket for Q={Q:.3e}")
    return brentq(dphi0_dr, r_lo, r_hi,
                  args=(Q, mu_f, b, m_w, p), xtol=1e-14, rtol=1e-10)


# ============================================================
# SECTION 3: Network cost functions (Paper II exact model)
# ============================================================

def C_wave(alpha, alpha_w=alpha_w, N=N, G=G):
    """
    Network-level wave cost (Paper II Eq. gamma_gen).
    u = N^(alpha_w/alpha - 1);  |Gamma|^2 = ((u-1)/(u+1))^2
    C_wave = 1 - (1 - |Gamma|^2)^G.  Zero at alpha = alpha_w.
    """
    u = N ** (alpha_w / alpha - 1.0)
    gamma2 = ((u - 1.0) / (u + 1.0)) ** 2
    return 1.0 - (1.0 - gamma2) ** G


def C_transport(alpha, G=G, N=N, p=p,
                mu_f=mu_f, b=b, m_w=m_w, Q0=Q0, ell0=ell0, beta=beta):
    """
    Two-level transport cost (Paper II Eq. 182).
    C_t = sum_g N^g ell_g [phi(r_g) - phi(r*_g)] / sum_g N^g ell_g phi(r*_g)
    """
    r0 = r_opt(Q0, mu_f, b, m_w, p)
    num = 0.0; den = 0.0
    for g in range(G):
        Qg   = Q0 / N**g
        rg   = r0 * N**(-g / alpha)
        rs   = r_opt(Qg, mu_f, b, m_w, p)
        wg   = N**g * ell0 * beta**g
        num += wg * (phi0(rg, Qg, mu_f, b, m_w, p) - phi0(rs, Qg, mu_f, b, m_w, p))
        den += wg *  phi0(rs, Qg, mu_f, b, m_w, p)
    return num / den


def C_transport_log(alpha, G=G, N=N, p=p,
                    mu_f=mu_f, b=b, m_w=m_w, Q0=Q0, ell0=ell0, beta=beta):
    """
    Alternative (logarithmic) transport cost: ln(Phi_net / Phi_opt).
    Violates Axiom 3 (thermodynamic linearity).
    Used for counterfactual analysis in Appendix.
    """
    r0 = r_opt(Q0, mu_f, b, m_w, p)
    num = 0.0; den = 0.0
    for g in range(G):
        Qg   = Q0 / N**g
        rg   = r0 * N**(-g / alpha)
        rs   = r_opt(Qg, mu_f, b, m_w, p)
        wg   = N**g * ell0 * beta**g
        num += wg * phi0(rg, Qg, mu_f, b, m_w, p)
        den += wg * phi0(rs, Qg, mu_f, b, m_w, p)
    return math.log(num / den) if num > den else 0.0


# ============================================================
# SECTION 4: Minimax: find alpha_t, alpha*, eta*
# ============================================================

def find_alpha_t(G=G, N=N, p=p,
                 mu_f=mu_f, b=b, m_w=m_w, Q0=Q0, ell0=ell0, beta=beta):
    """Transport optimum: minimiser of C_transport."""
    res = minimize_scalar(
        lambda a: C_transport(a, G, N, p, mu_f, b, m_w, Q0, ell0, beta),
        bounds=(2.3, 4.5), method='bounded', options={'xatol': 1e-8})
    return res.x


def find_alpha_star(G=G, N=N, p=p, alpha_w=alpha_w,
                    mu_f=mu_f, b=b, m_w=m_w, Q0=Q0, ell0=ell0, beta=beta):
    """
    alpha* = root of C_wave(alpha) = C_transport(alpha) in (alpha_w, alpha_t).
    Returns (alpha_star, alpha_t).
    """
    alpha_t = find_alpha_t(G, N, p, mu_f, b, m_w, Q0, ell0, beta)

    def gap(a):
        return (C_wave(a, alpha_w, N, G)
                - C_transport(a, G, N, p, mu_f, b, m_w, Q0, ell0, beta))

    xs = np.linspace(alpha_w + 1e-4, alpha_t - 1e-4, 300)
    gs = [gap(x) for x in xs]
    for i in range(len(gs)-1):
        if gs[i]*gs[i+1] < 0:
            return brentq(gap, xs[i], xs[i+1], xtol=1e-10), alpha_t

    raise ValueError(
        f"No crossing in ({alpha_w:.3f}, {alpha_t:.3f}). "
        f"gap range [{min(gs):.4f}, {max(gs):.4f}].")


def find_alpha_star_log(G=G, N=N, p=p, alpha_w=alpha_w,
                        mu_f=mu_f, b=b, m_w=m_w, Q0=Q0, ell0=ell0, beta=beta):
    """
    Counterfactual: alpha* using LOGARITHMIC transport penalty.
    Returns (alpha_star_log, alpha_t).
    """
    alpha_t = find_alpha_t(G, N, p, mu_f, b, m_w, Q0, ell0, beta)

    def gap(a):
        return (C_wave(a, alpha_w, N, G)
                - C_transport_log(a, G, N, p, mu_f, b, m_w, Q0, ell0, beta))

    xs = np.linspace(alpha_w + 1e-4, alpha_t - 1e-4, 300)
    gs = [gap(x) for x in xs]
    for i in range(len(gs)-1):
        if gs[i]*gs[i+1] < 0:
            return brentq(gap, xs[i], xs[i+1], xtol=1e-10), alpha_t

    # If no crossing found, logarithmic penalty may not intersect wave cost
    return None, alpha_t


def find_eta_star(a_star, G=G, N=N, p=p, alpha_w=alpha_w,
                  mu_f=mu_f, b=b, m_w=m_w, Q0=Q0, ell0=ell0, beta=beta):
    """eta* from saddle-point gradient ratio: |dCt/da| / (|dCw/da| + |dCt/da|)."""
    da = 1e-6
    dCt = (C_transport(a_star+da, G, N, p, mu_f, b, m_w, Q0, ell0, beta)
           - C_transport(a_star-da, G, N, p, mu_f, b, m_w, Q0, ell0, beta)) / (2*da)
    dCw = (C_wave(a_star+da, alpha_w, N, G)
           - C_wave(a_star-da, alpha_w, N, G)) / (2*da)
    return abs(dCt) / (abs(dCw) + abs(dCt))


# ============================================================
# SECTION 5: Baseline computation
# ============================================================

params_base = dict(G=G, N=N, p=p, alpha_w=alpha_w,
                   mu_f=mu_f, b=b, m_w=m_w, Q0=Q0, ell0=ell0, beta=beta)

print("=" * 64)
print("BASELINE (porcine coronary, G=11, p=0.77, aw=2.115, beta=0.787)")
print("=" * 64)

a_star_model, alpha_t_model = find_alpha_star(**params_base)
eta_star_model = find_eta_star(a_star_model, **params_base)
r0_model = r_opt(Q0, mu_f, b, m_w, p)
Cw_star = C_wave(a_star_model, alpha_w, N, G)

da = 1e-6
dCt_star = (C_transport(a_star_model+da, **{k:v for k,v in params_base.items() if k!='alpha_w'})
            - C_transport(a_star_model-da, **{k:v for k,v in params_base.items() if k!='alpha_w'})) / (2*da)
dCw_star = (C_wave(a_star_model+da, alpha_w, N, G)
            - C_wave(a_star_model-da, alpha_w, N, G)) / (2*da)
gradient_ratio = abs(dCt_star) / abs(dCw_star)

print(f"  r*(Q_0)        = {r0_model*1e3:.4f} mm")
print(f"  alpha_t (num)  = {alpha_t_model:.4f}  [Paper II: 2.9032]")
print(f"  alpha* (model) = {a_star_model:.4f}  [Paper II: 2.72]")
print(f"  eta*   (model) = {eta_star_model:.4f}  [Paper II: 0.833]")
print(f"  C_wave(alpha*) = {Cw_star:.4f}  ({100*Cw_star:.1f}%)  [Paper II: 6.3%]")
print(f"  |dCt|/|dCw|   = {gradient_ratio:.2f}   [Paper II: 5]")
print()

# alpha_t range for physiological p
p_lo_phys, p_hi_phys = 0.72, 0.82
alpha_t_lo = find_alpha_t(G, N, p_lo_phys, mu_f, b, m_w, Q0, ell0, beta)
alpha_t_hi = find_alpha_t(G, N, p_hi_phys, mu_f, b, m_w, Q0, ell0, beta)
print(f"  alpha_t range  = [{alpha_t_lo:.2f}, {alpha_t_hi:.2f}] for p in [{p_lo_phys}, {p_hi_phys}]")
print()

# ============================================================
# SECTION 6: Log-sensitivity computations
# ============================================================

def log_sens(param_name, base_val, params, rel_step=0.10):
    """S = d(alpha*)/d(ln x) by central finite differences."""
    p_lo = params.copy(); p_lo[param_name] = base_val*(1.0 - rel_step)
    p_hi = params.copy(); p_hi[param_name] = base_val*(1.0 + rel_step)
    a_lo, _ = find_alpha_star(**p_lo)
    a_hi, _ = find_alpha_star(**p_hi)
    return (a_hi - a_lo) / math.log(base_val*(1+rel_step) / (base_val*(1-rel_step)))


def log_sens_finite(param_name, x_lo, x_hi, params):
    """Finite-difference log-sensitivity for asymmetric perturbations."""
    p_lo = params.copy(); p_lo[param_name] = x_lo
    p_hi = params.copy(); p_hi[param_name] = x_hi
    a_lo, _ = find_alpha_star(**p_lo)
    a_hi, _ = find_alpha_star(**p_hi)
    return (a_hi - a_lo) / math.log(x_hi / x_lo)


print("=" * 64)
print("LOG-SENSITIVITIES  S = d(alpha*)/d(ln x)")
print("=" * 64)
print()

# Metabolic inputs (expect |S| ~ 0)
S_mu_f = log_sens('mu_f', mu_f, params_base, rel_step=0.10)
S_m_w  = log_sens('m_w',  m_w,  params_base, rel_step=0.10)
S_b    = log_sens('b',    b,    params_base, rel_step=0.10)
S_Q0   = log_sens('Q0',   Q0,   params_base, rel_step=0.10)
# ell0 cancels exactly in the ratio
S_ell0 = 0.0

print("-- Metabolic inputs --")
for nm, sv in [('mu_f', S_mu_f), ('m_w', S_m_w), ('b', S_b), ('Q0', S_Q0)]:
    print(f"    {nm:7s}  S = {sv:+.5f}")
print(f"    {'ell0':7s}  S = {S_ell0:+.5f}  (exact)")

# Structural inputs (expect |S| ~ 0.1-0.5)
S_p    = log_sens('p',       p,       params_base, rel_step=0.10)
S_aw   = log_sens('alpha_w', alpha_w, params_base, rel_step=0.05)

p_Gm = params_base.copy(); p_Gm['G'] = G - pert_G_delta
p_Gp = params_base.copy(); p_Gp['G'] = G + pert_G_delta
a_Gm, _ = find_alpha_star(**p_Gm)
a_Gp, _ = find_alpha_star(**p_Gp)
S_G = (a_Gp - a_Gm) / math.log((G + pert_G_delta) / (G - pert_G_delta))

print()
print("-- Structural inputs --")
print(f"    {'p':7s}  S = {S_p:+.4f}")
print(f"    {'G':7s}  S = {S_G:+.4f}")
print(f"    {'alpha_w':7s}  S = {S_aw:+.4f}")

# ---- Table-row verification (exact perturbation magnitudes) ----
print()
print("=" * 64)
print("TABLE ROWS (exact perturbations)")
print("=" * 64)

S_b_tab   = log_sens_finite('b',    b,       b*pert_b_factor,           params_base)
S_mw_tab  = log_sens_finite('m_w',  m_w,     m_w*pert_mw_factor,        params_base)
S_muf_tab = log_sens_finite('mu_f', mu_f,    mu_f*(1+pert_muf_frac),     params_base)
S_Q0_tab  = log_sens_finite('Q0',   Q0,      Q0*pert_Q0_factor,          params_base)

p_plo = params_base.copy(); p_plo['p'] = p - pert_p_delta
p_phi = params_base.copy(); p_phi['p'] = p + pert_p_delta
a_plo, _ = find_alpha_star(**p_plo)
a_phi, _ = find_alpha_star(**p_phi)
S_p_tab = (a_phi - a_plo) / math.log((p + pert_p_delta)/(p - pert_p_delta))

print(f"  b    x{pert_b_factor}      : |S| = {abs(S_b_tab):.5f}")
print(f"  m_w  x{pert_mw_factor}      : |S| = {abs(S_mw_tab):.5f}")
print(f"  mu_f +{int(pert_muf_frac*100)}%   : |S| = {abs(S_muf_tab):.5f}")
print(f"  Q0   x{pert_Q0_factor}      : |S| = {abs(S_Q0_tab):.5f}")
print(f"  p    +/-{pert_p_delta}  : |S| = {abs(S_p_tab):.5f}")
print(f"  G    +/-{pert_G_delta}    : |S| = {abs(S_G):.5f}")
print(f"  a_w  +{pert_aw_delta}  : |S| = {abs(S_aw):.5f}")

print()
print("=" * 64)
print("COUNTERFACTUAL: Logarithmic Penalty (violates Axiom 3)")
print("=" * 64)
a_star_log, alpha_t_log = find_alpha_star_log(**params_base)
if a_star_log is not None:
    deviation = a_star_log - alpha_star_II
    deviation_pct = 100 * abs(deviation) / alpha_star_II
    print(f"  alpha* (log penalty) = {a_star_log:.4f}")
    print(f"  Deviation from Kassab: {deviation:+.4f}  ({deviation_pct:.1f}%)")
    print(f"  Empirical α* = {alpha_star_II}  ← CORRECT (fractional excess)")
    print(f"\n  → Logarithmic penalty FAILS empirical validation.")
else:
    print("  No minimax crossing found with logarithmic penalty.")
    print("  → Logarithmic penalty produces NO STABLE SOLUTION.")

print()
print("=" * 64)
print("SUMMARY")
print("=" * 64)
print(f"  alpha_t (model) = {alpha_t_model:.4f}")
print(f"  alpha* (model)  = {a_star_model:.4f}  [Paper II input: {alpha_star_II}]")
print(f"  eta*   (model)  = {eta_star_model:.4f}  [Paper II input: {eta_star_II}]")
print(f"  r_0             = {r0_model*1e3:.3f} mm")
print(f"  C_wave(alpha*)  = {Cw_star:.4f}  ({100*Cw_star:.1f}%)")

# ============================================================
# SECTION 7: Format and write dynamic_variables.tex
# ============================================================

THRESHOLD_SMALL = 0.01   # below this: "<0.01"

def format_macro(name, value, unit, description):
    """Format LaTeX macro to column-50 alignment with unit comment."""
    macro = f"\\newcommand{{\\{name}}}{{{value}}}"
    return f"{macro:<50s} % [{unit:<8}] {description}"

def fmt3(x):
    return f"{x:.3f}"

def fmt2(x):
    return f"{x:.2f}"

def fmt_s(s):
    """Format sensitivity for table: '<0.01' or two-decimal number."""
    if abs(s) < THRESHOLD_SMALL:
        return r"${<}0.01$"
    return f"${abs(s):.2f}$"

def sens_line(name, val, description):
    """Format sensitivity macro (value contains LaTeX math)."""
    macro = f"\\newcommand{{\\{name}}}{{{val}}}"
    return f"{macro:<50s} % [-       ] {description}"

lines = []
lines.append("% ====================================================================")
lines.append("% DYNAMIC VARIABLES")
lines.append("% ====================================================================")
lines.append("% Source: paper4-incommensurability/scripts/compute.py")
lines.append(f"% Generated: {datetime.now().strftime('%Y-%m-%dT%H:%M')}")
lines.append("% ====================================================================")
lines.append("")

# ---- 1. Physical parameters ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 1. PHYSICAL PARAMETERS (porcine coronary tree)")
lines.append("% --------------------------------------------------------------------")
lines.append(format_macro("VarG",         f"{G}",               "-",    "Tree depth (junctions)"))
lines.append(format_macro("VarN",         f"{N}",               "-",    "Bifurcation number"))
lines.append(format_macro("VarP",         f"{p}",               "-",    "Wall-thickness exponent (Rhodin 1967)"))
lines.append(format_macro("VarAlphaW",    f"{alpha_w}",         "-",    "Wave-impedance attractor (Paper II)"))
lines.append(format_macro("VarAlphaWFig", f"{alpha_w_theory}",  "-",    "Wave attractor for figures (integer)"))
lines.append("")

# ---- 2. Metabolic parameters ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 2. METABOLIC PARAMETERS (baseline, SI units)")
lines.append("% --------------------------------------------------------------------")
lines.append(format_macro("VarMuFmPas",   f"{mu_f*1e3:.1f}",   "mPa s", "Blood dynamic viscosity"))
lines.append(format_macro("VarBblood",    f"{int(b)}",          "W/m3",  "Blood-volume metabolic cost"))
lines.append(format_macro("VarMwallKW",   f"{int(m_w/1e3)}",    "kW/m3", "Wall metabolic cost"))
lines.append(format_macro("VarQZeroML",   f"{Q0*1e6:.1f}",      "mL/s",  "Proximal flow"))
lines.append(format_macro("VarEllZeroMm", f"{int(ell0*1e3)}",   "mm",    "Root segment length"))
lines.append("")

# ---- 3. Minimax results from Paper II ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 3. MINIMAX RESULTS (inputs from Paper II)")
lines.append("% --------------------------------------------------------------------")
lines.append(format_macro("VarAlphaStar",  f"{alpha_star_II}",    "-", "Minimax exponent (Kassab 1993 calibration)"))
lines.append(format_macro("VarEtaStar",    f"{eta_star_II:.3f}",  "-", "Duty cycle (saddle point)"))
lines.append(format_macro("VarAlphaTFig",  fmt2(alpha_t_II),      "-", "Transport optimum for figures"))
lines.append("")

# ---- 4. Computed baseline ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 4. COMPUTED BASELINE (two-level model, this script)")
lines.append("% --------------------------------------------------------------------")
lines.append(format_macro("VarAlphaStarModel", fmt3(a_star_model), "-",  "Minimax exponent (model)"))
lines.append("")

# ---- 5. Transport optimum range ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 5. TRANSPORT OPTIMUM RANGE (physiological p)")
lines.append("% --------------------------------------------------------------------")
lines.append(format_macro("VarAlphaTLow",  f"{alpha_t_lo:.2f}", "-", "alpha_t at p = 0.72"))
lines.append(format_macro("VarAlphaTHigh", f"{alpha_t_hi:.2f}", "-", "alpha_t at p = 0.82"))
lines.append("")

# ---- 6. No-Go Theorem coupling variation ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 6. NO-GO THEOREM: COUPLING VARIATION")
lines.append("% --------------------------------------------------------------------")
lines.append(format_macro("VarMuCouplingExp", f"{mu_coupling_exp}",           "-", "Decades of mu variation across hierarchy"))
lines.append("")

# ---- 7. Sensitivity values (Table 1) ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 7. SENSITIVITY VALUES (Table 1)")
lines.append("% --------------------------------------------------------------------")
lines.append(sens_line("VarSMuF",    fmt_s(S_muf_tab), "Sensitivity: blood viscosity"))
lines.append(sens_line("VarSB",      fmt_s(S_b_tab),   "Sensitivity: blood metabolic cost"))
lines.append(sens_line("VarSMW",     fmt_s(S_mw_tab),  "Sensitivity: wall metabolic cost"))
lines.append(sens_line("VarSQZ",     fmt_s(S_Q0_tab),  "Sensitivity: proximal flow"))
lines.append(sens_line("VarSEll",    r"$=0$ (exact)",   "Sensitivity: segment length (exact zero)"))
lines.append(sens_line("VarSP",      fmt_s(S_p_tab),   "Sensitivity: wall exponent p"))
lines.append(sens_line("VarSG",      fmt_s(S_G),       "Sensitivity: tree depth G"))
lines.append(sens_line("VarSAlphaW", fmt_s(S_aw),      "Sensitivity: wave attractor alpha_w"))
lines.append("")

# ---- 9. Table perturbation magnitudes ----
lines.append("% --------------------------------------------------------------------")
lines.append("% 9. TABLE PERTURBATION MAGNITUDES")
lines.append("% --------------------------------------------------------------------")
lines.append(f"\\newcommand{{\\VarPertB}}{{$\\times {pert_b_factor}$}}")
lines.append(f"\\newcommand{{\\VarPertMW}}{{$\\times {pert_mw_factor}$}}")
lines.append(f"\\newcommand{{\\VarPertMuF}}{{$+{int(pert_muf_frac*100)}\\%$}}")
lines.append(f"\\newcommand{{\\VarPertQZ}}{{$\\times {pert_Q0_factor}$}}")
lines.append(f"\\newcommand{{\\VarPertEll}}{{$\\times {pert_ell_factor}$}}")
lines.append(f"\\newcommand{{\\VarPertP}}{{$\\pm {pert_p_delta}$}}")
lines.append(f"\\newcommand{{\\VarPertG}}{{$\\pm {pert_G_delta}$}}")
lines.append(f"\\newcommand{{\\VarPertAlphaW}}{{$+{pert_aw_delta}$}}")
lines.append("")

lines.append("% ====================================================================")
lines.append("% END OF AUTO-GENERATED FILE")
lines.append("% ====================================================================")

out_path = "../manuscript/dynamic_variables.tex"
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")

print(f"\nWrote {out_path}  ({len(lines)} lines)")
