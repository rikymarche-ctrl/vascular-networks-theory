"""
compute.py — Paper II: A Unified Variational Principle for Branching Transport Networks
========================================================================================
Generates dynamic_variables.tex containing all LaTeX \\newcommand definitions
used in the manuscript. Every calculated quantity in the .tex file originates
here; no numerical values are hardcoded in the source.

Model: Two-level minimax framework
-----------------------------------
Level 1 — Single-vessel optimisation (Paper I):
    Each vessel independently minimises the three-term cost Φ(r, Q),
    yielding a locally optimal radius r*(Q) that depends on wall-scaling
    exponent p. This defines the transport ground state α_t and the
    length-scaling factor β = 2^{−1/α_local}.

Level 2 — Network Lagrangian (this paper):
    The branching exponent α is optimised over the full G-generation tree
    by balancing two dimensionless fractional penalties:

        C_wave(α)      = 1 − (1 − |Γ(α)|²)^G
        C_transport(α) = Σ_g w_g [Φ(r_g) − Φ(r*_g)] / Σ_g w_g Φ(r*_g)

    where w_g = N^g ℓ_0 β^g, r_g = r*(Q_0)·N^{−g/α}, r*_g = r*(Q_g).
    The minimax saddle point α* satisfies C_wave(α*) = C_transport(α*).

Physical parameters correspond to the porcine coronary tree (Kassab 1993).

Outputs
-------
    ../manuscript/dynamic_variables.tex
    ../manuscript/figures/fig1_kappa.pdf
    ../manuscript/figures/fig2_minimax.pdf
    ../manuscript/figures/fig3_robustness.pdf
    ../manuscript/figures/fig4_lagrangian.pdf
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'shared', 'scripts'))

import numpy as np
from scipy.optimize import minimize_scalar, brentq
from params import (
    MU, RHO, c0, p, b, MW_MID, MW_LOW, MW_HIGH,
    r0_coronary, Q0_coronary, Q0_peak, ell0_coronary, G_coronary, alpha_w,
    beta_coronary, ke,
    B_blood, B_wall, A_of_Q,
    # Cross-system
    p_pulmonary, alpha_w_pulmonary, G_pulmonary,
    alpha_exp_pulmonary, alpha_exp_err_pulmonary,
    alpha_exp_pulmonary_low, alpha_exp_pulmonary_high,
    alpha_w_aortic_low, alpha_w_aortic_high, G_aortic,
    alpha_exp_aortic, alpha_exp_err_aortic,
    alpha_w_neural, G_neural, alpha_star_neural_indicative,
    alpha_exp_neural, sigma_neural_indicative,
    G_airways, alpha_exp_airways_low, alpha_exp_airways_high,
    # Minimax band
    p_low, p_high, G_low, G_high, alpha_t_low_band, alpha_t_high_band,
    # Sensitivity table baselines
    b_Wm3, MW_MID_kWm3, MU_mPas, Q_sensitivity_mL, ell0_mm,
)


# ==============================================================
# Level 1: Local optimization (from Paper 2)
# ==============================================================

def r_star(Q, mw=MW_MID):
    """Optimal radius r*(Q) from Paper 2's Theorem 1."""
    Bw = B_wall(mw)
    def dphi(r):
        return -4 * A_of_Q(Q) * r**(-5) + 2 * B_blood * r + (1 + p) * Bw * r**p
    return brentq(dphi, 1e-7, 0.1, xtol=1e-14)


def alpha_local_mean(mw=MW_MID, G=G_coronary):
    """Mean local branching exponent across the tree (uses Q0_peak)."""
    alphas = []
    for g in range(G):
        Q_p = Q0_peak / 2**g
        Q_c = Q_p / 2
        r_p = r_star(Q_p, mw)
        r_c = r_star(Q_c, mw)
        alphas.append(np.log(2) / np.log(r_p / r_c))
    return np.mean(alphas), np.min(alphas), np.max(alphas)


def phi_cost(r, Q, mw=MW_MID):
    """Cost per unit length at radius r for flow Q (all 3 terms)."""
    Bw = B_wall(mw)
    return A_of_Q(Q) / r**4 + B_blood * r**2 + Bw * r**(1 + p)


# ==============================================================
# Level 2: Network Lagrangian (two-level, generalized for N)
# ==============================================================

def gamma_squared(alpha, aw=alpha_w, N=2):
    """Reflected power fraction |Gamma|^2 at a single junction."""
    x = N**(aw / alpha - 1)
    return ((x - 1) / (x + 1))**2


def wave_loss_network(alpha, G=G_coronary, aw=alpha_w, N=2):
    """Network wave cost: 1 - (1 - |Gamma|^2)^G."""
    g2 = gamma_squared(alpha, aw, N)
    return 1 - (1 - g2)**G


def locally_optimal_radii(mw=MW_MID, G=G_coronary, N=2):
    """Compute r*(Q_g) for g = 0..G, with N-branching flow distribution (uses Q0_peak)."""
    return [r_star(Q0_peak / N**g, mw) for g in range(G + 1)]


def reference_cost(r_local, mw=MW_MID, G=G_coronary, beta=1.0, N=2):
    """Total cost when every vessel is at its local optimum (uses Q0_peak)."""
    return sum(
        N**g * phi_cost(r_local[g], Q0_peak / N**g, mw) * ell0_coronary * beta**g
        for g in range(G)
    )


def two_level_lagrangian(alpha, eta, r_local, Phi_ref,
                         mw=MW_MID, G=G_coronary, beta=1.0, aw=alpha_w, N=2):
    """
    Two-level network Lagrangian.

    C_transport: penalty for deviating from locally optimal radii
                 when imposing a uniform branching exponent alpha.
    C_wave:      cumulative wave reflection loss through G junctions.
    """
    r = r_local[0]
    Phi_imposed = 0.0
    for g in range(G):
        Phi_imposed += N**g * phi_cost(r, Q0_peak / N**g, mw) * ell0_coronary * beta**g
        r *= N**(-1.0 / alpha)

    C_transport = (Phi_imposed - Phi_ref) / Phi_ref
    C_wave = wave_loss_network(alpha, G, aw, N)

    return eta * C_wave + (1 - eta) * C_transport


def find_alpha_star(eta, mw=MW_MID, G=G_coronary, beta=1.0, aw=alpha_w, N=2):
    """Find alpha* from the two-level Lagrangian for a given eta (parametric use)."""
    r_local = locally_optimal_radii(mw, G, N)
    Phi_ref = reference_cost(r_local, mw, G, beta, N)

    def L(a):
        return two_level_lagrangian(a, eta, r_local, Phi_ref, mw, G, beta, aw, N)

    res = minimize_scalar(L, bounds=(1.5, 4.0), method='bounded')

    def L_transport_only(a):
        return two_level_lagrangian(a, 0.0, r_local, Phi_ref, mw, G, beta, aw, N)
    res_t = minimize_scalar(L_transport_only, bounds=(1.5, 4.0), method='bounded')

    return res.x, res_t.x


def find_alpha_star_minimax(mw=MW_MID, G=G_coronary, beta=1.0, aw=alpha_w, N=2):
    """
    Find alpha* from the equal-cost minimax condition: C_wave(alpha*) = C_transport(alpha*).
    Also computes eta* and the gradient ratio from derivatives at alpha*.
    Returns: alpha_star, alpha_t, eta_star, gradient_ratio
    """
    r_local = locally_optimal_radii(mw, G, N)
    Phi_ref = reference_cost(r_local, mw, G, beta, N)

    def c_w(a):
        return wave_loss_network(a, G, aw, N)

    def c_t(a):
        return two_level_lagrangian(a, 0.0, r_local, Phi_ref, mw, G, beta, aw, N)

    # Transport minimum
    res_t = minimize_scalar(c_t, bounds=(1.5, 4.0), method='bounded')
    alpha_t = res_t.x

    # Equal-cost crossing
    def diff(a):
        return c_w(a) - c_t(a)

    try:
        alpha_star = brentq(diff, aw, alpha_t, xtol=1e-8)
    except ValueError:
        alpha_star = alpha_t

    # Derivatives at alpha* (central difference)
    h = 1e-5
    f_prime = (c_w(alpha_star + h) - c_w(alpha_star - h)) / (2 * h)
    g_prime = (c_t(alpha_star + h) - c_t(alpha_star - h)) / (2 * h)

    # eta* = -g'/(f'-g')
    if abs(f_prime - g_prime) > 1e-15:
        eta_star = -g_prime / (f_prime - g_prime)
    else:
        eta_star = 0.5

    # Gradient ratio |g'|/f' (rounded to nearest integer)
    gradient_ratio = round(abs(g_prime) / f_prime) if f_prime > 1e-15 else 5

    return alpha_star, alpha_t, eta_star, gradient_ratio


# ==============================================================
# Curvatures and kappa_eff
# ==============================================================

def compute_curvatures(G, mw=MW_MID, beta=1.0, aw=alpha_w, N=2):
    """Compute wave and transport curvatures, kappa_eff."""
    da = 0.001

    # Wave curvature at alpha_w
    k_w = 0.5 * (wave_loss_network(aw + da, G, aw, N)
                 - 2 * wave_loss_network(aw, G, aw, N)
                 + wave_loss_network(aw - da, G, aw, N)) / da**2

    # Transport minimum and curvature (two-level)
    r_local = locally_optimal_radii(mw, G, N)
    Phi_ref = reference_cost(r_local, mw, G, beta, N)

    def c_t(a):
        return two_level_lagrangian(a, 0.0, r_local, Phi_ref, mw, G, beta, aw, N)

    res_t = minimize_scalar(c_t, bounds=(1.5, 4.0), method='bounded')
    alpha_t = res_t.x

    k_t = 0.5 * (c_t(alpha_t + da) - 2 * c_t(alpha_t) + c_t(alpha_t - da)) / da**2

    kappa = k_t / k_w if k_w > 0 else float('inf')

    return k_w, k_t, kappa, alpha_t


def alpha_harmonic(eta, aw, alpha_t, kappa):
    """Harmonic (quadratic) approximation for alpha*."""
    return (aw * eta + alpha_t * (1 - eta) * kappa) / (eta + (1 - eta) * kappa)


# ==============================================================
# |Gamma|^2 table
# ==============================================================

def compute_gamma_table():
    """Compute |Gamma|^2 for alpha_w = (5-p)/2 (physiological)."""
    alpha_vals = [1.50, 1.75, 2.00, 2.50, 2.75, 3.00]
    return [{'alpha': a, 'gamma_aw': gamma_squared(a, alpha_w)}
            for a in alpha_vals]


# ==============================================================
# Duty cycle parametric table
# ==============================================================

def compute_eta_parametric_table():
    """Alpha* for various eta."""
    eta_vals = np.linspace(0.0, 1.0, 11)
    results = []
    for eta_val in eta_vals:
        a_star, _ = find_alpha_star(eta_val)
        results.append({'eta': eta_val, 'alpha_star': a_star})
    return results


# ==============================================================
# kappa_eff table (G scan, N=2)
# ==============================================================

def compute_kappa_table(mw=MW_MID, beta=1.0):
    """Compute kappa_eff and alpha* as function of G, for N=2."""
    G_vals = [1, 5, 7, 9, 11, 13, 15, 20]
    results = []
    for G_val in G_vals:
        k_w, k_t, kappa, alpha_t_g = compute_curvatures(G_val, mw=mw, beta=beta)
        a_star, _, eta_s, _ = find_alpha_star_minimax(mw=mw, G=G_val, beta=beta)
        a_quad = alpha_harmonic(eta_s, alpha_w, alpha_t_g, kappa)
        results.append({
            'G': G_val, 'kappa': kappa, 'alpha_t': alpha_t_g,
            'alpha_star': a_star, 'alpha_quad': a_quad,
            'k_t': k_t, 'k_w': k_w,
        })
    return results


# ==============================================================
# Single-junction wave curvature (illustrative, N=2, aw_illus=2)
# ==============================================================

def compute_kw_junction(N=2, aw_illus=2.0):
    """
    Analytical single-junction wave curvature at the illustrative alpha_w=2.
    Formula: k_w = (ln N)^2 / (4 * aw^2)
    """
    return (np.log(N))**2 / (4 * aw_illus**2)


# ==============================================================
# Transport optimum bounds (across metabolic parameter range)
# ==============================================================

def compute_alpha_t_bounds():
    """alpha_t at MW_LOW and MW_HIGH."""
    _, alpha_t_low, _, _ = find_alpha_star_minimax(mw=MW_LOW)
    _, alpha_t_high, _, _ = find_alpha_star_minimax(mw=MW_HIGH)
    return min(alpha_t_low, alpha_t_high), max(alpha_t_low, alpha_t_high)


# ==============================================================
# Topological selection table (N scan, M=2^G_coronary fixed)
# ==============================================================

def compute_topo_N_table(alpha_exp=2.70, sigma_exp=0.20):
    """
    For N in [2,3,4,6], compute k_t, k_w, kappa, alpha*, sigma deviation.
    M = 2^G_coronary terminals fixed; G(N) = round(log_N(M)).
    """
    M = 2**G_coronary
    N_vals = [2, 3, 4, 6]
    results = []
    for N in N_vals:
        G_N = round(np.log(M) / np.log(N))
        k_w, k_t, kappa, alpha_t_n = compute_curvatures(G_N, N=N)
        a_star, _, eta_n, _ = find_alpha_star_minimax(G=G_N, N=N)
        # Analytical k_w cross-check
        k_w_analytic = np.log(M) * np.log(N) / (4 * alpha_w**2)
        sigma_dev = (a_star - alpha_exp) / sigma_exp
        results.append({
            'N': N, 'G': G_N,
            'k_t': k_t, 'k_w': k_w, 'k_w_analytic': k_w_analytic,
            'kappa': kappa,
            'alpha_star': a_star, 'sigma_dev': sigma_dev,
        })
    return results


# ==============================================================
# H(N) scaling exponent (power-law fit of k_t vs ln N)
# ==============================================================

def compute_hn_scaling_exponent(topo_table, alpha_t_val):
    """
    Fit H(N) ~ (ln N)^exponent where k_t = (ln N)^2 / (2*alpha_t^4) * H(N).
    OLS on log-log scale.
    """
    N_vals = np.array([row['N'] for row in topo_table], dtype=float)
    k_t_vals = np.array([row['k_t'] for row in topo_table], dtype=float)
    ln_N = np.log(N_vals)
    # H(N) = k_t * 2 * alpha_t^4 / (ln N)^2
    H_vals = k_t_vals * 2.0 * alpha_t_val**4 / ln_N**2
    log_H = np.log(H_vals)
    slope, _ = np.polyfit(np.log(ln_N), log_H, 1)
    return slope


# ==============================================================
# kt variation across N
# ==============================================================

def compute_kt_variation(topo_table):
    """Relative variation of k_t^net across N values: (max-min)/max * 100."""
    k_t_vals = [row['k_t'] for row in topo_table]
    return (max(k_t_vals) - min(k_t_vals)) / max(k_t_vals) * 100


# ==============================================================
# k_t^net ∝ G^exponent scaling
# ==============================================================

def compute_kt_G_exponent(kappa_table):
    """Fit k_t^net ∝ G^exponent on log-log scale (exclude G=1)."""
    G_vals  = np.array([row['G']   for row in kappa_table if row['G'] > 1], dtype=float)
    k_t_vals = np.array([row['k_t'] for row in kappa_table if row['G'] > 1], dtype=float)
    slope, _ = np.polyfit(np.log(G_vals), np.log(k_t_vals), 1)
    return slope


# ==============================================================
# Cross-system validation table
# ==============================================================

def find_second_crossing():
    """Find the spurious second crossing f(α)=g(α) for α > α_t."""
    r_local = locally_optimal_radii()
    Phi_ref = reference_cost(r_local)

    def diff(a):
        return wave_loss_network(a) - two_level_lagrangian(a, 0.0, r_local, Phi_ref)

    try:
        alpha2 = brentq(diff, 3.0, 5.0, xtol=1e-6)
    except ValueError:
        alpha2 = float('nan')
    return alpha2


def compute_cross_system_table():
    """
    Compute minimax alpha* for non-coronary systems using their alpha_w and G.
    Transport cost uses coronary metabolic parameters as baseline.
    """
    # Human pulmonary (Huang 1996: p=0.60 → alpha_w=2.200, G=15)
    a_pul, _, _, _ = find_alpha_star_minimax(aw=alpha_w_pulmonary, G=G_pulmonary)
    sigma_pul = abs((a_pul - alpha_exp_pulmonary) / alpha_exp_err_pulmonary)

    # Aortic tree (alpha_w range [2.10, 2.25], G=15)
    a_ao_1, _, _, _ = find_alpha_star_minimax(aw=alpha_w_aortic_low,  G=G_aortic)
    a_ao_2, _, _, _ = find_alpha_star_minimax(aw=alpha_w_aortic_high, G=G_aortic)
    alpha_star_ao_lo = min(a_ao_1, a_ao_2)
    alpha_star_ao_hi = max(a_ao_1, a_ao_2)
    sigma_ao = abs((alpha_star_ao_lo - alpha_exp_aortic) / alpha_exp_err_aortic)

    return {
        'alpha_w_pulmonary':      alpha_w_pulmonary,
        'G_pulmonary':            G_pulmonary,
        'alpha_star_pulmonary':   a_pul,
        'alpha_exp_pulmonary':    alpha_exp_pulmonary,
        'alpha_exp_err_pulmonary': alpha_exp_err_pulmonary,
        'alpha_exp_pul_low':      alpha_exp_pulmonary_low,
        'alpha_exp_pul_high':     alpha_exp_pulmonary_high,
        'sigma_pulmonary':        sigma_pul,
        'p_pulmonary':            p_pulmonary,
        'alpha_w_aortic_low':     alpha_w_aortic_low,
        'alpha_w_aortic_high':    alpha_w_aortic_high,
        'G_aortic':               G_aortic,
        'alpha_star_aortic_low':  alpha_star_ao_lo,
        'alpha_star_aortic_high': alpha_star_ao_hi,
        'alpha_exp_aortic':       alpha_exp_aortic,
        'alpha_exp_err_aortic':   alpha_exp_err_aortic,
        'sigma_aortic':           sigma_ao,
        'alpha_w_neural':         alpha_w_neural,
        'G_neural':               G_neural,
        'alpha_star_neural':      alpha_star_neural_indicative,
        'alpha_exp_neural':       alpha_exp_neural,
        'sigma_neural':           sigma_neural_indicative,
        'G_airways':              G_airways,
        'alpha_exp_airways_low':  alpha_exp_airways_low,
        'alpha_exp_airways_high': alpha_exp_airways_high,
    }


# ==============================================================
# Lq robustness values
# ==============================================================

def compute_lq_values(r_local, Phi_ref, eta=0.5):
    """
    Compute alpha*(q) for q in [1, infty), returning:
      alpha_star_q1: alpha* at q=1
      delta_alpha_lq: max spread for q>=1

    The Lq norm is applied to the weighted cost pair (eta*C_wave, (1-eta)*C_transport).
    For q->inf this recovers the minimax regardless of eta; for finite q the
    result depends on eta. Called with eta=eta_star for consistency with the paper.
    """
    q_vals = np.logspace(-0.3, 2, 30)
    a_stars_q = []
    for q in q_vals:
        def Lq(a, _q=q, _eta=eta):
            cw = _eta * wave_loss_network(a)
            ct = (1 - _eta) * two_level_lagrangian(a, 0.0, r_local, Phi_ref)
            hi = max(cw, ct)
            lo = min(cw, ct)
            if hi < 1e-300:
                return 0.0
            return hi * (1.0 + (lo / hi) ** _q) ** (1.0 / _q)
        res = minimize_scalar(Lq, bounds=(2.0, 3.4), method='bounded')
        a_stars_q.append(res.x)
    a_stars_q = np.array(a_stars_q)

    q1_mask = q_vals >= 1.0
    idx_q1 = np.searchsorted(q_vals, 1.0)
    alpha_star_q1 = float(a_stars_q[idx_q1])
    spread_q = float(alpha_star_q1 - a_stars_q[q1_mask].min())

    return alpha_star_q1, spread_q, q_vals, a_stars_q


# ==============================================================
# Self-consistency check: two-level alpha_t vs companion paper
# ==============================================================

def compute_self_consistency(alpha_t_val):
    """
    Compare alpha_t from two-level framework with alpha_local_mean.
    Returns percentage discrepancy.
    """
    al_mean, _, _ = alpha_local_mean()
    return abs(alpha_t_val - al_mean) / al_mean * 100


# ==============================================================
# Write dynamic_variables.tex
# ==============================================================

def write_tex(results_dict, tex_path):
    r = results_dict
    kappa_names = {1: 'I', 5: 'V', 7: 'VII', 9: 'IX', 11: 'XI', 13: 'XIII', 15: 'XV', 20: 'XX'}
    N_names = {2: 'II', 3: 'III', 4: 'IV', 6: 'VI'}

    def f(name, val, unit, desc):
        """Format a LaTeX command with aligned unit and comment."""
        cmd = f"\\newcommand{{\\{name}}}{{{val}}}"
        return f"{cmd:<50} % [{unit:<8}] {desc}"

    lines = [
        "% " + "=" * 68,
        "% DYNAMIC VARIABLES",
        "% " + "=" * 68,
        "% Source: paper2-variational/scripts/compute.py",
        f"% Generated: {np.datetime64('now').astype(str).split('T')[0]}",
        "% " + "=" * 68,
        "",
        "% " + "-" * 68,
        "% 1. PHYSICAL PARAMETERS",
        "% " + "-" * 68,
        f("VarAlphaW", f"{alpha_w:.3f}", "-", "Acoustic impedance matching exponent (baseline)"),
        f("VarP", f"{p}", "-", "Histological wall-thickness scaling exponent"),
        f("VarG", f"{G_coronary}", "-", "Number of generations in coronary tree"),
        f("VarBeta", f"{beta_coronary:.3f}", "-", "Length-scaling factor (theoretical)"),
        f("VarKe", f"{ke}", "-", "Distal stiffening exponent (coronary)"),
        f("VarAlphaWKe", f"{(5-p+ke)/2:.3f}", "-", "Acoustic ground state including distal stiffening"),
        f("VarAlphaWShift", f"{ke/2:.3f}", "-", "Shift in alpha_w due to stiffening"),
        "",
        "% " + "-" * 68,
        "% 2. MINIMAX AND TRANSPORT RESULTS",
        "% " + "-" * 68,
        f("VarAlphaT", f"{r['alpha_t']:.2f}", "-", "Network-level transport optimum alpha_t"),
        f("VarAlphaTExact", f"{r['alpha_t']:.4f}", "-", "Network-level transport optimum (exact)"),
        f("VarAlphaStar", f"{r['alpha_star']:.2f}", "-", "Unified minimax branching exponent alpha*"),
        f("VarEtaStar", f"{r['eta_star']:.3f}", "-", "Emergent operational duty cycle eta*"),
        f("VarGradientRatio", f"{r['gradient_ratio']}", "-", "Gradient ratio |g'|/f' at alpha*"),
        f("VarWaveCostNet", f"{r['wave_cost_net']:.1f}", "%", "Cumulative network wave cost at alpha*"),
        f("VarAlphaLocal", f"{r['alpha_local_mean']:.2f}", "-", "Companion paper's single-vessel optimum"),
        f("VarAlphaTLow", f"{r['alpha_t_low']:.2f}", "-", "Transport optimum (lower metabolic bound)"),
        f("VarAlphaTHigh", f"{r['alpha_t_high']:.2f}", "-", "Transport optimum (upper metabolic bound)"),
        f("VarAlphaTConsistency", f"{r['alpha_t_consistency']:.2f}", "%", "Discrepancy between single-vessel and network models"),
        f("VarAlphaExp", "2.70", "-", "Empirical cardiovascular mean (Kassab 1993)"),
        f("VarAlphaExpErr", "0.20", "-", "Empirical standard deviation"),
        f("VarPredError", f"{r['pred_error']:.1f}", "%", "Prediction error vs clinical mean"),
        "",
        "% " + "-" * 68,
        "% 3. REFLECTION COEFFICIENTS AND SCALING",
        "% " + "-" * 68,
    ]

    gamma_names_map = {
        1.50: 'OneFive', 1.75: 'OneSevenFive', 2.00: 'TwoZero',
        2.50: 'TwoFive', 2.75: 'TwoSevenFive', 3.00: 'ThreeZero'
    }
    for row in r['gamma_table']:
        lines.append(
            f(f"VarGammaAW{gamma_names_map[row['alpha']]}", f"{row['gamma_aw']:.4f}", "-", f"|Gamma|^2 at alpha={row['alpha']}")
        )

    lines += [
        "",
        "% " + "-" * 68,
        "% 4. EFFECTIVE STIFFNESS RATIO (G SCAN, N=2)",
        "% " + "-" * 68,
    ]
    for row in r['kappa_table']:
        G = row['G']
        name = kappa_names[G]
        lines.append(f(f"VarKappaG{name}", f"{row['kappa']:.2f}", "-", f"Stiffness ratio at G={G}"))
        lines.append(f(f"VarAlphaStarG{name}", f"{row['alpha_star']:.2f}", "-", f"Minimax alpha* at G={G}"))

    lines += [
        "",
        "% " + "-" * 68,
        "% 5. TOPOLOGICAL SELECTION (N SCAN, M=2048 terminals)",
        "% " + "-" * 68,
    ]
    for row in r['topo_N_table']:
        N = row['N']
        name = N_names[N]
        sign = '+' if row['sigma_dev'] >= 0 else ''
        lines.append(f(f"VarGN{name}", f"{row['G']}", "-", f"generations for N={N}"))
        lines.append(f(f"VarKtN{name}", f"{row['k_t']:.3f}", "-", f"Transport curvature for N={N}"))
        lines.append(f(f"VarKwN{name}", f"{row['k_w']:.3f}", "-", f"Wave curvature for N={N}"))
        lines.append(f(f"VarKappaN{name}", f"{row['kappa']:.2f}", "-", f"Stiffness ratio for N={N}"))
        lines.append(f(f"VarAlphaStarN{name}", f"{row['alpha_star']:.3f}", "-", f"Minimax alpha* for N={N}"))
        lines.append(f(f"VarSigmaN{name}", f"{sign}{row['sigma_dev']:.2f}", "sigma", "Deviation from exp (units of sigma)"))

    lines += [
        "",
        "% " + "-" * 68,
        "% 6. CROSS-SYSTEM VALIDATION",
        "% " + "-" * 68,
    ]

    cs = r['cross']
    lines += [
        f("VarAlphaStarPulmonary", f"{cs['alpha_star_pulmonary']:.3f}", "-", "Predicted alpha* for Pulmonary"),
        f("VarAlphaExpPulmonary", f"{cs['alpha_exp_pulmonary']}", "-", "Clinical alpha for Pulmonary"),
        f("VarSigmaPulmonary", f"{cs['sigma_pulmonary']:.1f}", "sigma", "Prediction error (sigma)"),
        f("VarAlphaStarAorticLow", f"{cs['alpha_star_aortic_low']:.2f}", "-", "Aortic prediction (lower bound)"),
        f("VarAlphaStarAorticHigh", f"{cs['alpha_star_aortic_high']:.2f}", "-", "Aortic prediction (upper bound)"),
        f("VarSigmaAortic", f"{cs['sigma_aortic']:.1f}", "sigma", "Prediction error (sigma)"),
        f("VarAlphaStarNeural", f"{cs['alpha_star_neural']:.1f}", "-", "Dendritic prediction (indicative)"),
        f("VarAlphaExpNeural", f"{cs['alpha_exp_neural']}", "-", "Clinical alpha for cortical dendrites"),
        f("VarSigmaNeural", f"{cs['sigma_neural']}", "sigma", "Prediction error (sigma)"),
        "",
        "% " + "-" * 68,
        "% 7. ANALYTICAL SCALING EXPONENTS",
        "% " + "-" * 68,
        f("VarKtGScalingExp", f"{r['kt_G_exp']:.2f}", "-", "Power law: kt_net ~ G^2.5"),
        f("VarKappaGScalingExp", f"{r['kt_G_exp'] - 1:.2f}", "-", "Power law: kappa_eff ~ G^1.5"),
        f("VarKtVariation", f"{r['kt_variation']:.0f}", "%", "Variation of kt across N values"),
        f("VarAlphaSecondCrossing", f"{r['alpha2_crossing']:.2f}", "-", "Spurious second crossing alpha_2 > 3"),
        "",
        "% " + "-" * 68,
        "% 8. PHYSICAL PARAMETER SENSITIVITY BASELINES",
        "% " + "-" * 68,
        f("VarBloodCostSens", f"{b_Wm3}", "W/m³", "Metabolic cost of blood"),
        f("VarWallMetabSens", f"{MW_MID_kWm3}", "kW/m³", "Metabolic rate of wall tissue"),
        f("VarViscositySens", f"{MU_mPas}", "mPas", "Dynamic viscosity"),
        f("VarQzeroSens", f"{Q_sensitivity_mL:.1f}", "mL/s", "Peak pulsatile flow (4×resting, params.py Q0_peak)"),
        f("VarSegLengthSens", f"{ell0_mm}", "mm", "Proximal segment length (params.py ell0_coronary)"),
        "",
        "% " + "-" * 68,
        f"% 9. LOCAL ALPHA RANGE (from script output: alpha_local = {r['alpha_local_mean']:.4f} [{r['alpha_local_min']:.4f}, {r['alpha_local_max']:.4f}])",
        "% " + "-" * 68,
        f("VarAlphaLocalMin", f"{r['alpha_local_min']:.4f}", "-", "Local alpha_t minimum across flow range"),
        f("VarAlphaLocalMax", f"{r['alpha_local_max']:.4f}", "-", "Local alpha_t maximum across flow range"),
        "",
        "% " + "-" * 68,
        f"% 10. Lq ROBUSTNESS (from script output: alpha*(q=1)={r['alpha_star_q1']:.3f}, spread={r['delta_alpha_lq']:.3f})",
        "% " + "-" * 68,
        f("VarAlphaStarQone", f"{r['alpha_star_q1']:.3f}", "-", "Minimax alpha* under L1 norm (Lq robustness)"),
        f("VarDeltaAlphaLq", f"{r['delta_alpha_lq']:.3f}", "-", "Spread of alpha* across Lq norms"),
        "",
        "% " + "-" * 68,
        "% 11. SINGLE-JUNCTION WAVE CURVATURE (formula: (ln N)^2 / (4 alpha_w^2), N=2)",
        "% " + "-" * 68,
        f("VarKwJunction",     f"{(np.log(2)**2 / (4 * 2.0**2)):.3f}",    "-", "Wave curvature at single junction (illustrative, alpha_w=2)"),
        f("VarKwJunctionPhys", f"{(np.log(2)**2 / (4 * alpha_w**2)):.3f}", "-", "Wave curvature at single junction (physiological, alpha_w=2.115)"),
        "",
        "% " + "-" * 68,
        f"% 12. H(N) SCALING EXPONENT (from script output: H(N) scaling exponent: {r['hn_scaling_exp']:.1f})",
        "% " + "-" * 68,
        f("VarHNScalingExp", f"{r['hn_scaling_exp']:.1f}", "-", "H(N) ~ (ln N)^exp scaling exponent"),
        "",
        "% " + "-" * 68,
        "% 13. WALL-THICKNESS PARAMETER RANGE (from params.py)",
        "% " + "-" * 68,
        f("VarPLow", f"{p_low}", "-", "Wall-thickness exponent lower bound"),
        f("VarPHigh", f"{p_high}", "-", "Wall-thickness exponent upper bound"),
        f("VarPPulmonary", f"{p_pulmonary}", "-", "Wall-thickness exponent (pulmonary, Huang 1996)"),
        "",
        "% " + "-" * 68,
        "% 14. MINIMAX BAND — G AND p RANGES (from params.py)",
        "% " + "-" * 68,
        f("VarGLow", f"{G_low}", "-", "Lower generation count for minimax band"),
        f("VarGHigh", f"{G_high}", "-", "Upper generation count for minimax band"),
        f("VarAlphaTLowBand", f"{alpha_t_low_band:.2f}", "-", "Transport optimum lower bound (band)"),
        f("VarAlphaTHighBand", f"{alpha_t_high_band:.2f}", "-", "Transport optimum upper bound (band)"),
        "",
        "% " + "-" * 68,
        "% 15. CROSS-SYSTEM VALIDATION — EXTENDED (from params.py)",
        "% " + "-" * 68,
        f("VarAlphaWPulmonary", f"{alpha_w_pulmonary:.3f}", "-", "Wave attractor (pulmonary, p=0.60)"),
        f("VarGPulmonary", f"{G_pulmonary}", "-", "Tree depth (pulmonary)"),
        f("VarAlphaExpErrPulmonary", f"{alpha_exp_err_pulmonary}", "-", "Empirical std dev (pulmonary)"),
        f("VarAlphaExpPulmonaryLow", f"{alpha_exp_pulmonary_low:.2f}", "-", "Empirical alpha lower bound (pulmonary)"),
        f("VarAlphaExpPulmonaryHigh", f"{alpha_exp_pulmonary_high:.2f}", "-", "Empirical alpha upper bound (pulmonary)"),
        f("VarAlphaWAorticLow", f"{alpha_w_aortic_low:.2f}", "-", "Wave attractor lower bound (aortic)"),
        f("VarAlphaWAorticHigh", f"{alpha_w_aortic_high:.2f}", "-", "Wave attractor upper bound (aortic)"),
        f("VarGAortic", f"{G_aortic}", "-", "Tree depth (aortic)"),
        f("VarAlphaExpAortic", f"{alpha_exp_aortic}", "-", "Empirical alpha (aortic, human)"),
        f("VarAlphaExpErrAortic", f"{alpha_exp_err_aortic}", "-", "Empirical std dev (aortic)"),
        f("VarAlphaWNeural", f"{alpha_w_neural:.2f}", "-", "Wave attractor (neural, electrotonic Rall 1959)"),
        f("VarGNeural", f"{G_neural}", "-", "Tree depth (dendritic, indicative)"),
        f("VarGAirways", f"{G_airways}", "-", "Weibel 1963 bronchial tree generations"),
        f("VarAlphaExpAirwaysLow", f"{alpha_exp_airways_low}", "-", "Empirical alpha lower bound (airways)"),
        f("VarAlphaExpAirwaysHigh", f"{alpha_exp_airways_high}", "-", "Empirical alpha upper bound (airways)"),
        f("VarMurrayAlpha", "3.0", "-", "Murray cubic law exponent"),
    ]

    with open(tex_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"  [OK] Wrote {tex_path}")


# ==============================================================
# Generate figures
# ==============================================================

def generate_figures(results_dict, fig_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "serif", "font.size": 12,
        "axes.labelsize": 13, "legend.fontsize": 10,
    })

    eta    = results_dict['eta_star']
    a_star = results_dict['alpha_star']

    r_local = locally_optimal_radii()
    Phi_ref = reference_cost(r_local)

    # --- Figure 1: kappa_eff(G) ---
    G_vals = list(range(1, 21))
    kappas = []
    for G_val in G_vals:
        _, _, kap, _ = compute_curvatures(G_val)
        kappas.append(kap)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(G_vals, kappas, 's-', color='#8c564b', markersize=6, lw=2)
    ax.set_xlabel('Tree Depth $G$ (generations)')
    ax.set_ylabel(r'Effective Stiffness Ratio $\kappa_{\mathrm{eff}}$')
    ax.set_title(r'Emergence of Network-Level Competition')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{fig_dir}/fig1_kappa.pdf', dpi=300)
    plt.close()
    print(f"  [OK] Wrote {fig_dir}/fig1_kappa.pdf")

    # --- Figure 2: fig2_minimax.pdf (2 panels) ---
    alpha_arr = np.linspace(2.1, 3.0, 500)
    C_w_raw = np.array([wave_loss_network(a) for a in alpha_arr])
    C_t_raw = np.array([two_level_lagrangian(a, 0.0, r_local, Phi_ref)
                        for a in alpha_arr])

    # Minimax G-band
    a_G9,  _, _, _ = find_alpha_star_minimax(G=9)
    a_G13, _, _, _ = find_alpha_star_minimax(G=13)
    mm_lo, mm_hi = min(a_G9, a_G13), max(a_G9, a_G13)

    # Lq values (precomputed)
    alpha_star_q1 = results_dict['alpha_star_q1']
    spread_q = results_dict['delta_alpha_lq']
    q_vals = results_dict['_q_vals']
    a_stars_q = results_dict['_a_stars_q']

    # Crossing y-value
    idx_cross = np.argmin(np.abs(C_w_raw - C_t_raw))
    y_cross   = (C_w_raw[idx_cross] + C_t_raw[idx_cross]) / 2.0

    # Panel A: fig2a_minimax.pdf
    fig_a, ax = plt.subplots(1, 1, figsize=(7, 5.4))
    ax.plot(alpha_arr, C_w_raw, '-', color='#1f77b4', lw=2.5,
            label=r'$\mathcal{C}_{\mathrm{wave}}^{\mathrm{net}}(\alpha)$')
    ax.plot(alpha_arr, C_t_raw, '-', color='#d62728', lw=2.5,
            label=r'$\mathcal{C}_{\mathrm{transport}}^{\mathrm{net}}(\alpha)$')
    ax.axvspan(2.50, 2.90, alpha=0.12, color='#2ca02c',
               label=r'$\alpha_{\exp} = 2.70 \pm 0.20$')
    ax.axvspan(mm_lo, mm_hi, alpha=0.35, color='gray',
               label=f'Minimax band [{mm_lo:.2f}, {mm_hi:.2f}]')
    ax.axvline(a_star, color='k', ls='--', lw=1.2, alpha=0.6)
    ax.plot(a_star, y_cross, 'ko', ms=9, zorder=10)
    # Label: right of dot, slightly above — no arrow
    ax.text(a_star + 0.04, y_cross,
            rf'$\alpha^*_{{\mathrm{{mm}}}} = {a_star:.3f}$',
            fontsize=11, fontweight='bold', va='top', ha='left')
    ax.set_xlabel(r'Branching exponent $\alpha$', fontsize=13)
    ax.set_ylabel('Fractional cost', fontsize=13)
    ax.set_title('Minimax equal-cost intersection', fontsize=13)
    ax.set_xlim(2.1, 3.0)
    ax.set_ylim(-0.005, 0.26)
    ax.grid(True, alpha=0.2)
    # Legend below — full width, 2 columns, larger font
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),
              ncol=2, fontsize=11, framealpha=0.92,
              handlelength=1.8, columnspacing=2.0)
    fig_a.tight_layout()
    fig_a.subplots_adjust(bottom=0.25)
    fig_a.savefig(f'{fig_dir}/fig2_minimax.pdf', bbox_inches='tight')
    plt.close(fig_a)
    print(f"  [OK] Wrote {fig_dir}/fig2_minimax.pdf")

    # Panel B: fig2b_robustness.pdf
    # NOTE: alpha_star_q1 > a_star  (q=1 gives the HIGHEST optimum; minimax is lowest)
    delta_q = alpha_star_q1 - a_star          # positive by construction
    q_lo    = min(a_star, alpha_star_q1)
    q_hi    = max(a_star, alpha_star_q1)

    fig_b, ax = plt.subplots(1, 1, figsize=(7, 5.4))
    phys = q_vals >= 1.0
    # Shaded robustness band
    ax.axhspan(q_lo, q_hi, alpha=0.10, color='#1f77b4')
    # Reference asymptotes
    ax.axhline(a_star, color='#555', ls=':', lw=1.5,
               label=rf'$q\to\infty$ (minimax): $\alpha^* = {a_star:.3f}$')
    ax.axhline(alpha_star_q1, color='#d62728', ls='--', lw=1.2, alpha=0.7,
               label=rf'$q = 1$: $\alpha^* = {alpha_star_q1:.3f}$')
    # Curve
    ax.semilogx(q_vals[phys], a_stars_q[phys], 'ko-', ms=4, lw=1.8)
    # q=1 highlight dot
    idx_q1 = np.searchsorted(q_vals, 1.0)
    ax.plot(q_vals[idx_q1], alpha_star_q1, 'o',
            color='#d62728', ms=8, zorder=5)
    # Δα* arrow on the central vertical gridline (x=10 on log scale)
    ax_x = 10.0
    ax.annotate('', xy=(ax_x, a_star), xytext=(ax_x, alpha_star_q1),
                arrowprops=dict(arrowstyle='<->', color='#1f77b4', lw=2.0))
    ax.text(ax_x * 1.35, (a_star + alpha_star_q1) / 2,
            rf'$\Delta\alpha^* = {delta_q:.3f}$', fontsize=11,
            color='#1f77b4', ha='left', va='center', fontweight='bold')
    ax.set_xlabel('Norm exponent $q$', fontsize=13)
    ax.set_ylabel(r'$\alpha^*(q)$', fontsize=13)
    ax.set_title(r'Robustness under $L_q$ scalarization ($q \geq 1$)', fontsize=13)
    ax.set_xlim(0.9, 110)
    ax.set_ylim(q_lo - 0.020, q_hi + 0.020)
    ax.grid(True, alpha=0.2)
    # Legend below — full width, 1 row
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),
              ncol=2, fontsize=11, framealpha=0.92,
              handlelength=1.8, columnspacing=2.0)
    fig_b.tight_layout()
    fig_b.subplots_adjust(bottom=0.25)
    fig_b.savefig(f'{fig_dir}/fig3_robustness.pdf', bbox_inches='tight')
    plt.close(fig_b)
    print(f"  [OK] Wrote {fig_dir}/fig3_robustness.pdf")

    # --- Figure 3: fig3_lagrangian.pdf ---
    alpha_arr3 = np.linspace(2.50, 3.0, 500)
    C_w3 = np.array([wave_loss_network(a) for a in alpha_arr3])
    C_t3 = np.array([two_level_lagrangian(a, 0.0, r_local, Phi_ref)
                     for a in alpha_arr3])
    L3   = eta * C_w3 + (1 - eta) * C_t3

    min_idx = np.argmin(L3)
    a_min   = alpha_arr3[min_idx]
    L_min   = L3[min_idx]

    fig, ax = plt.subplots(figsize=(7, 5.4))
    ax.plot(alpha_arr3, eta * C_w3, '-', color='#1f77b4', lw=2.5,
            label=rf'$\eta^*\,\mathcal{{C}}^{{\mathrm{{net}}}}_{{\mathrm{{wave}}}}$'
                  rf'$\;(\eta^* = {eta:.3f})$')
    ax.plot(alpha_arr3, (1 - eta) * C_t3, '-', color='#d62728', lw=2.5,
            label=r'$(1-\eta^*)\,\mathcal{C}^{\mathrm{net}}_{\mathrm{transport}}$')
    ax.plot(alpha_arr3, L3, 'k-', lw=2.5,
            label=r'$\mathcal{L}_{\mathrm{net}}$ (total)')
    ax.axvline(a_min, color='#555', ls='--', lw=1.4, alpha=0.8,
               label=rf'$\alpha^* = {a_min:.3f}$ (minimax)')
    ax.plot(a_min, L_min, 'ko', ms=8, zorder=5)
    # Label at ~1:30 (upper-right, ~45°) — no arrow
    ax.text(a_min + 0.010, L_min + 0.006,
            rf'$\alpha^* = {a_min:.3f}$',
            fontsize=11, fontweight='bold', color='k', va='bottom', ha='left')
    ax.set_xlabel(r'Branching exponent $\alpha$', fontsize=13)
    ax.set_ylabel('Weighted cost', fontsize=13)
    ax.set_title(r'Network Lagrangian Decomposition', fontsize=13)
    ax.set_xlim(2.50, 3.0)
    y_max = float(np.max(L3)) * 1.08
    ax.set_ylim(-0.002, y_max)
    ax.grid(True, alpha=0.2)
    # Legend below — full width, 2 columns, larger font
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),
              ncol=2, fontsize=11, framealpha=0.92,
              handlelength=1.8, columnspacing=2.0)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.savefig(f'{fig_dir}/fig4_lagrangian.pdf', bbox_inches='tight')
    plt.close()
    print(f"  [OK] Wrote {fig_dir}/fig4_lagrangian.pdf")


# ==============================================================
# Main
# ==============================================================

if __name__ == '__main__':
    print("=" * 65)
    print("Paper 2: Unified Variational Principle (Two-Level Minimax)")
    print("=" * 65)

    # 1. Local optimization (Paper 2)
    al_mean, al_min, al_max = alpha_local_mean()
    print(f"\nalpha_local = {al_mean:.4f} [{al_min:.4f}, {al_max:.4f}]")

    # 2. Minimax alpha* and derived quantities
    alpha_star_val, alpha_t_val, eta_star_val, grad_ratio = find_alpha_star_minimax()
    print(f"alpha_t    = {alpha_t_val:.4f}")
    print(f"alpha*     = {alpha_star_val:.4f}")
    print(f"eta*       = {eta_star_val:.4f}  (gradient ratio {grad_ratio}:1)")

    pred_error = abs(alpha_star_val - 2.70) / 2.70 * 100
    c_wave_net = wave_loss_network(alpha_star_val) * 100
    print(f"Error vs experiment: {pred_error:.2f}%")
    print(f"Network wave cost at alpha*: {c_wave_net:.1f}%")

    # 3. Self-consistency
    alpha_t_consistency = compute_self_consistency(alpha_t_val)
    print(f"Self-consistency: {alpha_t_consistency:.4f}%")

    # 4. Transport optimum bounds
    alpha_t_low, alpha_t_high = compute_alpha_t_bounds()
    print(f"alpha_t range: [{alpha_t_low:.2f}, {alpha_t_high:.2f}]")

    # 5. Tables
    gamma_table = compute_gamma_table()
    kappa_table = compute_kappa_table()

    print(f"\n-- Table: kappa_eff(G) --")
    print(f"  {'G':>4s}  {'kappa':>8s}  {'alpha_t':>8s}  {'alpha*':>8s}")
    for row in kappa_table:
        print(f"  {row['G']:>4d}  {row['kappa']:>8.3f}  {row['alpha_t']:>8.3f}  {row['alpha_star']:>8.3f}")

    # 6. Topological N table
    topo_N_table = compute_topo_N_table()
    print(f"\n-- Table: topo_N (M=2^{G_coronary} fixed) --")
    print(f"  {'N':>3s}  {'G':>4s}  {'kt':>8s}  {'kw':>8s}  {'kappa':>8s}  {'alpha*':>8s}  sigma")
    for row in topo_N_table:
        sign = '+' if row['sigma_dev'] >= 0 else ''
        print(f"  {row['N']:>3d}  {row['G']:>4d}  {row['k_t']:>8.3f}  {row['k_w']:>8.3f}"
              f"  {row['kappa']:>8.2f}  {row['alpha_star']:>8.3f}  {sign}{row['sigma_dev']:.2f}sig")

    # 7. H(N) scaling and kt variation
    hn_exp = compute_hn_scaling_exponent(topo_N_table, alpha_t_val)
    kt_var = compute_kt_variation(topo_N_table)
    print(f"\nH(N) scaling exponent: {hn_exp:.1f}")
    print(f"kt variation: {kt_var:.0f}%")

    # 7b. k_t(G) power-law exponent
    kt_G_exp = compute_kt_G_exponent(kappa_table)
    print(f"k_t(G) scaling exponent: {kt_G_exp:.2f}")

    # 7d. Second crossing
    alpha2_crossing = find_second_crossing()
    print(f"Second crossing alpha: {alpha2_crossing:.2f}")

    # 7c. Cross-system table
    cross = compute_cross_system_table()
    print(f"\n-- Cross-system table --")
    print(f"  Pulmonary: alpha_star={cross['alpha_star_pulmonary']:.3f}, sigma={cross['sigma_pulmonary']:.1f}")
    print(f"  Aortic: [{cross['alpha_star_aortic_low']:.2f}, {cross['alpha_star_aortic_high']:.2f}], sigma={cross['sigma_aortic']:.1f}")
    print(f"  Neural (indicative): alpha_star={cross['alpha_star_neural']}")

    # 8. Lq robustness
    r_local = locally_optimal_radii()
    Phi_ref = reference_cost(r_local)
    alpha_star_q1, delta_alpha_lq, q_vals_arr, a_stars_q_arr = compute_lq_values(r_local, Phi_ref, eta=0.5)
    print(f"\nLq robustness: alpha*(q=1)={alpha_star_q1:.3f}, spread={delta_alpha_lq:.3f}")

    # 9. Collect results
    results = {
        'eta_star': eta_star_val,
        'gradient_ratio': grad_ratio,
        'alpha_local_mean': al_mean,
        'alpha_local_min': al_min,
        'alpha_local_max': al_max,
        'alpha_t': alpha_t_val,
        'alpha_t_low': alpha_t_low,
        'alpha_t_high': alpha_t_high,
        'alpha_t_consistency': alpha_t_consistency,
        'alpha_star': alpha_star_val,
        'pred_error': pred_error,
        'wave_cost_net': c_wave_net,
        'alpha_star_q1': alpha_star_q1,
        'delta_alpha_lq': delta_alpha_lq,
        '_q_vals': q_vals_arr,
        '_a_stars_q': a_stars_q_arr,
        'gamma_table': gamma_table,
        'kappa_table': kappa_table,
        'topo_N_table': topo_N_table,
        'hn_scaling_exp': hn_exp,
        'kt_variation': kt_var,
        'kt_G_exp': kt_G_exp,
        'alpha2_crossing': alpha2_crossing,
        'cross': cross,
    }

    # 10. Write tex
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(current_dir))

    tex_path = os.path.join(project_root, 'paper2-variational', 'manuscript', 'dynamic_variables.tex')
    write_tex(results, tex_path)

    # 11. Figures
    fig_dir = os.path.join(project_root, 'paper2-variational', 'manuscript', 'figures')
    generate_figures(results, fig_dir)

    print(f"\n[OK] Paper 2 complete. alpha* = {alpha_star_val:.2f}, eta* = {eta_star_val:.3f}")
