"""
compute_paper1.py -- Compute all numerical results for Paper 1.

Uses the TWO-LEVEL MODEL:
  Level 1 (Paper 2): Each vessel locally optimizes r*(Q) including wall cost.
  Level 2 (Paper 1): Network Lagrangian balances transport deviation
                     from local optima against wave matching cost.

Generates: dynamic_variables.tex for Paper 1
           fig1_kappa.pdf, fig2_lagrangian.pdf
"""

import numpy as np
from scipy.optimize import minimize_scalar, brentq
from params import *


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
    """Mean local branching exponent across the tree."""
    alphas = []
    for g in range(G):
        Q_p = Q0_coronary / 2**g
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
# Level 2: Network Lagrangian (two-level)
# ==============================================================

def gamma_squared(alpha, aw=alpha_w):
    """Reflected power fraction |Gamma|^2 at a single junction."""
    x = 2**(aw / alpha - 1)
    return ((x - 1) / (x + 1))**2


def wave_loss_network(alpha, G=G_coronary, aw=alpha_w):
    """Network wave cost: 1 - (1 - |Gamma|^2)^G."""
    g2 = gamma_squared(alpha, aw)
    return 1 - (1 - g2)**G


def locally_optimal_radii(mw=MW_MID, G=G_coronary):
    """Compute r*(Q_g) for g = 0..G."""
    return [r_star(Q0_coronary / 2**g, mw) for g in range(G + 1)]


def reference_cost(r_local, mw=MW_MID, G=G_coronary, beta=1.0):
    """Total cost when every vessel is at its local optimum."""
    return sum(
        2**g * phi_cost(r_local[g], Q0_coronary / 2**g, mw) * ell0_coronary * beta**g
        for g in range(G)
    )


def two_level_lagrangian(alpha, eta, r_local, Phi_ref,
                         mw=MW_MID, G=G_coronary, beta=1.0, aw=alpha_w):
    """
    Two-level network Lagrangian.
    
    C_transport: penalty for deviating from locally optimal radii
                 when imposing a uniform branching exponent alpha.
    C_wave:      cumulative wave reflection loss through G junctions.
    """
    # Imposed radii from uniform alpha, starting from r*(Q_0)
    r = r_local[0]
    Phi_imposed = 0.0
    for g in range(G):
        Phi_imposed += 2**g * phi_cost(r, Q0_coronary / 2**g, mw) * ell0_coronary * beta**g
        r *= 2**(-1.0 / alpha)
    
    C_transport = (Phi_imposed - Phi_ref) / Phi_ref
    C_wave = wave_loss_network(alpha, G, aw)
    
    return eta * C_wave + (1 - eta) * C_transport


def find_alpha_star(eta, mw=MW_MID, G=G_coronary, beta=1.0, aw=alpha_w):
    """Find alpha* from the two-level Lagrangian."""
    r_local = locally_optimal_radii(mw, G)
    Phi_ref = reference_cost(r_local, mw, G, beta)
    
    def L(a):
        return two_level_lagrangian(a, eta, r_local, Phi_ref, mw, G, beta, aw)
    
    res = minimize_scalar(L, bounds=(1.5, 4.0), method='bounded')
    
    # Also find transport-only minimum
    def L_transport_only(a):
        return two_level_lagrangian(a, 0.0, r_local, Phi_ref, mw, G, beta, aw)
    res_t = minimize_scalar(L_transport_only, bounds=(1.5, 4.0), method='bounded')
    
    return res.x, res_t.x


# ==============================================================
# Parametric duty cycle (Phenomenological)
# ==============================================================

def get_baseline_eta():
    """Baseline duty cycle, corresponding to the minimax η* ≈ 0.833."""
    return 0.833


# ==============================================================
# Curvatures and kappa_eff
# ==============================================================

def compute_curvatures(G, eta=None, mw=MW_MID, beta=1.0, aw=alpha_w):
    """Compute wave and transport curvatures, kappa_eff."""
    da = 0.001
    
    # Wave curvature at alpha_w
    k_w = 0.5 * (wave_loss_network(aw + da, G, aw)
                 - 2 * wave_loss_network(aw, G, aw)
                 + wave_loss_network(aw - da, G, aw)) / da**2
    
    # Transport minimum and curvature (two-level)
    r_local = locally_optimal_radii(mw, G)
    Phi_ref = reference_cost(r_local, mw, G, beta)
    
    def c_t(a):
        return two_level_lagrangian(a, 0.0, r_local, Phi_ref, mw, G, beta, aw)
    
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
    """Compute |Gamma|^2 for both alpha_w = 2 and alpha_w = 2.115."""
    alpha_vals = [1.50, 1.75, 2.00, 2.50, 2.75, 3.00]
    return [{'alpha': a, 'gamma_aw2': gamma_squared(a, 2.0),
             'gamma_aw2115': gamma_squared(a, alpha_w)} for a in alpha_vals]


# ==============================================================
# Duty cycle table
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
# kappa_eff table
# ==============================================================

def compute_kappa_table(eta, mw=MW_MID, beta=1.0):
    """Compute kappa_eff and alpha* as function of G."""
    G_vals = [1, 5, 7, 9, 11, 15, 20]
    results = []
    for G_val in G_vals:
        k_w, k_t, kappa, alpha_t = compute_curvatures(G_val, mw=mw, beta=beta)
        a_star, _ = find_alpha_star(eta, mw, G_val, beta)
        a_quad = alpha_harmonic(eta, alpha_w, alpha_t, kappa)
        results.append({
            'G': G_val, 'kappa': kappa, 'alpha_t': alpha_t,
            'alpha_star': a_star, 'alpha_quad': a_quad,
        })
    return results


# ==============================================================
# Write dynamic_variables.tex
# ==============================================================

def write_tex(results_dict, tex_path):
    r = results_dict
    lines = [
        "% DYNAMICALLY GENERATED -- DO NOT EDIT MANUALLY",
        "% Source: branching_scripts/compute_paper1.py",
        "% Model: Two-level (Paper 2 local + Paper 1 network wave)",
        f"% Generated: {np.datetime64('now')}",
        "",
        "% -- Physical parameters --",
        f"\\newcommand{{\\VarAlphaW}}{{{alpha_w:.3f}}}",
        f"\\newcommand{{\\VarP}}{{{p}}}",
        f"\\newcommand{{\\VarG}}{{{G_coronary}}}",
        f"\\newcommand{{\\VarBeta}}{{{beta_coronary}}}",
        "",
        "% -- Duty cycle --",
        f"\\newcommand{{\\VarEtaStar}}{{{r['eta']:.2f}}}",
    ]
    
    lines += [
        "",
        "% -- Two-level model results --",
        f"\\newcommand{{\\VarAlphaLocal}}{{{r['alpha_local_mean']:.2f}}}",
        f"\\newcommand{{\\VarAlphaT}}{{{r['alpha_t']:.2f}}}",
        f"\\newcommand{{\\VarAlphaStar}}{{{r['alpha_star']:.2f}}}",
        f"\\newcommand{{\\VarAlphaExp}}{{2.70}}",
        f"\\newcommand{{\\VarAlphaExpErr}}{{0.20}}",
        f"\\newcommand{{\\VarPredError}}{{{r['pred_error']:.1f}}}",
        f"\\newcommand{{\\VarWaveCostNet}}{{{r['wave_cost_net']:.1f}}}",
    ]
    
    # |Gamma|^2 table
    lines.append("")
    lines.append("% -- Reflection coefficients --")
    gamma_names = {1.50: 'OneFive', 1.75: 'OneSevenFive', 2.00: 'TwoZero', 2.50: 'TwoFive', 2.75: 'TwoSevenFive', 3.00: 'ThreeZero'}
    for row in r['gamma_table']:
        lines.append(f"\\newcommand{{\\VarGammaAW{gamma_names[row['alpha']]}}}{{{row['gamma_aw2115']:.4f}}}")
    
    # Kappa table
    lines.append("")
    lines.append("% -- Effective stiffness ratio --")
    kappa_names = {1: 'I', 5: 'V', 7: 'VII', 9: 'IX', 11: 'XI', 15: 'XV', 20: 'XX'}
    for row in r['kappa_table']:
        G = row['G']
        name = kappa_names[G]
        lines.append(f"\\newcommand{{\\VarKappaG{name}}}{{{row['kappa']:.2f}}}")
        lines.append(f"\\newcommand{{\\VarAlphaStarG{name}}}{{{row['alpha_star']:.2f}}}")
        lines.append(f"\\newcommand{{\\VarAlphaQuadG{name}}}{{{row['alpha_quad']:.2f}}}")
    
    with open(tex_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"  [OK] Wrote {tex_path}")


# ==============================================================
# Generate figures
# ==============================================================

def generate_figures(results_dict, fig_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.optimize import minimize_scalar

    plt.rcParams.update({
        "font.family": "serif", "font.size": 12,
        "axes.labelsize": 13, "legend.fontsize": 10,
    })

    eta    = results_dict['eta']
    a_star = results_dict['alpha_star']
    alpha_t_val = results_dict['alpha_t']

    r_local = locally_optimal_radii()
    Phi_ref = reference_cost(r_local)

    # --- Figure 1: kappa_eff(G) --- (unchanged) ---
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

    # --- Figure 2: fig2_minimax.png (2 panels) ---
    alpha_arr = np.linspace(2.1, 3.0, 500)
    C_w_raw = np.array([wave_loss_network(a) for a in alpha_arr])
    C_t_raw = np.array([two_level_lagrangian(a, 0.0, r_local, Phi_ref)
                        for a in alpha_arr])

    # Minimax G-band: alpha* for G in [9, 13]
    a_G9,  _ = find_alpha_star(eta, G=9)
    a_G13, _ = find_alpha_star(eta, G=13)
    mm_lo, mm_hi = min(a_G9, a_G13), max(a_G9, a_G13)

    # Panel B: alpha*(q) under L_q scalarization (fix eta, vary q)
    # L_q(alpha) = (C_w^q + C_t^q)^(1/q): equal-weight norm.
    # q=1  → min of sum   C_w+C_t  (minimum to the right of crossing)
    # q→∞  → min of max(C_w, C_t) = minimax equal-cost crossing = a_star
    # Numerically stable via factoring out the larger term.
    q_vals = np.logspace(-0.3, 2, 30)   # 0.5 … 100
    a_stars_q = []
    for q in q_vals:
        def Lq(a, _q=q):
            cw = wave_loss_network(a)
            ct = two_level_lagrangian(a, 0.0, r_local, Phi_ref)
            hi = max(cw, ct)
            lo = min(cw, ct)
            if hi < 1e-300:
                return 0.0
            return hi * (1.0 + (lo / hi) ** _q) ** (1.0 / _q)
        res = minimize_scalar(Lq, bounds=(2.0, 3.4), method='bounded')
        a_stars_q.append(res.x)
    a_stars_q = np.array(a_stars_q)

    q1_mask = q_vals >= 1.0
    a_star_q1 = a_stars_q[np.searchsorted(q_vals, 1.0)]
    spread_q  = float(a_stars_q[q1_mask].max() - a_stars_q[q1_mask].min())

    # Crossing y-value
    idx_cross = np.argmin(np.abs(C_w_raw - C_t_raw))
    y_cross   = (C_w_raw[idx_cross] + C_t_raw[idx_cross]) / 2.0

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A
    ax = axes[0]
    ax.plot(alpha_arr, C_w_raw, '-', color='#1f77b4', lw=2.5,
            label=r'$\mathcal{C}_{\mathrm{wave}}^{\mathrm{net}}(\alpha)$')
    ax.plot(alpha_arr, C_t_raw, '-', color='#d62728', lw=2.5,
            label=r'$\mathcal{C}_{\mathrm{transport}}^{\mathrm{net}}(\alpha)$')
    ax.axvspan(2.50, 2.90, alpha=0.15, color='#2ca02c',
               label=r'$\alpha_{\exp} = 2.70 \pm 0.20$')
    ax.axvspan(mm_lo, mm_hi, alpha=0.40, color='gray',
               label=f'Minimax band [{mm_lo:.2f}, {mm_hi:.2f}]')
    ax.axvline(a_star, color='k', ls='--', lw=1, alpha=0.7)
    ax.plot(a_star, y_cross, 'ko', ms=8)
    ax.annotate(rf'$\alpha^*_{{\mathrm{{mm}}}} = {a_star:.3f}$',
                xy=(a_star, y_cross),
                xytext=(a_star - 0.30, y_cross + 0.05),
                fontsize=12, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='k'))
    ax.set_xlabel(r'Branching exponent $\alpha$')
    ax.set_ylabel('Fractional cost')
    ax.set_title('(A) Minimax equal-cost intersection')
    ax.set_xlim(2.1, 3.0)
    ax.set_ylim(-0.005, 0.26)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.2)

    # Panel B
    ax = axes[1]
    # Only plot q >= 1 (physical range)
    phys = q_vals >= 1.0
    ax.semilogx(q_vals[phys], a_stars_q[phys], 'ko-', ms=5, lw=1.8)
    ax.axhline(a_star, color='#888', ls=':', lw=1.5,
               label=rf'$q\to\infty$ (minimax): $\alpha^* = {a_star:.3f}$')
    # Bracket showing the spread
    ax.annotate('', xy=(80, a_star), xytext=(80, a_star_q1),
                arrowprops=dict(arrowstyle='<->', color='#1f77b4', lw=1.5))
    ax.text(95, (a_star + a_star_q1) / 2,
            rf'$\Delta\alpha^* = {spread_q:.3f}$', fontsize=9,
            color='#1f77b4', ha='left', va='center')
    # Mark q=1
    ax.plot(q_vals[np.searchsorted(q_vals, 1.0)], a_star_q1, 'o',
            color='#d62728', ms=7, zorder=5)
    ax.text(1.15, a_star_q1 + 0.003,
            rf'$q=1$: $\alpha^* = {a_star_q1:.3f}$', fontsize=9, color='#d62728')
    ax.set_xlabel('Norm exponent $q$')
    ax.set_ylabel(r'$\alpha^*(q)$')
    ax.set_title(r'(B) Robustness under $L_q$ scalarization ($q \geq 1$)')
    ax.set_xlim(0.9, 110)
    ax.set_ylim(a_star - 0.005, a_star_q1 + 0.015)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/fig2_minimax.pdf', bbox_inches='tight')
    plt.close()
    print(f"  [OK] Wrote {fig_dir}/fig2_minimax.pdf")

    # --- Figure 3: fig3_lagrangian.pdf (Lagrangian component view) ---
    alpha_arr3 = np.linspace(2.50, 3.0, 500)
    C_w3 = np.array([wave_loss_network(a) for a in alpha_arr3])
    C_t3 = np.array([two_level_lagrangian(a, 0.0, r_local, Phi_ref)
                     for a in alpha_arr3])
    L3   = eta * C_w3 + (1 - eta) * C_t3

    min_idx = np.argmin(L3)
    a_min   = alpha_arr3[min_idx]
    L_min   = L3[min_idx]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(alpha_arr3, eta * C_w3, '-', color='#1f77b4', lw=2.5,
            label=rf'$\eta^*\,\mathcal{{C}}^{{\mathrm{{net}}}}_{{\mathrm{{wave}}}}$'
                  rf'$\;(\eta^* = {eta:.2f})$')
    ax.plot(alpha_arr3, (1 - eta) * C_t3, '-', color='#d62728', lw=2.5,
            label=r'$(1-\eta^*)\,\mathcal{C}^{\mathrm{net}}_{\mathrm{transport}}$')
    ax.plot(alpha_arr3, L3, 'k-', lw=2.5,
            label=r'$\mathcal{L}_{\mathrm{net}}$ (total)')

    # Mark minimax minimum with vertical dashed line + dot
    ax.axvline(a_min, color='#555', ls='--', lw=1.4, alpha=0.8,
               label=rf'$\alpha^* = {a_min:.3f}$ (minimax)')
    ax.plot(a_min, L_min, 'ko', ms=8, zorder=5)
    ax.text(a_min + 0.012, L_min + 0.004, rf'$\alpha^* = {a_star:.3f}$',
            fontsize=11, fontweight='bold', color='k')

    ax.set_xlabel(r'Branching exponent $\alpha$')
    ax.set_ylabel('Weighted cost')
    ax.set_title(r'Network Lagrangian Decomposition')
    ax.set_xlim(2.50, 3.0)
    # Use actual data range — no artificial cap that causes clipping artefacts
    y_max = float(np.max(L3)) * 1.08
    ax.set_ylim(-0.002, y_max)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(f'{fig_dir}/fig3_lagrangian.pdf', bbox_inches='tight')
    plt.close()
    print(f"  [OK] Wrote {fig_dir}/fig3_lagrangian.pdf")


# ==============================================================
# Main
# ==============================================================

if __name__ == '__main__':
    print("=" * 65)
    print("Paper 1: Unified Variational Principle (Two-Level Model)")
    print("=" * 65)
    
    # 1. Duty cycle
    eta_baseline = get_baseline_eta()
    print(f"\neta baseline = {eta_baseline:.2f}")
    
    # 2. Local optimization (Paper 2)
    al_mean, al_min, al_max = alpha_local_mean()
    print(f"alpha_local = {al_mean:.4f} [{al_min:.4f}, {al_max:.4f}]")
    
    # 3. Two-level Lagrangian
    alpha_star_val, alpha_t_val = find_alpha_star(eta_baseline)
    print(f"alpha_t(two-level) = {alpha_t_val:.4f}")
    print(f"alpha* = {alpha_star_val:.4f}")
    pred_error = abs(alpha_star_val - 2.7) / 2.7 * 100
    print(f"Error vs experiment: {pred_error:.2f}%")
    
    # Wave cost at alpha*
    c_wave_net = wave_loss_network(alpha_star_val) * 100 # In percentage
    print(f"Network wave reflection cost at alpha*: {c_wave_net:.1f}%")

    # 4. Q-dependence of alpha_t
    print("\nLocal alpha_t variation across generations:")
    for g in [0, 5, 10]:
        Q_g = Q0_coronary / 2**g
        r_p = r_star(Q_g)
        r_c = r_star(Q_g/2)
        al_g = np.log(2) / np.log(r_p / r_c)
        print(f"  Gen {g}: Q={Q_g:.2e} m^3/s -> alpha_t = {al_g:.4f}")
    
    # 5. Tables
    gamma_table = compute_gamma_table()
    eta_table = compute_eta_parametric_table()
    kappa_table = compute_kappa_table(eta_baseline)
    
    print(f"\n-- Table: kappa_eff(G) --")
    print(f"  {'G':>4s}  {'kappa':>8s}  {'alpha_t':>8s}  {'alpha*':>8s}  {'alpha_q':>8s}")
    for row in kappa_table:
        print(f"  {row['G']:>4d}  {row['kappa']:>8.3f}  {row['alpha_t']:>8.3f}  {row['alpha_star']:>8.3f}  {row['alpha_quad']:>8.3f}")
    
    # 6. Collect results
    results = {
        'eta': eta_baseline,
        'alpha_local_mean': al_mean,
        'alpha_t': alpha_t_val,
        'alpha_star': alpha_star_val,
        'pred_error': pred_error,
        'wave_cost_net': c_wave_net,
        'gamma_table': gamma_table,
        'eta_table': eta_table,
        'kappa_table': kappa_table,
    }
    
    # 7. Write tex
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(current_dir))
    
    tex_path = os.path.join(project_root, 'paper2-variational', 'manuscript', 'dynamic_variables.tex')
    write_tex(results, tex_path)
    
    # 8. Add extra info for paper 2
    a_star_60, _ = find_alpha_star(0.60)
    a_star_90, _ = find_alpha_star(0.90)
    print(f"\nExtra info for eta sensitivity (G = 11):")
    print(f"alpha* at eta=0.60 is: {a_star_60:.4f}")
    print(f"alpha* at eta=0.90 is: {a_star_90:.4f}")
    
    # 9. Figures
    fig_dir = os.path.join(project_root, 'paper2-variational', 'manuscript', 'figures')
    generate_figures(results, fig_dir)
    
    print(f"\n[OK] Paper 1 complete. alpha* = {alpha_star_val:.2f}")
