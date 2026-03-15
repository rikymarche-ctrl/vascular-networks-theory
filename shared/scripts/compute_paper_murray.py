"""
compute_paper_murray.py — Compute all numerical results for the Paper: Beyond Murray's Law.

Generates: dynamic_variables.tex for the Paper
           fig_alpha_scale.pdf

All outputs are deterministic and reproducible.
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from params import *


# ═══════════════════════════════════════════════════════════════
# Core functions
# ═══════════════════════════════════════════════════════════════

def dphi_dr(r, Q, mw=MW_MID):
    """∂Φ/∂r for the three-term cost function."""
    Bw = B_wall(mw)
    return -4 * A_of_Q(Q) * r**(-5) + 2 * B_blood * r + (1 + p) * Bw * r**p


def r_star(Q, mw=MW_MID):
    """Optimal radius r*(Q) from Theorem 1."""
    return brentq(dphi_dr, 1e-7, 0.5, args=(Q, mw), xtol=1e-14)


def alpha_star(Q, mw=MW_MID):
    """Local branching exponent alpha*(Q) at a symmetric bifurcation."""
    r0 = r_star(Q, mw)
    r1 = r_star(Q / 2, mw)
    return np.log(2) / np.log(r0 / r1)


def phi_at_optimum(r, mw=MW_MID):
    """Φ*(r) = cost per unit length at the optimal radius (Lemma 9)."""
    Bw = B_wall(mw)
    return 1.5 * B_blood * r**2 + (5 + p) / 4 * Bw * r**(1 + p)


def angle_symmetric(Q, mw=MW_MID):
    """Optimal symmetric bifurcation half-angle θ (Theorem 9)."""
    r0 = r_star(Q, mw)
    r1 = r_star(Q / 2, mw)
    cos_theta = phi_at_optimum(r0, mw) / (2 * phi_at_optimum(r1, mw))
    return np.degrees(np.arccos(cos_theta))


def T_wall_ratio(N, p_val=p):
    """T_wall(N) / T_wall(2) for N-selection (Theorem 10)."""
    exp = (1 - p_val) / (5 + p_val)
    return (1 + N**exp) / (1 + 2**exp)


# ═══════════════════════════════════════════════════════════════
# Compute all Paper 1 results
# ═══════════════════════════════════════════════════════════════

def compute_all():
    results = {}

    # Table 2: alpha* and r* for three m_w values
    for mw_val, label in [(MW_LOW, 'low'), (MW_MID, 'mid'), (MW_HIGH, 'high')]:
        a = alpha_star(Q0_coronary, mw_val)
        r = r_star(Q0_coronary, mw_val) * 1e3  # mm
        results[f'alpha_mw_{label}'] = a
        results[f'rstar_mw_{label}'] = r

    # Murray limit
    results['alpha_murray'] = alpha_star(Q0_coronary, mw=1e-10)

    # Bounds
    results['bound_lower'] = wall_alpha_limit(p)   # (5+p)/2
    results['bound_upper'] = murray_alpha()         # 3.0

    # Classification table
    results['class_impedance'] = classification_alpha(0)     # 2.0
    results['class_davinci']   = classification_alpha(1)     # 2.5
    results['class_wall']      = classification_alpha(1 + p) # 2.885
    results['class_murray']    = classification_alpha(2)     # 3.0

    # Angle bounds
    # Murray limit angle
    cos_M = 2**(-1/3)
    results['angle_murray_half'] = np.degrees(np.arccos(cos_M))
    results['angle_murray_full'] = 2 * results['angle_murray_half']

    # Wall limit angle
    cos_W = 2**((p - 3) / (5 + p))
    results['angle_wall_half'] = np.degrees(np.arccos(cos_W))
    results['angle_wall_full'] = 2 * results['angle_wall_half']

    # Actual angle for mid m_w
    theta_mid = angle_symmetric(Q0_coronary, MW_MID)
    results['angle_mid_half'] = theta_mid
    results['angle_mid_full'] = 2 * theta_mid

    # Scale dependence over 4 decades
    Q_range = np.logspace(-8, -4, 50)
    alphas = [alpha_star(Q, MW_MID) for Q in Q_range]
    results['alpha_scale_min'] = min(alphas)
    results['alpha_scale_max'] = max(alphas)
    results['alpha_scale_range'] = max(alphas) - min(alphas)

    # Delta-alpha (the corrected abstract claim)
    alpha_best = results['alpha_mw_mid']
    results['delta_alpha_murray'] = murray_alpha() - 2.7
    results['delta_alpha_ours'] = alpha_best - 2.7

    # N-selection
    results['N2_exponent'] = (1 - p) / (5 + p)
    for N in [2, 3, 4, 5]:
        results[f'T_ratio_N{N}'] = T_wall_ratio(N)

    # Cross-network: pulmonary
    results['alpha_pulmonary_limit'] = wall_alpha_limit(p_pulmonary)

    return results


# ═══════════════════════════════════════════════════════════════
# Generate LaTeX dynamic variables
# ═══════════════════════════════════════════════════════════════

def write_tex(results, path):
    lines = [
        "% DYNAMICALLY GENERATED — DO NOT EDIT MANUALLY",
        "% Source: shared/scripts/compute_paper_murray.py",
        f"% Generated: {np.datetime64('today', 'D')}",
        "",
        "% Physical parameters",
        f"\\newcommand{{\\ParamMu}}{{{MU*1e3:.1f}}}                    % mPa·s",
        f"\\newcommand{{\\ParamB}}{{{b:.0f}}}                          % W/m³",
        f"\\newcommand{{\\ParamMwLow}}{{{MW_LOW/1e3:.0f}}}             % kW/m3",
        f"\\newcommand{{\\ParamMwMid}}{{{MW_MID/1e3:.0f}}}             % kW/m3",
        f"\\newcommand{{\\ParamMwHigh}}{{{MW_HIGH/1e3:.0f}}}           % kW/m3",
        f"\\newcommand{{\\ParamCzero}}{{{c0}}}",
        f"\\newcommand{{\\ParamP}}{{{p}}}",
        f"\\newcommand{{\\ParamQ}}{{{Q0_coronary*1e6:.1f}}}            % mL/s",
        f"\\newcommand{{\\ParamRzero}}{{{r0_coronary*1e3:.1f}}}        % mm",
        "",
        "% Table 2: Predicted alpha*",
        f"\\newcommand{{\\AlphaStarLow}}{{{results['alpha_mw_low']:.3f}}}",
        f"\\newcommand{{\\AlphaStarMid}}{{{results['alpha_mw_mid']:.3f}}}",
        f"\\newcommand{{\\AlphaStarHigh}}{{{results['alpha_mw_high']:.3f}}}",
        f"\\newcommand{{\\RstarLow}}{{{results['rstar_mw_low']:.3f}}}",
        f"\\newcommand{{\\RstarMid}}{{{results['rstar_mw_mid']:.3f}}}",
        f"\\newcommand{{\\RstarHigh}}{{{results['rstar_mw_high']:.3f}}}",
        f"\\newcommand{{\\AlphaMurrayLimit}}{{{results['alpha_murray']:.6f}}}",
        "",
        "% Bounds",
        f"\\newcommand{{\\BoundLower}}{{{results['bound_lower']:.3f}}}",
        f"\\newcommand{{\\BoundUpper}}{{{results['bound_upper']:.3f}}}",
        "",
        "% Classification table",
        f"\\newcommand{{\\ClassImpedance}}{{{results['class_impedance']:.3f}}}",
        f"\\newcommand{{\\ClassDaVinci}}{{{results['class_davinci']:.3f}}}",
        f"\\newcommand{{\\ClassWall}}{{{results['class_wall']:.3f}}}",
        f"\\newcommand{{\\ClassMurray}}{{{results['class_murray']:.3f}}}",
        "",
        "% Angle bounds",
        f"\\newcommand{{\\AngleMurrayFull}}{{{results['angle_murray_full']:.1f}}}",
        f"\\newcommand{{\\AngleWallFull}}{{{results['angle_wall_full']:.1f}}}",
        "",
        "% Delta-alpha",
        f"\\newcommand{{\\DeltaAlphaMurray}}{{{results['delta_alpha_murray']:.2f}}}",
        f"\\newcommand{{\\DeltaAlphaOurs}}{{{abs(results['delta_alpha_ours']):.2f}}}",
        "",
        "% Scale range",
        f"\\newcommand{{\\AlphaScaleRange}}{{{results['alpha_scale_range']:.3f}}}",
        "",
        "% Cross-network",
        f"\\newcommand{{\\AlphaPulmonaryLimit}}{{{results['alpha_pulmonary_limit']:.3f}}}",
    ]
    
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    
    print(f"  [OK] Wrote {path}")


# ═══════════════════════════════════════════════════════════════
# Generate alpha-scale figure
# ═══════════════════════════════════════════════════════════════

def generate_figure(fig_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    plt.rcParams.update({
        "font.family": "serif", "font.size": 12,
        "axes.labelsize": 12, "legend.fontsize": 10,
    })

    fig, ax = plt.subplots(figsize=(6, 4.5))
    Q_range = np.logspace(-8, -4, 200)

    for mw_val, color, label in [
        (MW_LOW,  '#2ca02c', f'$m_w = {MW_LOW/1e3:.0f}$ kW/m3'),
        (MW_MID,  '#1f77b4', f'$m_w = {MW_MID/1e3:.0f}$ kW/m3'),
        (MW_HIGH, '#d62728', f'$m_w = {MW_HIGH/1e3:.0f}$ kW/m3'),
    ]:
        alphas = [alpha_star(Q, mw_val) for Q in Q_range]
        ax.plot(Q_range * 1e6, alphas, color=color, lw=2.5, label=label)

    # Asymptotic limits
    mu = murray_alpha()
    wl = wall_alpha_limit()
    ax.axhline(mu, color='#555', ls=':', lw=1.2)
    ax.axhline(wl, color='#555', ls='--', lw=1.2)
    ax.text(1.1e2, mu + 0.0015, r'Murray ($\alpha=3$)', fontsize=9, color='#555', ha='right')
    ax.text(1.1e-2, wl - 0.0025, f'Wall limit ({wl:.3f})', fontsize=9, color='#555',
            va='top')

    # Experimental range: shade the overlap zone and annotate below the wall limit
    ax.axhspan(2.865, 2.885, alpha=0.20, color='#2ca02c')  # thin green band at wall limit
    ax.text(0.99, 0.04,
            r'$\alpha_{\exp} = 2.70 \pm 0.20$ lies below',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            color='#336633',
            bbox=dict(boxstyle='round,pad=0.2', fc='#e8f5e9', alpha=0.8))

    ax.set_xscale('log')
    ax.set_xlabel('Flow $Q$ (mL/s)')
    ax.set_ylabel(r'Local branching exponent $\alpha^*(Q)$')
    ax.set_title(r'Local Optimum: $m_w$ sensitivity')
    ax.set_ylim(2.865, 3.015)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"  [OK] Wrote {fig_path}")


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 60)
    print("Paper 1: Beyond Murray's Law (Generalized Cost Functions)")
    print("=" * 60)

    results = compute_all()

    # Print verification table
    print("\n-- Table 2 verification --")
    print(f"  {'m_w':>10s}  {'r* (mm)':>10s}  {'alpha*':>10s}  {'Delta_a':>10s}")
    for label, mw_name in [('5 kW/m3', 'low'), ('20 kW/m3', 'mid'), ('35 kW/m3', 'high')]:
        a = results[f'alpha_mw_{mw_name}']
        r = results[f'rstar_mw_{mw_name}']
        print(f"  {label:>10s}  {r:>10.3f}  {a:>10.4f}  {a-2.7:>+10.3f}")
    print(f"  {'Murray':>10s}  {'---':>10s}  {results['alpha_murray']:>10.6f}")

    print(f"\n-- Bounds: ({results['bound_lower']:.3f}, {results['bound_upper']:.3f})")
    print(f"-- Angles: {results['angle_murray_full']:.1f} < 2theta < {results['angle_wall_full']:.1f}")
    print(f"-- Delta-alpha: Murray={results['delta_alpha_murray']:.2f}, Ours={results['delta_alpha_ours']:.3f}")
    print(f"-- Scale range: {results['alpha_scale_range']:.4f} over 4 decades")

    # Write LaTeX
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(current_dir))
    tex_path = os.path.join(project_root, 'paper1-murray', 'manuscript', 'dynamic_variables.tex')
    write_tex(results, tex_path)

    # Generate figure
    fig_path = os.path.join(project_root, 'paper1-murray', 'manuscript', 'figures', 'fig_alpha_scale.pdf')
    generate_figure(fig_path)

    print("\n[OK] Paper 1 complete.")
