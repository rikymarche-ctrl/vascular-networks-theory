"""
compute.py — Paper I: Beyond Murray's Law
==========================================
Generates dynamic_variables.tex containing all LaTeX \\newcommand definitions
used in the manuscript, together with the manuscript figure. Every calculated
quantity in the .tex file originates here; no numerical values are hardcoded
in the source.

Model: Three-term single-vessel cost minimisation
--------------------------------------------------
The metabolic cost per unit length of a cylindrical vessel of radius r
carrying volumetric flow Q is:

    Φ(r, Q) = (8μQ²/πr⁴)  +  bπr²  +  2π m_w c_0 r^{1+p}
               viscous          blood       wall tissue

where the wall term is derived from integrating the volumetric metabolic
rate m_w over the thin-wall cross-section 2πr·h(r), with h(r) = c_0 r^p
taken from histological measurements (Rhodin 1967, p ≈ 0.77).

Strict convexity of Φ guarantees a unique optimal radius r*(Q) at every
flow (Theorem 1). The branching exponent α*(Q) is defined by the symmetric
bifurcation condition r*(Q)^α = 2·r*(Q/2)^α and satisfies the strict
bounds (5+p)/2 < α*(Q) < 3 for all Q > 0 (Theorem 4).

All parameters are drawn from independent literature sources; no value is
fitted to morphometric branching-exponent data.

Empirical inputs (literature sources, not derived here)
-------------------------------------------------------
    μ    = 3.5 mPa·s        blood dynamic viscosity        (Caro 1978)
    b    = 1500 W/m³         blood-volume metabolic cost    (Murray 1926; Taber 1998)
    m_w  ∈ [5, 35] kW/m³    wall tissue metabolic rate     (Paul 1980)
    c_0  = 0.041 m^{1−p}    wall-thickness prefactor       (Rhodin 1967)
    p    = 0.77              wall-thickness exponent        (Rhodin 1967; Kassab 1993)
    Q_0  = 1.3 mL/s          proximal coronary flow         (Kassab 1993)
    r_0  = 1.5 mm            proximal coronary radius       (Kassab 1993)

Outputs
-------
    ../manuscript/dynamic_variables.tex
    ../manuscript/figures/fig_alpha_scale.pdf
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

def format_macro(name, value, unit, description):
    """Helper to format LaTeX macro definitions consistently."""
    macro = f"\\newcommand{{\\{name}}}{{{value}}}"
    return f"{macro:<50s} % [{unit:<8}] {description}"

def write_tex(results, path):
    lines = [
        "% ====================================================================",
        "% DYNAMIC VARIABLES",
        "% ====================================================================",
        "% Source: paper1-murray/scripts/compute.py",
        f"% Generated: {np.datetime64('today', 'D')}",
        "% ====================================================================\n",
        
        "% --------------------------------------------------------------------",
        "% 1. PHYSICAL PARAMETERS",
        "% --------------------------------------------------------------------",
        format_macro("ParamMu", f"{MU*1e3:.1f}", "mPa s", "Dynamic viscosity of blood"),
        format_macro("ParamB", f"{b:.0f}", "W/m³", "Metabolic cost per unit blood volume"),
        format_macro("ParamMwLow", f"{MW_LOW/1e3:.0f}", "kW/m³", "Lower bound: Metabolic rate of wall tissue"),
        format_macro("ParamMwMid", f"{MW_MID/1e3:.0f}", "kW/m³", "Reference: Metabolic rate of wall tissue"),
        format_macro("ParamMwHigh", f"{MW_HIGH/1e3:.0f}", "kW/m³", "Upper bound: Metabolic rate of wall tissue"),
        format_macro("ParamCzero", f"{c0}", "m^{1-p}", "Histological wall thickness pre-factor"),
        format_macro("ParamP", f"{p}", "-", "Histological wall thickness exponent"),
        format_macro("ParamQ", f"{Q0_coronary*1e6:.1f}", "mL/s", "Inlet volumetric flow rate for coronary artery"),
        format_macro("ParamRzero", f"{r0_coronary*1e3:.1f}", "mm", "Inlet radius for coronary artery"),
        "",
        "% --------------------------------------------------------------------",
        "% 2. PREDICTED EXPONENTS AND GEOMETRY (TABLE 2)",
        "% --------------------------------------------------------------------",
        format_macro("AlphaStarLow", f"{results['alpha_mw_low']:.3f}", "-", "Optimum branching exponent at MW_LOW"),
        format_macro("AlphaStarMid", f"{results['alpha_mw_mid']:.3f}", "-", "Optimum branching exponent at MW_MID"),
        format_macro("AlphaStarHigh", f"{results['alpha_mw_high']:.3f}", "-", "Optimum branching exponent at MW_HIGH"),
        format_macro("RstarLow", f"{results['rstar_mw_low']:.3f}", "mm", "Optimum parent radius at MW_LOW"),
        format_macro("RstarMid", f"{results['rstar_mw_mid']:.3f}", "mm", "Optimum parent radius at MW_MID"),
        format_macro("RstarHigh", f"{results['rstar_mw_high']:.3f}", "mm", "Optimum parent radius at MW_HIGH"),
        format_macro("AlphaMurrayLimit", f"{results['alpha_murray']:.6f}", "-", "Branching exponent in the pure Murray limit ($m_w \\to 0$)"),
        "",
        "% --------------------------------------------------------------------",
        "% 3. THEORETICAL BOUNDS",
        "% --------------------------------------------------------------------",
        format_macro("BoundLower", f"{results['bound_lower']:.3f}", "-", "Strict lower bound for the exponent $\\alpha^* > (5+p)/2$"),
        format_macro("BoundUpper", f"{results['bound_upper']:.3f}", "-", "Strict upper bound for the exponent $\\alpha^* < 3$"),
        "",
        "% --------------------------------------------------------------------",
        "% 4. SINGLE-TERM CLASSIFICATION TABLE",
        "% --------------------------------------------------------------------",
        format_macro("ClassImpedance", f"{results['class_impedance']:.3f}", "-", "Exponent for impedance matching ($\\gamma=0$)"),
        format_macro("ClassDaVinci", f"{results['class_davinci']:.3f}", "-", "Exponent for Da Vinci surface law ($\\gamma=1$)"),
        format_macro("ClassWall", f"{results['class_wall']:.3f}", "-", "Exponent for pure wall cost ($\\gamma=1+p$)"),
        format_macro("ClassMurray", f"{results['class_murray']:.3f}", "-", "Exponent for pure Murray volume cost ($\\gamma=2$)"),
        "",
        "% --------------------------------------------------------------------",
        "% 5. OPTIMAL BIFURCATION ANGLES",
        "% --------------------------------------------------------------------",
        format_macro("AngleMurrayFull", f"{results['angle_murray_full']:.1f}", "deg", "Full bifurcation angle $2\\theta$ in Murray limit"),
        format_macro("AngleWallFull", f"{results['angle_wall_full']:.1f}", "deg", "Full bifurcation angle $2\\theta$ in Wall limit"),
        "",
        "% --------------------------------------------------------------------",
        "% 6. MACROSCOPIC COMPARISONS",
        "% --------------------------------------------------------------------",
        format_macro("DeltaAlphaMurray", f"{results['delta_alpha_murray']:.2f}", "-", "Discrepancy: Murray law vs. Empirical mean (3.0 - 2.7)"),
        format_macro("DeltaAlphaOurs", f"{abs(results['delta_alpha_ours']):.2f}", "-", "Remaining Discrepancy: Our model vs. Empirical mean"),
        format_macro("AlphaScaleRange", f"{results['alpha_scale_range']:.3f}", "-", "Variation of $\\alpha^*$ over 4 decades of flow scale"),
        format_macro("AlphaPulmonaryLimit", f"{results['alpha_pulmonary_limit']:.3f}", "-", "Theoretical wall limit for the pulmonary network"),
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
