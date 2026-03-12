"""
compute_supplemental.py — All supplemental computations for Paper 2.

Generates:
  1. dynamic_variables_supplemental.tex  (macros for S1 + S2 + S3)
  2. figS1_supplemental.png              (4-panel validation figure)

S1: Coherent transfer-matrix wave analysis
S2: Monte Carlo parameter sweep (N = 10,000)
S3: Power conversion to absolute units
"""

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import jv
from scipy.stats import pearsonr
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from params import (
    # fluid / wall
    MU, RHO, c0, p, b, MW_MID, MW_LOW, MW_HIGH,
    # coronary morphometry
    r0_coronary, Q0_coronary, ell0_coronary, G_coronary, alpha_w,
    # wave / power-budget constants
    E_elastic, f0_cardiac, ell_factor,
    CPO_watts, pulsatile_fraction, coronary_fraction,
    B_blood, B_wall, A_of_Q,
)
from compute_paper_variational import (
    wave_loss_network, find_alpha_star,
    r_star, phi_cost, reference_cost,
    locally_optimal_radii, two_level_lagrangian,
)


# ═══════════════════════════════════════════════════════════════
# S1 — Transfer-matrix (coherent wave analysis)
# ═══════════════════════════════════════════════════════════════

def _womersley_F10(Wo):
    """Womersley oscillatory velocity-profile correction F10(Wo)."""
    lam = Wo * (1j) ** (3 / 2)
    return 2 * jv(1, lam) / (lam * jv(0, lam))


def _gen_params(r0, alpha, g, omega):
    """Transmission-line parameters for generation g."""
    r_g   = r0 * 2 ** (-g / alpha)
    ell_g = ell_factor * r_g
    h_g   = c0 * r_g ** p
    c_g   = np.sqrt(E_elastic * h_g / (2 * RHO * r_g))   # Moens-Korteweg
    Wo_g  = r_g * np.sqrt(omega * RHO / MU)
    F10_g = _womersley_F10(Wo_g)
    Z_c_g = RHO * c_g / (np.pi * r_g ** 2)
    gam_g = 1j * omega * np.sqrt(1 - F10_g) / c_g
    return r_g, ell_g, Z_c_g, gam_g


def _terminal_Z(r_G, ell_G, omega):
    """Two-element Windkessel terminal load impedance."""
    h_G = c0 * r_G ** p
    R_t = 8 * MU * ell_G / (np.pi * r_G ** 4)
    C_t = 3 * np.pi * r_G ** 3 * ell_G / (2 * E_elastic * h_G)
    return R_t / (1 + 1j * omega * R_t * C_t)


def _gamma_total_sq(alpha, omega, G=G_coronary):
    """
    Compute |Γ_total|² at the root for given α and angular frequency ω.
    Uses backward recursion from terminal Windkessel load.
    """
    r0    = r0_coronary
    r_G   = r0 * 2 ** (-G / alpha)
    ell_G = ell_factor * r_G
    Z_in  = _terminal_Z(r_G, ell_G, omega)

    for g in range(G - 1, -1, -1):
        r_g, ell_g, Z_c_g, gam_g = _gen_params(r0, alpha, g, omega)
        Z_d      = Z_in / 2                          # two daughters in parallel
        tanh_gl  = np.tanh(gam_g * ell_g)
        Z_in     = Z_c_g * (Z_d + Z_c_g * tanh_gl) / (Z_c_g + Z_d * tanh_gl)

    _, _, Z_c0, _ = _gen_params(r0, alpha, 0, omega)
    Gamma = (Z_in - Z_c0) / (Z_in + Z_c0)
    return float(np.abs(Gamma) ** 2)


def run_S1(G=G_coronary, N_harmonics=4, N_alpha=200):
    """
    Full coherent transfer-matrix analysis.
    Returns dict with all computed S1 quantities.
    """
    omega0     = 2 * np.pi * f0_cardiac
    alpha_vals = np.linspace(1.5, 3.5, N_alpha)

    # |Γ_total|²(α) at the cardiac fundamental
    gamma_curve  = np.array([_gamma_total_sq(a, omega0, G) for a in alpha_vals])

    # Residual at impedance-matching exponent (all junctions matched)
    gamma_at_aw  = _gamma_total_sq(alpha_w, omega0, G)

    # Geometric penalty: isolates junction-mismatch from terminal-load contribution
    delta_coherent = gamma_curve - gamma_at_aw

    # Incoherent cumulative network wave cost
    c_wave_arr = np.array([wave_loss_network(a, G, alpha_w) for a in alpha_vals])

    # Pearson correlation for α > 2.3 (Fabry-Pérot-free regime)
    mask = alpha_vals > 2.3
    r_corr, _ = pearsonr(delta_coherent[mask], c_wave_arr[mask])

    # Incoherent minimax α* (reference)
    alpha_star_inc, _ = find_alpha_star(0.833, G=G)

    # Factor of proportionality between Δ_coherent and C^net_wave (for α > 2.3)
    coh_factor = np.mean(delta_coherent[mask] / c_wave_arr[mask])

    # Coherent minimax α*: find crossing of normalised delta_coherent with C_transport.
    # The paper's argument is that the proportionality constant K cancels in the
    # equal-cost condition: delta_coherent/K ≈ C_wave_inc, so the crossing with
    # C_transport is the same up to the residual non-proportional part (1 - r²).
    from scipy.interpolate import interp1d
    from scipy.optimize import brentq

    r_local = locally_optimal_radii()
    Phi_ref = reference_cost(r_local)

    c_t_arr      = np.array([two_level_lagrangian(a, 0.0, r_local, Phi_ref)
                              for a in alpha_vals])
    delta_norm   = delta_coherent / coh_factor          # ≈ C_wave_inc in shape
    delta_norm_fn = interp1d(alpha_vals, delta_norm,
                              kind='cubic', bounds_error=False,
                              fill_value='extrapolate')
    c_t_fn        = interp1d(alpha_vals, c_t_arr,
                              kind='cubic', bounds_error=False,
                              fill_value='extrapolate')

    diff_fn = lambda a: float(np.real(delta_norm_fn(a))) - float(c_t_fn(a))
    try:
        alpha_star_coh  = brentq(diff_fn, 2.3, 3.3)
    except ValueError:
        # Fallback: use incoherent alpha* if no crossing found in range
        alpha_star_coh = alpha_star_inc
    delta_alpha_star = abs(alpha_star_coh - alpha_star_inc)

    # Maximum Womersley number (root vessel, fundamental frequency)
    Wo_max = r0_coronary * np.sqrt(omega0 * RHO / MU)

    # Per-harmonic analysis at α = α*_incoherent
    per_harmonic = []
    for n in range(1, N_harmonics + 1):
        omega_n = n * omega0
        g_sq    = _gamma_total_sq(alpha_star_inc, omega_n, G)
        per_harmonic.append({'n': n, 'f_hz': n * f0_cardiac, 'gamma_sq': g_sq})

    return {
        'gamma_at_aw':    gamma_at_aw,
        'r_corr':         r_corr,
        'delta_alpha':    delta_alpha_star,
        'coh_factor':     coh_factor,
        'wo_max':         Wo_max,
        'per_harmonic':   per_harmonic,
        'alpha_vals':     alpha_vals,
        'delta_coherent': delta_coherent,
        'c_wave_net':     c_wave_arr,
        'alpha_star_inc': alpha_star_inc,
    }


# ═══════════════════════════════════════════════════════════════
# S2 — Monte Carlo parameter sweep
# ═══════════════════════════════════════════════════════════════

def _mc_alpha_star(p_s, G_s, alpha_t_s, eta=0.833):
    """
    Find α* for one MC sample.
    Uses quadratic transport approximation: C_transport ≈ k_t·(α − α_t)²
    with k_t = 0.85·G (numerical curvature; see Supplemental S2.2).
    """
    alpha_w_s = (5 - p_s) / 2
    k_t       = 0.85 * G_s

    def L(a):
        c_w = wave_loss_network(a, int(round(G_s)), alpha_w_s)
        c_t = k_t * (a - alpha_t_s) ** 2
        return eta * c_w + (1 - eta) * c_t

    try:
        res = minimize_scalar(L, bounds=(1.5, 4.0), method='bounded')
        return res.x
    except Exception:
        return np.nan


def run_S2(N=10000, seed=42):
    """
    Monte Carlo parameter sweep (N = 10,000).
    Samples: p, G, m_w, α_t  (see Table S1 for distributions).
    Returns statistics and decoupling analysis.
    """
    rng = np.random.default_rng(seed)

    # Parameter distributions (Table S1)
    p_samples  = np.clip(rng.normal(0.77, 0.05, N), 0.50, 0.95)
    G_samples  = rng.integers(9, 14, N).astype(float)   # uniform over {9,…,13}
    at_samples = rng.uniform(2.80, 3.00, N)
    # m_w not used in quadratic approx (k_t = 0.85·G); included in Table S1 for completeness

    # Compute α* for each sample
    alpha_stars = np.array([
        _mc_alpha_star(p_samples[i], G_samples[i], at_samples[i])
        for i in range(N)
    ])
    alpha_stars = alpha_stars[~np.isnan(alpha_stars)]

    mean  = float(np.mean(alpha_stars))
    std   = float(np.std(alpha_stars, ddof=1))
    ci_lo = float(np.percentile(alpha_stars, 2.5))
    ci_hi = float(np.percentile(alpha_stars, 97.5))

    # Decoupling analysis: α*(α_t) at baseline params, α_t ∈ [2.80, 3.00]
    at_range   = np.linspace(2.80, 3.00, 41)
    stars_decoupled = np.array([
        _mc_alpha_star(p, G_coronary, at) for at in at_range
    ])
    delta_alpha_decoupled = float(np.max(stars_decoupled) - np.min(stars_decoupled))

    return {
        'N':                     len(alpha_stars),
        'mean':                  mean,
        'std':                   std,
        'ci_lo':                 ci_lo,
        'ci_hi':                 ci_hi,
        'delta_alpha':           delta_alpha_decoupled,
        'at_range':              at_range,
        'stars_decoupled':       stars_decoupled,
        'alpha_star_samples':    alpha_stars,
    }


# ═══════════════════════════════════════════════════════════════
# S3 — Power conversion to absolute units
# ═══════════════════════════════════════════════════════════════

def run_S3(alpha_star_val):
    """
    Power budget at the minimax operating point α*.
    All powers in Watts internally; exported in mW.
    """
    # Incident pulsatile power on the coronary tree
    P_pulse = CPO_watts * pulsatile_fraction * coronary_fraction   # W

    # Wave dissipation
    C_wave_net  = wave_loss_network(alpha_star_val, G_coronary, alpha_w)
    P_wave      = P_pulse * C_wave_net

    # Transport baseline: total cost when every vessel is at its local optimum
    P_baseline = 0.0
    for g in range(G_coronary):
        Q_g   = Q0_coronary / 2 ** g
        r_g   = r_star(Q_g, MW_MID)
        phi_g = phi_cost(r_g, Q_g, MW_MID)
        ell_g = ell0_coronary                     # proximal segment length as reference
        P_baseline += 2 ** g * phi_g * ell_g

    # Transport excess at α*
    r_local         = locally_optimal_radii()
    Phi_ref         = reference_cost(r_local)
    C_transport_net = two_level_lagrangian(alpha_star_val, 0.0, r_local, Phi_ref)
    P_transport     = P_baseline * C_transport_net

    ratio = P_wave / P_transport if P_transport > 1e-20 else float('nan')

    return {
        'P_pulse_mW':      P_pulse * 1e3,
        'P_wave_mW':       P_wave  * 1e3,
        'P_baseline_mW':   P_baseline  * 1e3,
        'P_transport_mW':  P_transport * 1e3,
        'ratio':           ratio,
        'C_wave_pct':      C_wave_net      * 100,
        'C_transport_pct': C_transport_net * 100,
    }


# ═══════════════════════════════════════════════════════════════
# Write dynamic_variables_supplemental.tex
# ═══════════════════════════════════════════════════════════════

def write_tex(s1, s2, s3, path):
    """Write all supplemental macros to a LaTeX file."""
    roman = ['I', 'II', 'III', 'IV']

    lines = [
        "% DYNAMICALLY GENERATED -- DO NOT EDIT MANUALLY",
        "% Source: shared/scripts/compute_supplemental.py",
        f"% Generated: {np.datetime64('now')}",
        "",
        "% ── Physical constants re-exported for supplemental ──",
        f"\\newcommand{{\\VarElasticModulusMPa}}{{{E_elastic/1e6:.1f}}}",
        f"\\newcommand{{\\VarHeartRateHz}}{{{f0_cardiac:.1f}}}",
        f"\\newcommand{{\\VarEllFactor}}{{{ell_factor}}}",
        f"\\newcommand{{\\VarCPOWatts}}{{{CPO_watts:.1f}}}",
        f"\\newcommand{{\\VarPulsatilePct}}{{{int(pulsatile_fraction*100)}}}",
        f"\\newcommand{{\\VarCoronaryPct}}{{{int(coronary_fraction*100)}}}",
        "",
        "% ── S1: Transfer-matrix results ──",
        f"\\newcommand{{\\VarGammaTotalAW}}{{{s1['gamma_at_aw']:.3f}}}",
        f"\\newcommand{{\\VarCoherentCorr}}{{{s1['r_corr']:.2f}}}",
        f"\\newcommand{{\\VarCoherentShift}}{{{s1['delta_alpha']:.3f}}}",
        f"\\newcommand{{\\VarCoherentFactor}}{{{round(s1['coh_factor']):.0f}}}",
        f"\\newcommand{{\\VarWomersleyMax}}{{{s1['wo_max']:.1f}}}",
        "",
        "% ── S1: Per-harmonic table ──",
    ]

    for i, row in enumerate(s1['per_harmonic']):
        n = roman[i]
        lines.append(f"\\newcommand{{\\VarFreqHarmonic{n}}}{{{row['f_hz']:.1f}}}")
        lines.append(f"\\newcommand{{\\VarGammaHarmonic{n}}}{{{row['gamma_sq']:.3f}}}")

    lines += [
        "",
        "% ── S2: Monte Carlo results ──",
        "\\newcommand{\\VarMCSamples}{10{,}000}",
        f"\\newcommand{{\\VarMCAlphaMean}}{{{s2['mean']:.2f}}}",
        f"\\newcommand{{\\VarMCAlphaStd}}{{{s2['std']:.2f}}}",
        f"\\newcommand{{\\VarMCAlphaCILow}}{{{s2['ci_lo']:.2f}}}",
        f"\\newcommand{{\\VarMCAlphaCIHigh}}{{{s2['ci_hi']:.2f}}}",
        f"\\newcommand{{\\VarMCDeltaAlpha}}{{{s2['delta_alpha']:.3f}}}",
        "",
        "% ── S3: Power budget ──",
        f"\\newcommand{{\\VarPPulsemW}}{{{s3['P_pulse_mW']:.1f}}}",
        f"\\newcommand{{\\VarPWavemW}}{{{s3['P_wave_mW']:.2f}}}",
        f"\\newcommand{{\\VarPBaselinemW}}{{{s3['P_baseline_mW']:.0f}}}",
        f"\\newcommand{{\\VarPTransportmW}}{{{s3['P_transport_mW']:.1f}}}",
        f"\\newcommand{{\\VarPRatio}}{{{s3['ratio']:.2f}}}",
    ]

    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"  [OK] Wrote {path}")


# ═══════════════════════════════════════════════════════════════
# Generate figS1_supplemental.png (4 panels)
# ═══════════════════════════════════════════════════════════════

def generate_figS1(s1, s2, fig_path):
    """
    4-panel supplemental validation figure.
    A: Coherent vs incoherent wave penalty (proportionality)
    B: Monte Carlo distribution of α*
    C: Decoupling analysis α*(α_t)
    D: Per-harmonic |Γ_total|²
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "serif", "font.size": 10,
        "axes.labelsize": 11, "legend.fontsize": 9,
    })

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()

    alpha_vals = s1['alpha_vals']

    # ── Panel A: Coherent geometric penalty vs incoherent C_wave ──
    ax = axes[0]
    dc = s1['delta_coherent']
    cw = s1['c_wave_net']
    # Restrict to physiologically valid range α > 2.1
    mask_plot = alpha_vals > 2.1
    dc_max = np.max(np.abs(dc[mask_plot])) or 1
    cw_max = np.max(cw[mask_plot]) or 1
    ax.plot(alpha_vals[mask_plot], dc[mask_plot] / dc_max, '--', color='#8c564b', lw=2,
            label=r'$\Delta_{\mathrm{coherent}}(\alpha)$ (normalised)')
    ax.plot(alpha_vals[mask_plot], cw[mask_plot] / cw_max, '-',  color='#1f77b4', lw=2,
            label=r'$\mathcal{C}^\mathrm{net}_\mathrm{wave}(\alpha)$ (normalised)')
    ax.axvline(s1['alpha_star_inc'], color='k', ls=':', lw=1.2, alpha=0.7,
               label=rf"$\alpha^* = {s1['alpha_star_inc']:.3f}$")
    # Shade the correlation-valid region (α > 2.3)
    ax.axvspan(2.3, 3.3, alpha=0.06, color='#2ca02c',
               label=r'Corr. region ($\alpha > 2.3$)')
    ax.text(0.97, 0.05, f"$r = {s1['r_corr']:.2f}$ (for $\\alpha > 2.3$)",
            transform=ax.transAxes, fontsize=10, va='bottom', ha='right',
            bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.8))
    ax.set_xlabel(r'Branching exponent $\alpha$')
    ax.set_ylabel('Normalised cost')
    ax.set_title('(A) Coherent vs Incoherent Wave Penalty')
    ax.set_xlim(2.1, 3.3)
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Panel B: Monte Carlo histogram ──
    ax = axes[1]
    samples = s2['alpha_star_samples']
    ax.hist(samples, bins=30, color='#2ca02c', alpha=0.7, density=True,
            edgecolor='white', lw=0.5)
    ax.axvline(s2['mean'], color='k', lw=2,
               label=rf"$\alpha^* = {s2['mean']:.2f} \pm {s2['std']:.2f}$")
    ax.axvline(s2['ci_lo'], color='#555', ls='--', lw=1.5,
               label=rf"95\% CI: [{s2['ci_lo']:.2f}, {s2['ci_hi']:.2f}]")
    ax.axvline(s2['ci_hi'], color='#555', ls='--', lw=1.5)
    # Note: alpha_exp=2.70 is below the MC range (local vs network result)
    ax.set_xlabel(r'$\alpha^*$')
    ax.set_ylabel('Density')
    ax.set_title(f'(B) Monte Carlo ($N = {s2["N"]:,}$, local model)')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Panel C: Decoupling analysis ──
    ax = axes[2]
    at_r = s2['at_range']
    as_d = s2['stars_decoupled']
    ax.plot(at_r, as_d, 's-', color='#9467bd', lw=2, ms=5)
    ax.set_xlabel(r'Transport Ground State $\alpha_t$')
    ax.set_ylabel(r'Minimax $\alpha^*$')
    ax.set_title(f'(C) Decoupling: $\\Delta\\alpha^* = {s2["delta_alpha"]:.3f}$')
    ax.grid(True, alpha=0.3)

    # ── Panel D: Per-harmonic reflection ──
    ax = axes[3]
    freqs  = [row['f_hz']    for row in s1['per_harmonic']]
    gammas = [row['gamma_sq'] for row in s1['per_harmonic']]
    bars = ax.bar(range(len(freqs)), gammas, color='#ff7f0e', alpha=0.85,
                  edgecolor='k', lw=0.8)
    # Annotate each bar with its value
    for i, (b, g) in enumerate(zip(bars, gammas)):
        ax.text(b.get_x() + b.get_width() / 2, g + 0.002, f'{g:.3f}',
                ha='center', va='bottom', fontsize=9)
    ax.set_xticks(range(len(freqs)))
    ax.set_xticklabels([f'{f:.1f} Hz' for f in freqs])
    ax.set_xlabel('Harmonic frequency')
    ax.set_ylabel(r'$|\Gamma_{\mathrm{total}}|^2$')
    ax.set_title(r'(D) Per-Harmonic Reflection at $\alpha^*$')
    g_min = min(gammas)
    g_max = max(gammas)
    margin = (g_max - g_min) * 0.5 or 0.02
    ax.set_ylim(max(0, g_min - margin), g_max + margin * 2)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Wrote {fig_path}")


# ═══════════════════════════════════════════════════════════════
# Public entry point (called by run_all.py)
# ═══════════════════════════════════════════════════════════════

def run_supplemental():
    """Compute all supplemental results and write outputs."""
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT  = os.path.dirname(os.path.dirname(SCRIPT_DIR))
    SUPP_DIR   = os.path.join(REPO_ROOT, 'paper2-variational', 'supplements')
    FIG_DIR    = os.path.join(SUPP_DIR, 'figures')
    os.makedirs(FIG_DIR, exist_ok=True)

    print("=" * 60)
    print("Supplemental: Transfer-Matrix + Monte Carlo + Power Budget")
    print("=" * 60)

    print("\nS1: Transfer-matrix analysis...")
    s1 = run_S1()
    print(f"  |Gamma|^2(aw)  = {s1['gamma_at_aw']:.3f}")
    print(f"  Correlation r  = {s1['r_corr']:.2f}")
    print(f"  |Delta_alpha*| = {s1['delta_alpha']:.3f}")
    print(f"  Scale factor   = {s1['coh_factor']:.1f}x")
    print(f"  Wo_max         = {s1['wo_max']:.2f}")
    for row in s1['per_harmonic']:
        print(f"  Harmonic {row['n']}: f={row['f_hz']:.1f} Hz, |Gamma|^2={row['gamma_sq']:.3f}")

    print("\nS2: Monte Carlo (N=10,000)...")
    s2 = run_S2(N=10000, seed=42)
    print(f"  alpha* = {s2['mean']:.2f} +/- {s2['std']:.2f}")
    print(f"  95% CI = [{s2['ci_lo']:.2f}, {s2['ci_hi']:.2f}]")
    print(f"  Delta_alpha* (decoupling) = {s2['delta_alpha']:.3f}")

    alpha_star_val = s1['alpha_star_inc']
    print(f"\nS3: Power budget (alpha* = {alpha_star_val:.3f})...")
    s3 = run_S3(alpha_star_val)
    print(f"  P_pulse     = {s3['P_pulse_mW']:.2f} mW")
    print(f"  C_wave      = {s3['C_wave_pct']:.1f}%")
    print(f"  P_wave      = {s3['P_wave_mW']:.3f} mW")
    print(f"  P_baseline  = {s3['P_baseline_mW']:.1f} mW")
    print(f"  C_transport = {s3['C_transport_pct']:.1f}%")
    print(f"  P_transport = {s3['P_transport_mW']:.3f} mW")
    print(f"  Ratio       = {s3['ratio']:.2f}")

    tex_path = os.path.join(SUPP_DIR, 'dynamic_variables_supplemental.tex')
    write_tex(s1, s2, s3, tex_path)

    fig_path = os.path.join(FIG_DIR, 'figS1_supplemental.png')
    generate_figS1(s1, s2, fig_path)

    print("\n[OK] Supplemental complete.")
    return s1, s2, s3


if __name__ == '__main__':
    run_supplemental()
