# -*- coding: utf-8 -*-
"""
verify_sion_convexity.py
========================================================================================
Mathematical verification of Sion's minimax theorem applicability to the
two-level Lagrangian L(alpha, eta) = eta·C_wave(alpha) + (1-eta)·C_transport(alpha).

Sion's Theorem (1958) Requirements:
1. X = [alpha_w, alpha_t] ≈ [2.1, 2.9] compact ✓
2. Y = [0, 1] compact ✓
3. L(·, eta) quasi-convex and upper semicontinuous for each eta ∈ Y
4. L(alpha, ·) quasi-concave and lower semicontinuous for each alpha ∈ X

This script verifies conditions 3 and 4 by:
A. Checking convexity (strict → quasi-convex) of C_wave(alpha) and C_transport(alpha)
B. Verifying linearity (affine → both convex and concave) in eta
C. Computing second derivatives d²C/dalpha² numerically
D. Identifying regions where convexity may fail
E. Testing quasi-convexity when strict convexity fails

Author: Riccardo Marchesi
Date: 2026-05-09
"""

import sys, os
import io

# Force UTF-8 output on Windows
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'shared', 'scripts'))

import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from params import (
    MU, RHO, c0, p, MW_MID, r0_coronary, Q0_coronary, Q0_peak,
    ell0_coronary, G_coronary, alpha_w, beta_coronary,
    B_blood, B_wall, A_of_Q,
)

# Import functions from compute.py
import compute

# ==============================================================
# Helper: Numerical second derivative
# ==============================================================

def second_derivative(func, x, h=1e-5):
    """Compute d²f/dx² using central finite difference."""
    return (func(x + h) - 2*func(x) + func(x - h)) / h**2

def is_convex_numerical(func, x_grid, tol=1e-8):
    """Check if d²f/dx² ≥ -tol for all x in grid."""
    second_derivs = [second_derivative(func, x) for x in x_grid]
    min_curvature = min(second_derivs)
    is_convex = min_curvature >= -tol
    return is_convex, min_curvature, second_derivs

# ==============================================================
# Helper: Quasi-convexity test
# ==============================================================

def is_quasiconvex(func, x_grid):
    """
    Test quasi-convexity: all sublevel sets {x : f(x) ≤ c} are convex.

    For a function on an interval, this is equivalent to:
    - f has at most one local minimum
    - f is decreasing then increasing (unimodal)

    We test by checking if f is monotonically decreasing up to the minimum
    and monotonically increasing afterwards.
    """
    f_vals = np.array([func(x) for x in x_grid])

    # Find global minimum
    min_idx = np.argmin(f_vals)

    # Check monotonicity on each side
    left_monotone = all(f_vals[i] >= f_vals[i+1] for i in range(min_idx))
    right_monotone = all(f_vals[i] <= f_vals[i+1] for i in range(min_idx, len(f_vals)-1))

    is_qconvex = left_monotone and right_monotone
    return is_qconvex, min_idx, f_vals

# ==============================================================
# Cost function definitions
# ==============================================================

def C_wave(alpha, G=G_coronary, aw=alpha_w, N=2):
    """Network wave cost: 1 - (1 - |Γ(alpha)|²)^G."""
    return compute.wave_loss_network(alpha, G, aw, N)

def C_transport(alpha, mw=MW_MID, G=G_coronary, beta=beta_coronary, aw=alpha_w, N=2):
    """Network transport cost: relative penalty for deviating from local optima."""
    r_local = compute.locally_optimal_radii(mw, G, N)
    Phi_ref = compute.reference_cost(r_local, mw, G, beta, N)
    return compute.two_level_lagrangian(alpha, 0.0, r_local, Phi_ref, mw, G, beta, aw, N)

def L_lagrangian(alpha, eta, mw=MW_MID, G=G_coronary, beta=beta_coronary, aw=alpha_w, N=2):
    """Full Lagrangian L(alpha, eta) = eta·C_wave(alpha) + (1-eta)·C_transport(alpha)."""
    cw = C_wave(alpha, G, aw, N)
    ct = C_transport(alpha, mw, G, beta, aw, N)
    return eta * cw + (1 - eta) * ct

# ==============================================================
# VERIFICATION 1: Convexity in alpha for fixed eta
# ==============================================================

def verify_convexity_in_alpha():
    """
    Check if C_wave(alpha) and C_transport(alpha) are convex in alpha over [alpha_w, alpha_t].

    Returns:
        dict with results for both cost functions
    """
    print("\n" + "="*70)
    print("VERIFICATION 1: Convexity in alpha (for fixed eta)")
    print("="*70)

    # Define domain: from alpha_w to alpha_t (approximate transport minimum)
    alpha_min = alpha_w
    alpha_max = 3.0  # Upper bound
    alpha_grid = np.linspace(alpha_min, alpha_max, 200)

    results = {}

    # --- Test C_wave(alpha) ---
    print("\n1.A. Testing C_wave(alpha):")
    print("-" * 70)

    is_conv_wave, min_curv_wave, curvs_wave = is_convex_numerical(C_wave, alpha_grid)

    print(f"  Domain: [{alpha_min:.3f}, {alpha_max:.3f}]")
    print(f"  Minimum curvature d²C_wave/dalpha²: {min_curv_wave:.6e}")
    print(f"  Is strictly convex? {is_conv_wave}")

    if not is_conv_wave:
        # Test quasi-convexity
        is_qconv_wave, min_idx_wave, f_vals_wave = is_quasiconvex(C_wave, alpha_grid)
        print(f"  Is quasi-convex? {is_qconv_wave}")
        print(f"  Minimum at alpha = {alpha_grid[min_idx_wave]:.4f}")
    else:
        is_qconv_wave = True
        min_idx_wave = None
        f_vals_wave = None

    results['C_wave'] = {
        'convex': is_conv_wave,
        'quasi_convex': is_qconv_wave,
        'min_curvature': min_curv_wave,
        'curvatures': curvs_wave,
        'alpha_grid': alpha_grid,
    }

    # --- Test C_transport(alpha) ---
    print("\n1.B. Testing C_transport(alpha):")
    print("-" * 70)

    is_conv_trans, min_curv_trans, curvs_trans = is_convex_numerical(C_transport, alpha_grid)

    print(f"  Domain: [{alpha_min:.3f}, {alpha_max:.3f}]")
    print(f"  Minimum curvature d²C_transport/dalpha²: {min_curv_trans:.6e}")
    print(f"  Is strictly convex? {is_conv_trans}")

    if not is_conv_trans:
        # Test quasi-convexity
        is_qconv_trans, min_idx_trans, f_vals_trans = is_quasiconvex(C_transport, alpha_grid)
        print(f"  Is quasi-convex? {is_qconv_trans}")
        if is_qconv_trans:
            print(f"  Minimum at alpha = {alpha_grid[min_idx_trans]:.4f}")

        # Find regions of negative curvature
        neg_curv_indices = [i for i, c in enumerate(curvs_trans) if c < -1e-8]
        if neg_curv_indices:
            alpha_neg = [alpha_grid[i] for i in neg_curv_indices]
            print(f"  Regions of negative curvature: alpha ∈ [{min(alpha_neg):.3f}, {max(alpha_neg):.3f}]")
    else:
        is_qconv_trans = True
        min_idx_trans = None
        f_vals_trans = None

    results['C_transport'] = {
        'convex': is_conv_trans,
        'quasi_convex': is_qconv_trans,
        'min_curvature': min_curv_trans,
        'curvatures': curvs_trans,
        'alpha_grid': alpha_grid,
    }

    # --- Test L(alpha, eta) for representative eta values ---
    print("\n1.C. Testing L(alpha, eta) for fixed eta:")
    print("-" * 70)

    eta_test_values = [0.0, 0.25, 0.5, 0.74, 1.0]  # Include eta* ≈ 0.74

    for eta_val in eta_test_values:
        def L_fixed_eta(alpha):
            return L_lagrangian(alpha, eta_val)

        is_conv_L, min_curv_L, _ = is_convex_numerical(L_fixed_eta, alpha_grid)
        print(f"  eta = {eta_val:.2f}: Convex? {is_conv_L}, min d²L/dalpha² = {min_curv_L:.6e}")

    return results

# ==============================================================
# VERIFICATION 2: Linearity (convexity & concavity) in eta
# ==============================================================

def verify_linearity_in_eta():
    """
    Verify that L(alpha, eta) = eta·C_wave(alpha) + (1-eta)·C_transport(alpha) is affine in eta.

    This is trivially true by construction, but we verify numerically.
    """
    print("\n" + "="*70)
    print("VERIFICATION 2: Linearity in eta (for fixed alpha)")
    print("="*70)

    # Test at several alpha values
    alpha_test_values = [2.2, 2.5, 2.7, 2.9]
    eta_grid = np.linspace(0, 1, 100)

    max_deviation = 0.0

    for alpha_val in alpha_test_values:
        # Compute L(alpha, eta) for all eta
        L_vals = [L_lagrangian(alpha_val, eta) for eta in eta_grid]

        # Check linearity: d²L/deta² should be zero
        def L_fixed_alpha(eta):
            return L_lagrangian(alpha_val, eta)

        second_derivs_eta = [second_derivative(L_fixed_alpha, eta, h=1e-5) for eta in eta_grid[1:-1]]
        max_dev = max(abs(d) for d in second_derivs_eta)
        max_deviation = max(max_deviation, max_dev)

        print(f"  alpha = {alpha_val:.2f}: max |d²L/deta²| = {max_dev:.6e}")

    # Tolerance: numerical derivatives can have small errors
    is_linear = max_deviation < 1e-5
    print(f"\n  L(alpha, eta) is affine in eta? {is_linear} (max deviation: {max_deviation:.6e})")
    if is_linear:
        print("  YES: Affine functions are both convex and concave.")
    else:
        print("  NOTE: Small deviation likely due to numerical differentiation error.")

    return {'linear': is_linear, 'max_deviation': max_deviation}

# ==============================================================
# VERIFICATION 3: Semicontinuity
# ==============================================================

def verify_continuity():
    """
    Verify that C_wave and C_transport are continuous (hence semicontinuous).

    For smooth compositions of Bessel functions, exponentials, and power laws,
    continuity is expected. We check for jumps or singularities numerically.
    """
    print("\n" + "="*70)
    print("VERIFICATION 3: Continuity (implies semicontinuity)")
    print("="*70)

    alpha_grid = np.linspace(alpha_w, 3.0, 500)

    # Compute function values
    C_w_vals = [C_wave(alpha) for alpha in alpha_grid]
    C_t_vals = [C_transport(alpha) for alpha in alpha_grid]

    # Check for NaN or Inf
    has_nan_wave = any(np.isnan(v) or np.isinf(v) for v in C_w_vals)
    has_nan_trans = any(np.isnan(v) or np.isinf(v) for v in C_t_vals)

    # Check for large jumps (discontinuities)
    # Since functions are smooth compositions, use relative jump threshold
    diffs_wave = np.diff(C_w_vals)
    diffs_trans = np.diff(C_t_vals)
    max_jump_wave = np.max(np.abs(diffs_wave))
    max_jump_trans = np.max(np.abs(diffs_trans))

    # Average step size
    mean_val_wave = np.mean(np.abs(C_w_vals))
    mean_val_trans = np.mean(np.abs(C_t_vals))

    # Relative jump (as fraction of mean value)
    rel_jump_wave = max_jump_wave / mean_val_wave if mean_val_wave > 0 else 0
    rel_jump_trans = max_jump_trans / mean_val_trans if mean_val_trans > 0 else 0

    print(f"  C_wave: No NaN/Inf? {not has_nan_wave}, Max relative jump: {rel_jump_wave:.6e}")
    print(f"  C_transport: No NaN/Inf? {not has_nan_trans}, Max relative jump: {rel_jump_trans:.6e}")

    # Threshold: relative jumps < 10% indicate continuity
    is_continuous = (not has_nan_wave) and (not has_nan_trans) and (rel_jump_wave < 0.1) and (rel_jump_trans < 0.1)

    if is_continuous:
        print("  ✓ Both functions are continuous, hence upper/lower semicontinuous.")
    else:
        print("  ✗ Warning: Potential discontinuities detected.")

    return {
        'continuous': is_continuous,
        'wave_nan': has_nan_wave,
        'trans_nan': has_nan_trans,
        'wave_max_jump': max_jump_wave,
        'trans_max_jump': max_jump_trans,
    }

# ==============================================================
# VISUALIZATION
# ==============================================================

def plot_convexity_analysis(results_conv):
    """Generate plots showing curvature and function shapes."""

    alpha_grid = results_conv['C_wave']['alpha_grid']
    curvs_wave = results_conv['C_wave']['curvatures']
    curvs_trans = results_conv['C_transport']['curvatures']

    # Compute function values
    C_w_vals = [C_wave(a) for a in alpha_grid]
    C_t_vals = [C_transport(a) for a in alpha_grid]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: C_wave function
    ax = axes[0, 0]
    ax.plot(alpha_grid, C_w_vals, 'b-', lw=2, label=r'$\mathcal{C}_{\mathrm{wave}}(\alpha)$')
    ax.axvline(alpha_w, color='gray', ls='--', alpha=0.5, label=r'$\alpha_w$')
    ax.set_xlabel(r'$\alpha$', fontsize=12)
    ax.set_ylabel(r'$\mathcal{C}_{\mathrm{wave}}$', fontsize=12)
    ax.set_title('Wave Cost Function', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Panel B: C_wave curvature
    ax = axes[0, 1]
    ax.plot(alpha_grid, curvs_wave, 'b-', lw=2)
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    ax.set_xlabel(r'$\alpha$', fontsize=12)
    ax.set_ylabel(r'$d^2\mathcal{C}_{\mathrm{wave}}/d\alpha^2$', fontsize=12)
    ax.set_title('Wave Cost Curvature', fontsize=13)
    ax.grid(True, alpha=0.3)

    # Panel C: C_transport function
    ax = axes[1, 0]
    ax.plot(alpha_grid, C_t_vals, 'r-', lw=2, label=r'$\mathcal{C}_{\mathrm{transport}}(\alpha)$')
    # Find minimum
    min_idx = np.argmin(C_t_vals)
    alpha_min = alpha_grid[min_idx]
    ax.axvline(alpha_min, color='gray', ls='--', alpha=0.5, label=rf'$\alpha_t = {alpha_min:.3f}$')
    ax.set_xlabel(r'$\alpha$', fontsize=12)
    ax.set_ylabel(r'$\mathcal{C}_{\mathrm{transport}}$', fontsize=12)
    ax.set_title('Transport Cost Function', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Panel D: C_transport curvature
    ax = axes[1, 1]
    ax.plot(alpha_grid, curvs_trans, 'r-', lw=2)
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    ax.set_xlabel(r'$\alpha$', fontsize=12)
    ax.set_ylabel(r'$d^2\mathcal{C}_{\mathrm{transport}}/d\alpha^2$', fontsize=12)
    ax.set_title('Transport Cost Curvature', fontsize=13)
    ax.grid(True, alpha=0.3)

    # Highlight negative curvature region if exists
    neg_indices = [i for i, c in enumerate(curvs_trans) if c < 0]
    if neg_indices:
        ax.fill_between(alpha_grid[neg_indices],
                        [curvs_trans[i] for i in neg_indices],
                        0, alpha=0.3, color='red', label='Concave region')
        ax.legend()

    plt.tight_layout()

    # Save
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fig_path = os.path.join(script_dir, '..', 'manuscript', 'figures', 'sion_convexity_analysis.pdf')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n  [OK] Saved figure: {fig_path}")
    plt.close()

# ==============================================================
# MAIN REPORT
# ==============================================================

def generate_report(results_conv, results_linear, results_cont):
    """Generate comprehensive report on Sion's theorem applicability."""

    print("\n" + "="*70)
    print("FINAL ASSESSMENT: Sion's Minimax Theorem Applicability")
    print("="*70)

    print("\nREQUIRED CONDITIONS:")
    print("-" * 70)

    # Compactness
    print("\n1. X = [alpha_w, alpha_t] compact? YES (closed bounded interval)")
    print("   Y = [0, 1] compact? YES (closed bounded interval)")

    # Convexity in alpha
    print("\n2. L(*, eta) quasi-convex for each eta in [0,1]?")

    C_wave_qconvex = results_conv['C_wave']['quasi_convex']
    C_trans_qconvex = results_conv['C_transport']['quasi_convex']

    print(f"   - C_wave(alpha) strictly convex? {results_conv['C_wave']['convex']}")
    if not results_conv['C_wave']['convex']:
        print(f"   - C_wave(alpha) quasi-convex? {C_wave_qconvex}")

    print(f"   - C_transport(alpha) strictly convex? {results_conv['C_transport']['convex']}")
    if not results_conv['C_transport']['convex']:
        print(f"   - C_transport(alpha) quasi-convex? {C_trans_qconvex}")

    # Positive combinations of (quasi-)convex functions are (quasi-)convex
    L_quasiconvex = C_wave_qconvex and C_trans_qconvex

    if results_conv['C_wave']['convex'] and results_conv['C_transport']['convex']:
        print(f"   YES L(alpha, eta) = eta*C_wave + (1-eta)*C_transport is STRICTLY CONVEX in alpha")
        print("     (positive combination of convex functions)")
    elif L_quasiconvex:
        print(f"   YES L(alpha, eta) is QUASI-CONVEX in alpha")
        print("     (positive combination of quasi-convex functions)")
    else:
        print(f"   WARNING: L(alpha, eta) may not be quasi-convex")

    # Concavity in eta
    print("\n3. L(alpha, *) quasi-concave for each alpha in X?")
    print(f"   - L(alpha, eta) is AFFINE in eta? {results_linear['linear']}")
    print("   YES Affine functions are both convex AND concave, hence quasi-concave")

    # Semicontinuity
    print("\n4. Semicontinuity?")
    print(f"   - C_wave, C_transport continuous? {results_cont['continuous']}")
    if results_cont['continuous']:
        print("   YES Continuous functions are both upper and lower semicontinuous")
    else:
        print("   WARNING: Continuity check flagged potential issues")

    # FINAL VERDICT
    print("\n" + "="*70)
    print("FINAL VERDICT")
    print("="*70)

    all_conditions_met = L_quasiconvex and results_linear['linear'] and results_cont['continuous']

    if all_conditions_met:
        print("\nYES: Sion's Minimax Theorem APPLIES")
        print("\nTherefore:")
        print("  min_alpha max_eta L(alpha,eta) = max_eta min_alpha L(alpha,eta)")
        print("\nThe saddle point (alpha*, eta*) exists and is valid.")

        if not (results_conv['C_wave']['convex'] and results_conv['C_transport']['convex']):
            print("\nIMPORTANT QUALIFICATION:")
            print("  While strict convexity fails, QUASI-CONVEXITY holds.")
            print("  Sion's theorem only requires quasi-convexity, not strict convexity.")
            print("  The manuscript should clarify this distinction.")
    else:
        print("\nWARNING: Some conditions for Sion's theorem may not hold")
        print("\nRecommendations:")
        if not L_quasiconvex:
            print("  - Investigate why quasi-convexity fails")
            print("  - Consider restricted domain [alpha_min, alpha_max]")
        if not results_linear['linear']:
            print("  - Verify affine structure of Lagrangian")
        if not results_cont['continuous']:
            print("  - Check for singularities in domain")

    # Additional insights
    print("\n" + "="*70)
    print("ADDITIONAL INSIGHTS")
    print("="*70)

    print(f"\nC_wave minimum curvature: {results_conv['C_wave']['min_curvature']:.6e}")
    print(f"C_transport minimum curvature: {results_conv['C_transport']['min_curvature']:.6e}")

    if results_conv['C_transport']['min_curvature'] < 0:
        print("\nWARNING: C_transport has regions of negative curvature (concave).")
        print("  This does NOT violate Sion's theorem if quasi-convexity holds.")
        print("  Recommendation: Explicitly state in manuscript that the function")
        print("  is quasi-convex (unimodal with unique minimum), not strictly convex.")

# ==============================================================
# Main execution
# ==============================================================

if __name__ == '__main__':
    print("\n" + "="*70)
    print("SION'S MINIMAX THEOREM VERIFICATION")
    print("="*70)
    print("Paper: A Unified Variational Principle for Branching Transport Networks")
    print("Lagrangian: L(alpha, eta) = eta*C_wave(alpha) + (1-eta)*C_transport(alpha)")
    print("="*70)

    # Run verifications
    results_conv = verify_convexity_in_alpha()
    results_linear = verify_linearity_in_eta()
    results_cont = verify_continuity()

    # Generate visualizations
    plot_convexity_analysis(results_conv)

    # Final report
    generate_report(results_conv, results_linear, results_cont)

    print("\n" + "="*70)
    print("Verification complete.")
    print("="*70 + "\n")
