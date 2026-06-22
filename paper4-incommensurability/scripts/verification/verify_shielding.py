r"""
verify_shielding.py
Topological Shielding Sensitivity Analysis — Paper IV Supplemental

This script verifies the biophysical logic of "Topological Shielding": the global
minimax attractor alpha* is structurally decoupled from the distally localized
thick-wall (Lamé) elastic correction in the terminal arterioles.

It is a thin presentation layer over the canonical engine in compute.py. The
saddle alpha* is found with the SAME eta-free marginal-balance criterion as the
main model (compute.find_minimax solves C_wave = C_visc): here we solve
C_wave_lame(alpha, M) = C_visc(alpha, M), with and without the thick-wall
correction. No fixed Lagrangian weight eta is imposed, so the reported shift
|Delta alpha*| is a pure consequence of the Lamé correction. The numbers printed
here are exactly those exported to manuscript/dynamic_variables.tex as the
\VarShield* macros and shown in Table~\ref{tab:shielding}.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'main'))
import numpy as np
import math

# Import the canonical engine (resolved dynamically at runtime via sys.path)
from compute import (  # type: ignore # pylint: disable=import-error
    nu, f_h0, r0_human, M0, N, G, A_ratio,
    admittance_factor, lame_factor, find_minimax_lame,
)

WOC = np.sqrt(3)   # critical Womersley number


def _wo0_of_M(M):
    """Fundamental aortic-root Womersley number Wo_0 (generation 0) at mass M."""
    f_h = f_h0 * (M / M0) ** (-0.25)
    r_root = r0_human * (M / M0) ** 0.375
    return r_root * math.sqrt(2 * np.pi * f_h / nu)


def _mass_for_wo0(wo0_target):
    """Body mass whose aortic-root Womersley number equals wo0_target."""
    return M0 * (wo0_target / _wo0_of_M(M0)) ** 4


def gamma_sq_junction(alpha, Wo, L_parent=1.0, L_child=1.0):
    """Squared junction reflection coefficient with optional thick-wall factors.

    Uses the same baseline asymmetry A_ratio as the canonical engine. Setting the
    Lamé factors L_parent = L_child = 1 recovers the thin-wall reflection.
    """
    r_scale = (1.0 + A_ratio ** alpha) ** (-1.0 / alpha)
    Yp = admittance_factor(Wo) / L_parent
    Y1 = (r_scale ** 2) * admittance_factor(Wo * r_scale) / L_child
    Y2 = ((A_ratio * r_scale) ** 2) * admittance_factor(Wo * A_ratio * r_scale) / L_child
    return np.abs((Yp - (Y1 + Y2)) / (Yp + (Y1 + Y2))) ** 2


# ── Main analysis ───────────────────────────────────────────────────────────

def run():
    print("=" * 64)
    print("TOPOLOGICAL SHIELDING - ETA-FREE SADDLE SENSITIVITY ANALYSIS")
    print("Saddle criterion: C_wave_lame(alpha, M) = C_visc(alpha, M)  [no eta]")
    print("Threshold: |Delta Alpha*| < 0.01")
    print("=" * 64)

    # Lame sanity check
    print("\nLame factor check (should be approx 1.02 and approx 1.10):")
    print(f"  h/r = 0.08 -> L = {lame_factor(0.08):.4f}")
    print(f"  h/r = 0.42 -> L = {lame_factor(0.42):.4f}")

    # Sensitivity across the distal physiological Womersley range
    print("\n{:<8} {:<10} {:<14} {:<14} {:<14} {}".format(
        "Wo_0", "M (g)", "Alpha* (thin)", "Alpha* (Lame)", "|Delta Alpha*|", "Result"))
    print("-" * 64)
    all_pass = True
    for wo_0 in [2.0, 4.0, 8.0]:
        M_t = _mass_for_wo0(wo_0)
        a_thin = find_minimax_lame(M_t, use_lame=False)
        a_lame = find_minimax_lame(M_t, use_lame=True)
        diff = abs(a_thin - a_lame)
        ok = diff < 0.01
        all_pass = all_pass and ok
        print(f"{wo_0:<8.1f} {M_t*1000:<10.2f} {a_thin:<14.4f} {a_lame:<14.4f} "
              f"{diff:<14.6f} {'OK' if ok else 'FAIL'}")

    print("-" * 64)
    print(f"Overall: {'SHIELDING CONFIRMED' if all_pass else 'SHIELDING FAILED'}")

    # Generation-by-generation reflection breakdown at Wo_0 = 4
    print("\nGeneration breakdown at Wo_0=4.0, alpha=2.72 (thin vs Lame):")
    print("{:<6} {:<8} {:<8} {:<14} {:<14}".format(
        "gen", "Wo", "h/r", "gamma^2 (thin)", "gamma^2 (Lame)"))
    gens = np.arange(G + 1)
    wo_g = 4.0 * N ** (-gens / 2.72)
    h_r = 0.08 + 0.34 * gens / G
    L = np.array([lame_factor(x) for x in h_r])
    for g in range(G):
        gs_t = gamma_sq_junction(2.72, wo_g[g], 1.0, 1.0)
        gs_l = gamma_sq_junction(2.72, wo_g[g], L[g], L[g + 1])
        print(f"{g:<6d} {wo_g[g]:<8.3f} {h_r[g]:<8.3f} {gs_t:<14.2e} {gs_l:<14.2e}")


if __name__ == "__main__":
    run()
