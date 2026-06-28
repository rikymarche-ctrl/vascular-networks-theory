r"""
verify_heterogeneity.py

This script performs a first-order sensitivity analysis to evaluate the robustness
of the global Minimax attractor—which is theoretically anchored by the Kinematic
Matching Criterion (\mathrm{Wo}_c = \sqrt{3})—under anatomical and morphometric
heterogeneities.

Specifically, it demonstrates how asymmetry and vessel tapering act as opposing
remodeling forces on the wave-reflection landscape:
  - Morphometric asymmetry (A < 1.0) perturbs the wave-reflection landscape and
    shifts the emergent branching exponent relative to the symmetric ground state.
  - Distal vessel tapering (taper < 1.0) acts as a stabilizing mechanism that
    keeps the network close to the wave-optimum.

All computations delegate directly to compute.py, which uses the rigorous
impedance-based wave reflection formulation with full harmonic spectrum summation
(n = 1..5 Fourier components of the physiological pulse).
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'main'))

# Import canonical minimax solvers from compute.py (multi-harmonic, physiological).
# These are the same functions used to generate all paper/supplement values.
from compute import (  # type: ignore # pylint: disable=import-error
    find_minimax_hetero,
    find_minimax_hetero_elastic,
)

# ---------------------------------------------------------------------------
# NOTE: No local redefinitions of C_wave_hetero, C_transport, or
# find_minimax_hetero.  All computations use the canonical engine in
# compute.py, which sums over PULSE_HARMONICS (n = 1..5) for physiological
# accuracy.  Any single-harmonic approximation would give different (incorrect)
# numerical results and is NOT representative of the paper's claims.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 70)
    print("STRUCTURAL HETEROGENEITY ANALYSIS (Impedance Formulation)")
    print("=" * 70)
    print("\nTesting asymmetry and taper effects on minimax alpha*")
    print("Parameters: M=70 kg, G=11, N=2")
    print("Wave cost: full harmonic spectrum (n=1..5) from compute.py\n")

    print("[1] BASELINE (Symmetric, No Taper)")
    alpha_base = find_minimax_hetero(70.0, A=1.0, taper=1.0)
    print(f"    alpha* = {alpha_base:.4f}")

    print("\n[2] ASYMMETRY (A = 0.85)")
    alpha_asym = find_minimax_hetero(70.0, A=0.85, taper=1.0)
    print(f"    alpha* = {alpha_asym:.4f}")
    delta_asym = alpha_asym - alpha_base
    print(f"    Shift: Delta_alpha = {delta_asym:+.4f}")

    print("\n[3] TAPER (2% per generation)")
    alpha_taper = find_minimax_hetero(70.0, A=1.0, taper=0.98)
    print(f"    alpha* = {alpha_taper:.4f}")
    delta_taper = alpha_taper - alpha_base
    print(f"    Shift: Delta_alpha = {delta_taper:+.4f}")

    print("\n[4] COMBINED (A=0.85 + 2% Taper)")
    alpha_full = find_minimax_hetero(70.0, A=0.85, taper=0.98)
    print(f"    alpha* = {alpha_full:.4f}")
    delta_full = alpha_full - alpha_base
    print(f"    Net Shift: Delta_alpha = {delta_full:+.4f}")

    print("\n[5] ELASTIC WALL CORRECTION (baseline)")
    alpha_elastic_base = find_minimax_hetero_elastic(70.0, A=1.0, taper=1.0, alpha_w=2.115)
    print(f"    alpha* = {alpha_elastic_base:.4f}")

    print("\n[6] ELASTIC + COMBINED HETEROGENEITY")
    alpha_elastic_full = find_minimax_hetero_elastic(70.0, A=0.85, taper=0.98, alpha_w=2.115)
    print(f"    alpha* = {alpha_elastic_full:.4f}")
    delta_elastic = alpha_elastic_full - alpha_elastic_base
    print(f"    Net Shift: Delta_alpha = {delta_elastic:+.4f}")

    print("\n" + "=" * 70)
    print("SUMMARY (canonical multi-harmonic values from compute.py)")
    print("=" * 70)
    print(f"Baseline:              alpha* = {alpha_base:.4f}")
    print(f"Asymmetry shift:       Delta  = {delta_asym:+.4f}")
    print(f"Taper shift:           Delta  = {delta_taper:+.4f}")
    print(f"Combined:              alpha* = {alpha_full:.4f}")
    print(f"Elastic (baseline):    alpha* = {alpha_elastic_base:.4f}")
    print(f"Elastic + Combined:    alpha* = {alpha_elastic_full:.4f}")
    print(f"\nEmpirical (Kassab):    alpha* ~ 2.72")
    print(f"Remaining gap after first-order rigid correction: {2.72 - alpha_full:.4f}")
    print("=" * 70)
