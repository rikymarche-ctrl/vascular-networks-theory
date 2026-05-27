r"""
verify_heterogeneity.py

This script performs a first-order sensitivity analysis to evaluate the robustness
of the global Minimax attractor—which is theoretically anchored by the Kinematic
Matching Criterion (\mathrm{Wo}_c = \sqrt{3})—under anatomical and morphometric
heterogeneities.

Specifically, it demonstrates how asymmetry and vessel tapering act as opposing
remodeling forces:
  - Morphometric asymmetry (A < 1.0) pulls the emergent branching exponent upward
    relative to the symmetric ground state.
  - Distal vessel tapering (taper < 1.0) pulls the exponent downward, acting as a
    vital stabilizing mechanism that keeps the network close to the wave-optimum.

Uses the rigorous impedance-based wave reflection formulation from compute.py.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'main'))

import numpy as np
from scipy.optimize import brentq
from scipy.special import jv
import math

# Import base parameters and functions (resolved dynamically at runtime via sys.path)
from compute import (  # type: ignore # pylint: disable=import-error
    rho, mu_f, nu, f_h0, M0, r0_human, N, G, Wo_c, p,
    b_blood, m_wall, phi0, r_opt, admittance_factor, propagation_const, lambda_L
)

def C_transport(alpha, Q0, ell0, beta_l):
    """Compute the fractional excess metabolic transport cost for a heterogeneous network.

    This function calculates the additional metabolic power required by the network
    when the branching exponent deviates from the Murray law optimum (α = 3). It
    integrates the local viscous dissipation, blood volume, and wall material costs
    across all generations.

    Args:
        alpha: Branching exponent (area-preserving law parameter)
        Q0: Volumetric flow rate at the root vessel (m³/s)
        ell0: Reference length of the root vessel (m)
        beta_l: Length scaling factor per generation

    Returns:
        Fractional excess transport cost C_T ∈ [0, ∞)
    """
    r0_s = r_opt(Q0)
    num, den = 0.0, 0.0
    for g in range(G):
        Qg = Q0 / N**g
        rs = r_opt(Qg)
        rg = r0_s * (N**(-g/alpha))
        wg = N**g * ell0 * beta_l**g
        num += wg * (phi0(rg, Qg) - phi0(rs, Qg))
        den += wg * phi0(rs, Qg)
    return num / den

def C_wave_hetero(alpha, M, A=1.0, taper=1.0):
    """Compute the cumulative wave reflection cost accounting for structural asymmetry and taper.

    Propagates pressure wave energy through a structurally heterogeneous bifurcation tree,
    evaluating both the attenuation along each parent segment and the reflections generated
    by junction impedance mismatch and continuous vessel tapering.

    Args:
        alpha: Branching exponent (area-preserving law parameter)
        M: Body mass of the organism (kg)
        A: Daughter vessel asymmetry ratio r₂/r₁ (default: 1.0, symmetric)
        taper: Radial taper ratio per generation (default: 1.0, no taper)

    Returns:
        Total reflected power (dimensionless cost) normalized to incident power
    """
    f_h = f_h0 * (M / M0)**(-0.25)
    omega = 2 * np.pi * f_h
    r_root = r0_human * (M / M0)**0.375

    P_at_node = 1.0
    reflected_power = 0.0
    curr_r = r_root
    
    for g in range(G):
        Wo_p = curr_r * math.sqrt(omega / nu)
        
        # Rigorous wave propagation attenuation along the parent segment
        kappa = np.real(propagation_const(Wo_p, omega))
        L_parent = lambda_L * curr_r
        P_at_node *= np.exp(-2.0 * kappa * L_parent)
        
        # Taper and asymmetry radii
        r_p_distal = curr_r * taper
        r1 = r_p_distal * (1.0 + A**alpha)**(-1.0 / alpha)
        r2 = A * r1
        
        # Admittances at the junction
        Yp = (r_p_distal**2) * admittance_factor(Wo_p * taper)
        Y1 = (r1**2) * admittance_factor(math.sqrt(omega / nu) * r1)
        Y2 = (r2**2) * admittance_factor(math.sqrt(omega / nu) * r2)
        
        # Junction reflection
        Gamma_junction = np.abs((Yp - (Y1 + Y2)) / (Yp + (Y1 + Y2)))**2
        
        # Taper reflection
        Y_prox = (curr_r**2) * admittance_factor(Wo_p)
        Y_dist = (r_p_distal**2) * admittance_factor(Wo_p * taper)
        Gamma_taper = np.abs((Y_prox - Y_dist) / (Y_prox + Y_dist))**2
        
        # Series reflector (incoherent)
        ref = Gamma_junction + (1.0 - Gamma_junction) * Gamma_taper
        
        reflected_power += P_at_node * ref
        P_at_node *= (1.0 - ref)
        
        # Step to child branch
        curr_r = r1
        
    return reflected_power

def find_minimax_hetero(M, A=1.0, taper=1.0):
    """Find the emergent minimax branching exponent alpha* for the given heterogeneous parameters.

    Solves the zero-sum balance equation residual(α) = C_wave(α) - C_transport(α) = 0
    under Womersley pulsatile flow to determine the robust minimax attractor.

    Args:
        M: Body mass of the organism (kg)
        A: Daughter vessel asymmetry ratio r₂/r₁ (default: 1.0, symmetric)
        taper: Radial taper ratio per generation (default: 1.0, no taper)

    Returns:
        Minimax branching exponent alpha* (or np.nan if no root is found)
    """
    Q0 = 1.3e-6 * (M / M0)**0.75
    ell0 = 0.015 * (M / M0)**(1.0/3.0)
    beta_l = 0.787

    def residual(a):
        return C_wave_hetero(a, M, A, taper) - C_transport(a, Q0, ell0, beta_l)

    try:
        return brentq(residual, 2.1, 2.95, xtol=1e-5)
    except ValueError as e:
        print(f"   Warning: No root for A={A}, taper={taper}: {e}")
        return np.nan

if __name__ == "__main__":
    print("=" * 70)
    print("STRUCTURAL HETEROGENEITY ANALYSIS (Impedance Formulation)")
    print("=" * 70)
    print("\nTesting asymmetry and taper effects on minimax alpha*")
    print("Parameters: M=70kg, G=11, N=2\n")

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

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Baseline:              alpha* = {alpha_base:.4f}")
    print(f"Asymmetry shift:       Delta = {delta_asym:+.4f}")
    print(f"Taper shift:           Delta = {delta_taper:+.4f}")
    print(f"Combined:              alpha* = {alpha_full:.4f}")
    print(f"\nEmpirical (Kassab):    alpha* ~ 2.72")
    print(f"Residual gap:          {2.72 - alpha_full:.4f}")
    print("=" * 70)
