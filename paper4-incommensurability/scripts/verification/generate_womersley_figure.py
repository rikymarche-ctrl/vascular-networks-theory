#!/usr/bin/env python3
"""
Generate Womersley Verification Figure for Supplemental Material

This script numerically verifies the critical Womersley number prediction
Wo_c = sqrt(3) for 3D branching networks by solving the phase constraint
tan(phi) = d - 1, where phi is the reactive phase lag.

Output:
    - supplements/figures/womersley_verification.png
    - Console output with numerical verification statistics
"""

import numpy as np
from scipy.special import jv
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import os


def F10(Wo):
    """
    Womersley function F10(z) = 2*J1(z)/(z*J0(z)).

    Args:
        Wo: Womersley number

    Returns:
        Complex F10 value with z = Wo * exp(i*3*pi/4)
    """
    if Wo < 1e-4:
        return 1.0 - (1j * Wo**2 / 8.0)
    z = Wo * np.exp(1j * 0.75 * np.pi)
    j0 = jv(0, z)
    j1 = jv(1, z)
    return (2 * j1) / (z * j0)


def phase_function(Wo):
    """
    Computes the reactive phase lag phi of the admittance factor.
    We negate the angle since the 3pi/4 convention yields conjugate values.

    Args:
        Wo: Womersley number

    Returns:
        Phase lag in radians
    """
    return -np.angle(1 / (1 - F10(Wo)))


def solve_critical_womersley(d=3):
    """
    Solve Q_inv = d - 1 for critical Womersley number.

    Args:
        d: Network dimension (default: 3 for 3D networks)

    Returns:
        Critical Womersley number Wo_c
    """
    target = d - 1
    objective = lambda w: np.tan(phase_function(w)) - target

    # The function is monotonic in the region of interest
    root = brentq(objective, 1.0, 3.0)
    return root


def generate_verification_plot(output_path="../../supplements/figures/womersley_verification.png"):
    """
    Generate and save the Womersley verification figure.

    Args:
        output_path: Path where to save the figure
    """
    # Numerical solution
    Wo_c_num = solve_critical_womersley(d=3)
    Wo_c_analytic = np.sqrt(3)

    # Print verification statistics
    print("=" * 60)
    print("WOMERSLEY NUMBER VERIFICATION (d=3)")
    print("=" * 60)
    print(f"Numerical Critical Womersley:  {Wo_c_num:.6f}")
    print(f"Analytical Prediction (sqrt3): {Wo_c_analytic:.6f}")
    print(f"Absolute Error:                {abs(Wo_c_num - Wo_c_analytic):.2e}")
    print(f"Relative Error:                {abs(Wo_c_num - Wo_c_analytic)/Wo_c_analytic:.2e}")
    print("=" * 60)

    # Generate plot data
    Wo_vals = np.linspace(1e-3, 5.0, 500)
    Q_inv = np.array([np.tan(phase_function(w)) for w in Wo_vals])

    # Create figure
    plt.figure(figsize=(8, 5))
    plt.plot(Wo_vals, Q_inv, 'b-', linewidth=2, label=r"$\mathcal{Q}^{-1}(\mathrm{Wo})$")
    plt.axhline(y=2, color='r', linestyle='--', linewidth=1.5, label=r"$\mathcal{Q}^{-1} = d-1=2$")
    plt.axvline(x=np.sqrt(3), color='g', linestyle=':', linewidth=1.5, label=r"$\mathrm{Wo}_c = \sqrt{3}$")

    # Mark intersection point
    plt.plot(Wo_c_num, 2, 'ko', markersize=8, zorder=5)

    plt.xlabel(r"Womersley Number $\mathrm{Wo}$", fontsize=12)
    plt.ylabel(r"$\mathcal{Q}^{-1}$", fontsize=12)
    plt.title("Numerical Verification of Critical Womersley Threshold", fontsize=13, fontweight='bold')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.ylim(0, 10)
    plt.xlim(0, 5)
    plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=10, frameon=True)
    plt.subplots_adjust(bottom=0.18)

    # Save figure: 600 dpi raster (editor requirement for linework) plus a
    # resolution-independent vector PDF for press-quality reproduction.
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    pdf_path = os.path.splitext(output_path)[0] + ".pdf"
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"\n[OK] Verification plot saved to: {output_path} (600 dpi)")
    print(f"[OK] Vector PDF saved to: {pdf_path}")
    print()


if __name__ == "__main__":
    generate_verification_plot()
