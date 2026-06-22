"""
Numerical verification of the WAVE critical Womersley threshold against the
exact Womersley/Bessel solution (TODO item A1).

Background
----------
The fluid threshold is the root of  Q^{-1}_{Y_L}(Wo) = d-1  and is already
verified elsewhere (Wo_c^fluid = sqrt(6/(d-1)); d=3 -> sqrt(3) ~ 1.732, exact
root ~ 1.740).

The WAVE threshold concerns the *characteristic* admittance Y_c. Two readings
were debated:
  (A) Y_c proportional to sqrt(1 - F10)            [Womersley/McDonald]
  (B) Y_c proportional to sqrt(Y_L) = sqrt(-i(1-F10))   [naive "sqrt(Y_L)"]
They differ by the factor sqrt(-i) = e^{-i pi/4}. This script computes the
inverse quality factor  Q^{-1} = Re(.)/|Im(.)|  for BOTH readings from the
*exact* Bessel functions and finds the root of  Q^{-1} = d-1, settling whether
the paper's quoted Wo_c^wave = 3/sqrt(2) ~ 2.121 is correct (reading A) or
whether it should be ~2.83 (reading B). It also checks the d=2 behaviour.

F10(z) = 2 J1(z) / (z J0(z)),  z = Wo * exp(-i pi/4)   (i^{3/2} convention,
consistent with the e^{i w t} time convention used in the manuscript).
"""

import cmath

try:
    import mpmath as mp
    mp.mp.dps = 30

    def besselj(n, z):
        return mp.besselj(n, z)

    def to_complex(z):
        return complex(z)
    BACKEND = "mpmath"
except Exception:  # pragma: no cover - fallback
    from scipy.special import jv as _jv

    def besselj(n, z):
        return _jv(n, z)

    def to_complex(z):
        return complex(z)
    BACKEND = "scipy"


def F10(Wo):
    """Exact Womersley function 2 J1(z)/(z J0(z)) with z = Wo e^{-i pi/4}."""
    z = Wo * cmath.exp(-1j * cmath.pi / 4)
    j0 = to_complex(besselj(0, z))
    j1 = to_complex(besselj(1, z))
    return 2.0 * j1 / (z * j0)


def qinv_YL(Wo):
    """Q^{-1} of the longitudinal admittance Y_L ∝ (1-F10)/(i) = -i(1-F10)."""
    g = 1.0 - F10(Wo)
    YL = -1j * g
    return abs(YL.real) / abs(YL.imag)


def qinv_Yc_readingA(Wo):
    """Reading A: Y_c ∝ sqrt(1 - F10)  (Womersley/McDonald)."""
    g = 1.0 - F10(Wo)
    Yc = cmath.sqrt(g)
    return abs(Yc.real) / abs(Yc.imag)


def qinv_Yc_readingB(Wo):
    """Reading B: Y_c ∝ sqrt(Y_L) = sqrt(-i(1-F10))."""
    g = 1.0 - F10(Wo)
    Yc = cmath.sqrt(-1j * g)
    return abs(Yc.real) / abs(Yc.imag)


def bisect_root(f, target, lo, hi, tol=1e-9, itmax=200):
    """Find Wo in [lo,hi] with f(Wo)=target, assuming f is monotone there."""
    flo = f(lo) - target
    fhi = f(hi) - target
    if flo * fhi > 0:
        return None  # no sign change -> no root in bracket
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        fm = f(mid) - target
        if abs(fm) < tol:
            return mid
        if flo * fm < 0:
            hi = mid
            fhi = fm
        else:
            lo = mid
            flo = fm
    return 0.5 * (lo + hi)


def report():
    print(f"[backend: {BACKEND}]")
    print("=" * 64)
    print("FLUID threshold  (Q^-1_{Y_L} = d-1)")
    for d, name in [(3, "3D"), (2, "2D")]:
        root = bisect_root(qinv_YL, d - 1, 0.3, 6.0)
        analytic = (6.0 / (d - 1)) ** 0.5
        print(f"  d={d} ({name}): exact root Wo_c^fluid = "
              f"{root:.4f}   (analytic sqrt(6/(d-1)) = {analytic:.4f})")
    print("-" * 64)
    print("WAVE threshold   (Q^-1_{Y_c} = d-1)")
    print("  Reading A:  Y_c ~ sqrt(1 - F10)   [paper / Womersley-McDonald]")
    for d in (3, 2):
        # scan to see monotonicity / sign change
        root = bisect_root(qinv_Yc_readingA, d - 1, 0.05, 12.0)
        val_lo = qinv_Yc_readingA(0.2)
        val_hi = qinv_Yc_readingA(10.0)
        msg = f"{root:.4f}" if root is not None else "no finite root in (0.05,12]"
        print(f"    d={d}: Wo_c^wave(A) = {msg}    "
              f"[Q^-1 at Wo=0.2 -> {val_lo:.3f}, at Wo=10 -> {val_hi:.3f}]")
    print("  paper value 3/sqrt(2) = 2.1213")
    print("  Reading B:  Y_c ~ sqrt(Y_L) = sqrt(-i(1-F10))   [naive]")
    for d in (3, 2):
        root = bisect_root(qinv_Yc_readingB, d - 1, 0.05, 12.0)
        msg = f"{root:.4f}" if root is not None else "no finite root in (0.05,12]"
        print(f"    d={d}: Wo_c^wave(B) = {msg}")
    print("=" * 64)


if __name__ == "__main__":
    report()
