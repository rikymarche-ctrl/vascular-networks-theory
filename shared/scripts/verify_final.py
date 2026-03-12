"""Final verification of Papers 1 & 2."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import compute_paper2 as p2
import compute_paper1 as p1

# Paper 1 (Murray)
r = p2.compute_all()
print("=== Paper 1 (Murray) ===")
print(f"  alpha*(mw=5)  = {r['alpha_mw_low']:.3f}")
print(f"  alpha*(mw=35) = {r['alpha_mw_high']:.3f}")
print(f"  delta_alpha   = {abs(r['delta_alpha_ours']):.2f}")

# Paper 2 (Variational)
eta = p1.optimal_eta()
al_mean, _, _ = p1.alpha_local_mean()
alpha_star, alpha_t = p1.find_alpha_star(eta)
print()
print("=== Paper 2 (Variational) ===")
print(f"  eta*       = {eta:.4f}")
print(f"  alpha_t    = {alpha_t:.3f}")
print(f"  alpha*     = {alpha_star:.3f}")
print(f"  error      = {abs(alpha_star - 2.7)/2.7 * 100:.1f}%")

# Sensitivity check (η)
print()
print("=== Sensitivity (eta) ===")
all_ok = True
for test_eta in [0.6, 0.7, 0.73, 0.8, 0.9]:
    a, _ = p1.find_alpha_star(test_eta)
    print(f"  eta={test_eta:.2f} -> alpha*={a:.3f}")

# Cross-checks
print()
print("=== Cross-Checks ===")
checks = []

# 1. alpha_local matches Paper 1 alpha* range
ok = 2.885 < al_mean < r['alpha_mw_low'] + 0.05
checks.append(("alpha_local vs Paper 1", ok, f"{al_mean:.3f} in [{2.885:.3f}, {r['alpha_mw_low']+0.05:.3f}]"))

# 2. alpha* in experimental range
ok = abs(alpha_star - 2.7) < 0.2
checks.append(("alpha* vs experiment", ok, f"|{alpha_star:.3f} - 2.70| = {abs(alpha_star-2.7):.3f} < 0.20"))

# 3. delta_alpha Murray > delta_alpha ours
ok = abs(r['delta_alpha_murray']) > abs(r['delta_alpha_ours'])
checks.append(("Murray gap > our gap", ok, f"{abs(r['delta_alpha_murray']):.2f} > {abs(r['delta_alpha_ours']):.2f}"))

# 4. alpha_w = (5-p)/2 (from params)
import params
ok = abs(params.alpha_w - 2.115) < 0.001
checks.append(("alpha_w = (5-p)/2", ok, f"{params.alpha_w:.3f} == 2.115"))

# 5. Bound: (5+p)/2 < alpha_local < 3
bound_lower = (5 + 0.77) / 2
ok = bound_lower < al_mean < 3.0
checks.append(("(5+p)/2 < alpha_t < 3", ok, f"{bound_lower:.3f} < {al_mean:.3f} < 3.000"))

for name, ok, detail in checks:
    status = "[PASS]" if ok else "[FAIL]"
    print(f"  {status} {name}: {detail}")
    if not ok:
        all_ok = False

print()
if all_ok:
    print("ALL CHECKS PASSED")
else:
    print("SOME CHECKS FAILED")
    sys.exit(1)
