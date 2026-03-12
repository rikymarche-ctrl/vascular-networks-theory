"""
run_all.py -- Master script for the branching-networks monorepo.

Regenerates all dynamic variables, figures, and runs verification checks.
Run from the repository root:  python shared/scripts/run_all.py
"""

import os
import sys
import importlib

# Ensure we can import from shared/scripts
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
sys.path.insert(0, SCRIPT_DIR)

# Output paths
P1_MANUSCRIPT = os.path.join(REPO_ROOT, 'paper1-murray', 'manuscript')
P2_MANUSCRIPT = os.path.join(REPO_ROOT, 'paper2-variational', 'manuscript')


def run_paper2():
    """Compute Paper 1 (Murray) values and write dynamic_variables.tex."""
    print("=" * 60)
    print("Paper 1 (Murray Cost Functions)")
    print("=" * 60)
    
    import compute_paper_murray as p2
    importlib.reload(p2)
    
    results = p2.compute_all()
    
    # Write to monorepo location
    tex_path = os.path.join(P1_MANUSCRIPT, 'dynamic_variables.tex')
    p2.write_tex(results, tex_path)
    
    fig_path = os.path.join(P1_MANUSCRIPT, 'figures', 'fig_alpha_scale.pdf')
    p2.generate_figure(fig_path)
    
    print(f"  alpha* = [{results['alpha_mw_low']:.3f}, {results['alpha_mw_high']:.3f}]")
    print(f"  Delta-alpha = {abs(results['delta_alpha_ours']):.2f}")
    print()
    return results


def run_paper1():
    """Compute Paper 2 (Variational) values and write dynamic_variables.tex."""
    print("=" * 60)
    print("Paper 2 (Variational Principle, Two-Level Model)")
    print("=" * 60)
    
    import compute_paper_variational as p1
    importlib.reload(p1)

    eta = p1.get_baseline_eta()
    al_mean, _, _ = p1.alpha_local_mean()
    alpha_star, alpha_t = p1.find_alpha_star(eta)
    pred_error = abs(alpha_star - 2.7) / 2.7 * 100
    wave_cost_net = p1.wave_loss_network(alpha_star) * 100

    gamma_table = p1.compute_gamma_table()
    eta_table = p1.compute_eta_parametric_table()
    kappa_table = p1.compute_kappa_table(eta)

    results = {
        'eta': eta,
        'alpha_local_mean': al_mean,
        'alpha_t': alpha_t,
        'alpha_star': alpha_star,
        'pred_error': pred_error,
        'wave_cost_net': wave_cost_net,
        'gamma_table': gamma_table,
        'eta_table': eta_table,
        'kappa_table': kappa_table,
    }
    
    tex_path = os.path.join(P2_MANUSCRIPT, 'dynamic_variables.tex')
    p1.write_tex(results, tex_path)
    
    fig_dir = os.path.join(P2_MANUSCRIPT, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    p1.generate_figures(results, fig_dir)
    
    print(f"  alpha_local = {al_mean:.3f}")
    print(f"  alpha_t     = {alpha_t:.3f}")
    print(f"  alpha*      = {alpha_star:.3f}")
    print(f"  error       = {pred_error:.1f}%")
    print()
    return results


def verify(p1_results, p2_results):
    """Cross-check consistency between papers."""
    print("=" * 60)
    print("VERIFICATION")
    print("=" * 60)
    
    checks = []
    
    # Paper 1: alpha* in [2.885, 3.000]
    a_low = p1_results['alpha_mw_low']
    a_high = p1_results['alpha_mw_high']
    ok = 2.885 < a_high < a_low < 3.0
    checks.append(('P1 alpha* bounds', ok, f'{a_high:.3f} < {a_low:.3f} < 3.000'))
    
    # Paper 2: alpha* in (2.0, 3.0)
    a_star = p2_results['alpha_star']
    ok = 2.0 < a_star < 3.0
    checks.append(('P2 alpha* range', ok, f'{a_star:.3f}'))
    
    # Paper 2: alpha* close to experiment
    ok = abs(a_star - 2.7) < 0.2
    checks.append(('P2 vs experiment', ok, f'|{a_star:.3f} - 2.70| = {abs(a_star-2.7):.3f} < 0.20'))
    
    # Cross-check: P2 alpha_t ~ P1 alpha_local
    al = p2_results['alpha_local_mean']
    ok = abs(al - a_low) < 0.1
    checks.append(('P2 alpha_local ~ P1 alpha*', ok, f'{al:.3f} ~ {a_low:.3f}'))
    
    # All eta values in [0,1]
    for row in p2_results['eta_table']:
        ok = 0.0 <= row['eta'] <= 1.0
        if not ok:
            checks.append((f"eta={row['eta']:.2f}", False, f"{row['eta']:.2f}"))
            break
    else:
        checks.append(('All eta in [0,1]', True, 'OK'))
    
    all_ok = True
    for name, ok, detail in checks:
        status = '[PASS]' if ok else '[FAIL]'
        print(f'  {status} {name}: {detail}')
        if not ok:
            all_ok = False
    
    print()
    if all_ok:
        print('  ALL CHECKS PASSED')
    else:
        print('  SOME CHECKS FAILED')
        sys.exit(1)


def run_supplemental_wrapper():
    """Compute all supplemental results (S1 transfer-matrix, S2 MC, S3 power budget)."""
    import compute_supplemental as cs
    importlib.reload(cs)
    return cs.run_supplemental()


if __name__ == '__main__':
    verify_mode = '--verify' in sys.argv

    p1_res   = run_paper2()
    p2_res   = run_paper1()
    supp_res = run_supplemental_wrapper()

    if verify_mode:
        verify(p1_res, p2_res)
    
    print()
    print("[OK] All computations complete.")
