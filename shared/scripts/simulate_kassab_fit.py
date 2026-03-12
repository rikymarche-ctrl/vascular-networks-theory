"""
simulate_kassab_fit.py
Simulates Kassab's morphometric extraction of the branching exponent alpha on a theoretically optimal, asymmetric tree.

Methodology:
1. Generate an asymmetric tree of G generations, distributing flow purely via mass conservation.
2. Calculate theoretically optimal radii r*(Q) for every segment independently according to Paper 1 (Three-term cost).
3. Apply Kassab's original global least squares fit on the node residuals  (d_parent^a - d_child1^a - d_child2^a)^2
4. Compare the globally fitted alpha_morph to the flow-weighted mean of alpha_local.
"""

import numpy as np
from scipy.optimize import minimize_scalar
import sys
import os

sys.path.append(r'C:\Users\ricca\Documenti\Progetti IA\branching-networks\shared\scripts')
from params import *
from compute_paper1 import r_star, alpha_local_mean

def generate_tree(Q0, f_asym, G):
    """
    Generates a list of all bifurcations in the tree via pure flow conservation.
    f_asym: flow asymmetry parameter (Q1/Q0).
    Returns list of tuples: (Q_parent, Q_child1, Q_child2)
    """
    nodes = []
    current_generation = [Q0]
    
    for _ in range(G):
        next_generation = []
        for q in current_generation:
            q1 = q * f_asym
            q2 = q * (1 - f_asym)
            nodes.append((q, q1, q2))
            next_generation.extend([q1, q2])
        current_generation = next_generation
        
    return nodes

def fit_morphometric_alpha(nodes, mw):
    """
    Applies global least squares to find the optimal alpha minimizing Kassab's residual.
    """
    # 1. First assign theoretically perfect radii to every flow
    radii_nodes = []
    for q0, q1, q2 in nodes:
        d0 = 2 * r_star(q0, mw)
        d1 = 2 * r_star(q1, mw)
        d2 = 2 * r_star(q2, mw)
        radii_nodes.append((d0, d1, d2))
        
    # 2. Define the Kassab Global Residual Sum of Squares globally
    def sum_squared_residuals(alpha):
        ssr = 0
        for d0, d1, d2 in radii_nodes:
            # We normalize by d0^alpha to prevent massive over-weighting of the root node
            residual = (d0**alpha - d1**alpha - d2**alpha) / (d0**alpha)
            ssr += residual**2
        return ssr
        
    # 3. Minimize
    res = minimize_scalar(sum_squared_residuals, bounds=(2.0, 4.0), method='bounded')
    return res.x, radii_nodes

def compute_weighted_local_mean(nodes, mw):
    """Computes the exact flow-weighted mean of the local theoretical alpha exponents."""
    sum_w_alpha = 0
    sum_w = 0
    
    for q0, q1, q2 in nodes:
        # Theoretical local alpha using the exact formula Alpha = log(2) / log(d0/d1) symmetrically
        # For asymmetric trees, local alpha_i solves d0^a = d1^a + d2^a
        d0 = 2 * r_star(q0, mw)
        d1 = 2 * r_star(q1, mw)
        d2 = 2 * r_star(q2, mw)
        
        def local_residual(a):
            return abs(d0**a - d1**a - d2**a) / (d0**a)
            
        a_local = minimize_scalar(local_residual, bounds=(2.0, 4.0), method='bounded').x
        
        # The weight implied by Kassab minimizing the relative residual variance 
        # is analytically derived from a Taylor expansion as the squared Shannon Entropy
        f_split = q1 / q0
        H_f = - (f_split * np.log(f_split) + (1 - f_split) * np.log(1 - f_split))
        w_i = (H_f / a_local)**2
        
        sum_w_alpha += w_i * a_local
        sum_w += w_i
        
    return sum_w_alpha / sum_w

if __name__ == '__main__':
    print("=" * 60)
    print("Simulating Kassab Morphometric Fit on Assymetric Theoretical Tree")
    print("=" * 60)
    
    G = 11  # Kassab Orders
    Q0 = Q0_coronary
    mw = MW_MID
    
    # Kassab reports significant asymmetry in coronary trees. We test across an asymmetry spectrum.
    # f = 0.5 is perfectly symmetric. f = 0.8 is highly asymmetric.
    f_values = [0.5, 0.6, 0.7, 0.8]
    
    print(f"Parameters: G={G}, Q0={Q0*1e6:.1f} mL/s, mw={mw/1e3:.0f} kW/m³")
    print("-" * 60)
    print(f"{'Asymmetry (f)':<15} | {'Alpha_morph (Fit)':<20} | {'Weighted Local Mean':<20}")
    print("-" * 60)
    
    for f in f_values:
        nodes = generate_tree(Q0, f, G)
        a_morph, _ = fit_morphometric_alpha(nodes, mw)
        a_mean = compute_weighted_local_mean(nodes, mw)
        
        print(f"f = {f:<11.2f} | {a_morph:<20.4f} | {a_mean:<20.4f}")
        
    print("=" * 60)
