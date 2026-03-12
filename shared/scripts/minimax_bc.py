import numpy as np
from scipy.optimize import root_scalar
import sys
import os

# Ensure we can import physical constants from params.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from params import *
from compute_paper2 import alpha_star

class VascularMinimaxEngine:
    """
    Motore Morfologico basato sul Principio Variazionale Unificato.
    Calcola l'impedenza periferica perfetta per un dato paziente.
    """
    def __init__(self, p_exponent=p, stiffness_ke=1.0, mw=MW_MID, G=11, Q_in=Q0_coronary):
        """
        Inizializza il motore con i parametri anatomici e fisiologici.
        
        Parametri:
        - p_exponent: Costo metabolico della parete (esponenziale di Marchesi, difetto anatomico)
        - stiffness_ke: Modulo di rigidità della parete (fattore patologico per alterare l'impedenza)
        - mw: Costo metabolico del tessuto parietale (W/m^3)
        - G: Profondità della rete (numero di generazioni)
        - Q_in: Flusso di ingresso alla rete (mL/s)
        """
        self.p = p_exponent
        self.ke = stiffness_ke
        self.mw = mw
        self.G = G
        self.Q_in = Q_in
        
        # 1. Attrattore statico locale (Viscoso + Metabolico) -> ~2.92
        self.alpha_t = alpha_star(self.Q_in, self.mw)
        
        # 2. Attrattore d'onda (Impedance Matching) -> ~2.115 per paziente sano
        # Z proporzionale a r^(-(5-p)/2). Con k_e = 1.0, alpha_w = (5-0.77)/2 = 2.115
        self.alpha_w = ((5 - self.p) / 2) * self.ke

    def cost_transport_net(self, alpha):
        """
        Penalità di Trasporto di Rete C_transport^{net}(alpha).
        Per il PoC, usiamo l'espansione di Taylor (approssimazione armonica)
        attorno al minimo teorico locale alpha_t. La concavità (peso) 
        andrebbe calibrata sulle derivate seconde esatte di Phi.
        """
        # Approssimazione quadratica: k_t * (alpha - alpha_t)^2
        k_t = 1.0 # Da sostituire con (d^2 Phi_net / d_alpha^2) esatto
        return k_t * (alpha - self.alpha_t)**2

    def cost_wave_net(self, alpha):
        """
        Penalità di Onde di Rete C_wave^{net}(alpha).
        L'energia dell'onda retrograda è minima quando la geometria
        rispetta il limite per la conservazione dell'area-adattiva alpha_w.
        """
        # Approssimazione quadratica: k_w * (alpha - alpha_w)^2
        k_w = 1.0 # Da sostituire con la severità della riflessione d'onda
        return k_w * (alpha - self.alpha_w)**2

    def lagrangian(self, alpha, eta):
        """
        L'Equazione di Rete (La Lagrangiana).
        Bilancia sforzo statico e pulsatile.
        L_net = eta * C_wave + (1 - eta) * C_transport
        """
        return eta * self.cost_wave_net(alpha) + (1 - eta) * self.cost_transport_net(alpha)

    def find_optimal_alpha(self):
        """
        Risolutore Minimax: Trova l'esponente ottimo alpha*.
        Il minimax stabilisce che al punto di ottimo le due penalità si eguagliano:
        C_wave_net(alpha) = C_transport_net(alpha)
        per annullare la dipendenza dalle fluttuazioni ambientali (eta).
        """
        def cost_difference(alpha):
            return self.cost_transport_net(alpha) - self.cost_wave_net(alpha)
        
        # Troviamo la radice (l'intersezione) tra alpha_w e alpha_t
        try:
            sol = root_scalar(cost_difference, bracket=[self.alpha_w, self.alpha_t], method='brentq')
            return sol.root
        except ValueError:
            # Fallback se non c'è radice nel bracket
            raise ValueError("Il solver non è riuscito a trovare un incrocio variazionale nel dominio atteso.")

    def find_optimal_eta(self, alpha_star):
        """
        Calcola la frazione di pulsatilità (duty cycle) implicita ottima \eta*.
        Si ottiene imponendo che dL/d(alpha) = 0.
        """
        # Derivata numerica semplice (Central Difference) per il PoC
        h = 1e-6
        d_C_trans = (self.cost_transport_net(alpha_star + h) - self.cost_transport_net(alpha_star - h)) / (2 * h)
        d_C_wave = (self.cost_wave_net(alpha_star + h) - self.cost_wave_net(alpha_star - h)) / (2 * h)
        
        # dL/d_alpha = eta * d_C_wave + (1-eta) * d_C_trans = 0
        # eta_star = d_C_trans / (d_C_trans - d_C_wave)
        if (d_C_trans - d_C_wave) == 0:
            return 0.5
            
        eta_star = d_C_trans / (d_C_trans - d_C_wave)
        return eta_star


# --- Simulazione e Test del Modulo ---
if __name__ == "__main__":
    print("="*60)
    print(" VASCULAR MINIMAX ENGINE - PROOF OF CONCEPT")
    print("="*60)
    
    # 1. Paziente Sano (Parametri fisiologici standard)
    sano = VascularMinimaxEngine(p_exponent=0.77, stiffness_ke=1.0)
    alpha_sano = sano.find_optimal_alpha()
    eta_sano = sano.find_optimal_eta(alpha_sano)
    
    print("\n[PAZIENTE SANO]")
    print(f"Attrattore Statico   (a_t): {sano.alpha_t:.3f} (dal paper 2, approx 2.92)")
    print(f"Attrattore Onda      (a_w): {sano.alpha_w:.3f} (dal paper 2, (5-p)/2)")
    print(f"-> Esponente Ottimo  (a^*): {alpha_sano:.3f}")
    print(f"-> Frazione Pulsatile(n^*): {eta_sano:.3f}")
    
    # 2. Paziente Patologico (Ipertensione / Irrigidimento)
    # k_e aumenta, alzando l'attrattore d'onda.
    stiffness_patho = 1.07
    paziente_pah = VascularMinimaxEngine(p_exponent=0.77, stiffness_ke=stiffness_patho)
    alpha_pah = paziente_pah.find_optimal_alpha()
    
    print("\n[PAZIENTE PATOLOGICO - Rigidità Aumentata del 7%]")
    print(f"Nuovo Attrattore Onda(a_w): {paziente_pah.alpha_w:.3f}")
    print(f"-> Nuovo Esponente   (a^*): {alpha_pah:.3f}")
    print(f"-> Delta Architettura     : +{alpha_pah - alpha_sano:.3f}")
    
    print("\nIl Motore Morfologico è pronto per essere agganciato al CFD per generare le Boundary Conditions.")
    print("="*60)
