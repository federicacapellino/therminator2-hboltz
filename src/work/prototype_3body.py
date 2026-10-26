"""Reference prototype for the three-body pair functions (spec Sec. 3.1, tests T1.1-T1.2).

Gas at rest, one channel a -> 1 2 3 with constant matrix element, Boltzmann parent,
branching ratio B = 1 and degeneracy g_a = 1. Units: GeV.
"""
import numpy as np
from scipy import integrate

T = 0.155


def lam(x, y, z):
    return x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x


def pstar(M, m1, m2):
    return np.sqrt(max(lam(M * M, m1 * m1, m2 * m2), 0.0)) / (2 * M)


def Phi2(M, m1, m2):
    return 4 * np.pi * pstar(M, m1, m2) / M


def Phi3_mass(ma, m1, m2, m3):
    f = lambda M: M * Phi2(ma, M, m3) * Phi2(M, m1, m2)
    return integrate.quad(f, m1 + m2, ma - m3, limit=200)[0]


def Phi3_dalitz(ma, m1, m2, m3):
    def width(s12):  # (s23)_max - (s23)_min = 4 p2 p3 in the 12 rest frame
        r = np.sqrt(s12)
        E2 = (s12 - m1**2 + m2**2) / (2 * r)
        E3 = (ma**2 - s12 - m3**2) / (2 * r)
        return 4 * np.sqrt(max(E2**2 - m2**2, 0)) * np.sqrt(max(E3**2 - m3**2, 0))
    AD = integrate.quad(width, (m1 + m2) ** 2, (ma - m3) ** 2, limit=200)[0]
    return 2 * np.pi**2 * AD / ma**2


f_parent = lambda e: np.exp(-e / T)            # Boltzmann parent distribution
tq, wq = np.polynomial.legendre.leggauss(64)   # Gauss-Legendre nodes on [-1, 1]


def CK_Cu(M, kap, ma, m3, D):
    """Spec Sec. 3.1 (LCS Eq. 18 with eps_P = (eK EP* + kap p3* t)/M)."""
    eK = np.sqrt(M * M + kap * kap)
    EP = (ma**2 + M**2 - m3**2) / (2 * M)
    p3 = np.sqrt(lam(ma**2, M**2, m3**2)) / (2 * M)
    fP = f_parent((eK * EP + kap * p3 * tq) / M)
    I0, I1 = np.sum(wq * fP), np.sum(wq * tq * fP)
    CK = 2 * np.pi * D * p3 / (kap * M**2) * (kap * EP * I0 + eK * p3 * I1)
    Cu = -2 * np.pi * D * p3**2 / (kap * eK) * I1
    return CK, Cu


if __name__ == "__main__":
    ma, m1, m2, m3 = 0.78266, 0.13957, 0.13957, 0.13498   # omega -> pi+ pi- pi0
    P3a, P3b = Phi3_mass(ma, m1, m2, m3), Phi3_dalitz(ma, m1, m2, m3)
    print(f"T1.1  Phi3: mass integral {P3a:.12g}, Dalitz area {P3b:.12g}, rel. diff {abs(P3a / P3b - 1):.1e}")

    D = 1.0 / ((2 * np.pi) ** 3 * P3a)                    # B = g_a = 1
    n_a = integrate.quad(lambda k: 4 * np.pi * k * k * f_parent(np.hypot(ma, k)), 0, 10)[0] / (2 * np.pi) ** 3
    worst = 0.0
    for M in [0.30, 0.40, 0.50, 0.60]:
        rho = 16 * np.pi**2 * pstar(M, m1, m2) * integrate.quad(
            lambda k: k * k * sum(CK_Cu(M, k, ma, m3, D)), 1e-6, 8, limit=200)[0]
        dalitz = n_a * M * Phi2(M, m1, m2) * Phi2(ma, M, m3) / P3a   # IM Eq. 9
        worst = max(worst, abs(rho / dalitz - 1))
        print(f"T1.2  M = {M:.2f} GeV: rho(M) from C = {rho:.10e}, Dalitz projection = {dalitz:.10e}")
    print(f"T1.2  worst relative deviation {worst:.1e}  ->  {'PASS' if worst < 1e-6 else 'FAIL'}")
    print("T1.6  small kappa (finite, positive):", [tuple(f"{c:.4e}" for c in CK_Cu(0.4, k, ma, m3, D)) for k in (1e-4, 1e-2, 0.5)])
