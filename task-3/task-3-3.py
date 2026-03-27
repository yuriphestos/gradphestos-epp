#!/usr/bin/env python3
"""
(* ====================================================== *)
(* GitHub: https://github.com/yuriphestos/gradphestos-epp *)
(* Graduate Level Elementary Particle Physics             *)
(* Random Tasks for Fun                                   *)
(* ====================================================== *)
(* Task 3                                                 *)
(* Task 3.3: Solar Neutrino-Electron Scattering           *)
(* ====================================================== *)
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import dblquad, quad

# ============================================================
# Boron-8 Neutrino Flux Shape Data
# Digitized from Bahcall tables
# Enu in MeV, shape values (to be normalized)
# ============================================================

data = [
    (2.17466, 1.0e6), (2.21575, 1.0e6), (2.28767, 1.0e6),
    (2.33904, 1.0e6), (2.40068, 1.0e6), (2.45205, 1.0e6),
    (2.50342, 1.0e6), (2.56507, 9.56098e5), (2.62671, 1.0e6),
    (2.68836, 1.04390e6), (2.73973, 1.04390e6), (2.79110, 1.08780e6),
    (2.85274, 1.13171e6), (2.93493, 1.21951e6), (2.99658, 1.21951e6),
    (3.05822, 1.30732e6), (3.10959, 1.39512e6), (3.18151, 1.48293e6),
    (3.24315, 1.57073e6), (3.32534, 1.65854e6), (3.40753, 1.79024e6),
    (3.47945, 1.92195e6), (3.55137, 1.96585e6), (3.62329, 2.05366e6),
    (3.68493, 2.14146e6), (3.76712, 2.31707e6), (3.83904, 2.49268e6),
    (3.93151, 2.58049e6), (4.00342, 2.75610e6), (4.08562, 2.88780e6),
    (4.13699, 3.01951e6), (4.21918, 3.19512e6), (4.31164, 3.32683e6),
    (4.37329, 3.45854e6), (4.45548, 3.67805e6), (4.53767, 3.80976e6),
    (4.63014, 3.98537e6), (4.70205, 4.20488e6), (4.77397, 4.38049e6),
    (4.85616, 4.55610e6), (4.91781, 4.68780e6), (4.97945, 4.81951e6),
    (5.06164, 5.03902e6), (5.14957, 5.0e6),    (5.25641, 5.13043e6),
    (5.40598, 5.21739e6), (5.55556, 5.39130e6), (5.72650, 5.43478e6),
    (5.89744, 5.60870e6), (6.04701, 5.69565e6), (6.21795, 5.82609e6),
    (6.34615, 5.91304e6), (6.53846, 6.0e6),    (6.66667, 6.04348e6),
    (6.83761, 6.08696e6), (7.00855, 6.13043e6), (7.15812, 6.17391e6),
    (7.32906, 6.21739e6), (7.47863, 6.26087e6), (7.64957, 6.26087e6),
    (7.79915, 6.21739e6), (7.94872, 6.26087e6), (8.07692, 6.21739e6),
    (8.20513, 6.17391e6), (8.31197, 6.13043e6), (8.46154, 6.08696e6),
    (8.58974, 5.95652e6), (8.73932, 5.82609e6), (8.88889, 5.73913e6),
    (8.99573, 5.60870e6), (9.08120, 5.39130e6), (9.20940, 5.26087e6),
    (9.29487, 5.04348e6), (9.38034, 4.86957e6), (9.44444, 4.65217e6),
    (9.52991, 4.43478e6), (9.59402, 4.26087e6), (9.61538, 4.04348e6),
    (9.72222, 3.82609e6), (9.78632, 3.65217e6), (9.89316, 3.30435e6),
    (9.93590, 3.08696e6), (10.0,    2.86957e6), (10.17241, 2.58049e6),
    (10.30172, 2.40488e6),(10.47414, 2.18537e6), (10.68966, 2.00976e6),
    (10.90517, 1.83415e6),(11.12069, 1.65854e6), (11.33621, 1.52683e6),
    (11.50862, 1.35122e6),(11.81034, 1.26341e6), (12.02586, 1.17561e6),
    (12.28448, 1.08780e6),
]

Enu_arr = np.array([d[0] for d in data])
flux_arr = np.array([d[1] for d in data])

# Raw interpolation of the shape
f_raw = interp1d(Enu_arr, flux_arr, kind='cubic', fill_value=0.0, bounds_error=False)

# Normalize shape to total B8 flux (BS05 SSM)
Phi_total = 5.69e6  # cm^-2 s^-1
shape_integral, _ = quad(f_raw, Enu_arr[0], Enu_arr[-1], limit=200)
norm = Phi_total / shape_integral

def flux(Enu):
    """dPhi/dEnu in cm^-2 s^-1 MeV^-1"""
    return norm * f_raw(Enu)

print("=" * 55)
print("  Flux normalization")
print("=" * 55)
print(f"  Shape integral     = {shape_integral:.4e}")
print(f"  Total B8 flux      = {Phi_total:.2e} cm^-2 s^-1")
print(f"  Norm factor        = {norm:.4f}")
print()

# ============================================================
# Couplings
# ============================================================

sw2 = 0.231
gL  = 2*sw2 - 1   # -0.538
gR  = 2*sw2        #  0.462

print("=" * 55)
print("  Weak couplings")
print("=" * 55)
print(f"  sin^2(theta_W) = {sw2}")
print(f"  gL = {gL:.3f}")
print(f"  gR = {gR:.3f}")
print()

# ============================================================
# Cross section ratio
# ============================================================

sig_ratio = (gL**2 + gR**2/3) / ((gL + 2)**2 + gR**2/3)
print("=" * 55)
print("  Cross-section comparison")
print("=" * 55)
print(f"  sigma(nu_mu) / sigma(nu_e) = {sig_ratio:.4f}")
print(f"  sigma(nu_e)  / sigma(nu_mu) = {1/sig_ratio:.2f}")
print()

# ============================================================
# Differential cross section dsigma/dEe
#
# dsigma/dEe = [GF^2 me / (2 pi)] * [(gL+2)^2 + gR^2 (1 - Ee/Enu)^2]
#
# Units: GF in GeV^-2, me in GeV
#        GF^2 * me = GeV^-5
#        Multiply by (hbar c)^2 = 0.3894e-27 GeV^2 cm^2
#        -> GeV^-3 * cm^2 = cm^2 / GeV^3 ... but we have GeV^-5 * GeV^2 = GeV^-3
#        Actually: [GF^2] = GeV^-4, [me] = GeV^1
#        GF^2 * me * (hbar c)^2 = GeV^-4 * GeV * GeV^2 * cm^2 = cm^2/GeV
#        Then multiply by 1e-3 to get cm^2/MeV (since Ee is in MeV)
# ============================================================

GF      = 1.1663788e-5   # GeV^-2
me_GeV  = 0.511e-3       # GeV
me_MeV  = 0.511          # MeV
hbarc2  = 0.3894e-27     # GeV^2 cm^2

prefactor = GF**2 * me_GeV * hbarc2 * 1e-3 / (2 * np.pi)  # cm^2 / MeV

print("=" * 55)
print("  Cross section")
print("=" * 55)
print(f"  Prefactor = {prefactor:.6e} cm^2/MeV")

def dsigma_dEe(Enu, Ee):
    """dsigma/dEe for nu_e - e scattering [cm^2/MeV]"""
    return prefactor * ((gL + 2)**2 + gR**2 * (1 - Ee/Enu)**2)

# Cross-check at Enu = 10 MeV
Tmax_10 = 2 * 10**2 / (me_MeV + 2*10)
sig_10, _ = quad(dsigma_dEe, 0, Tmax_10, args=(10,))  # note: quad(f, a, b, args) -> f(Ee, Enu)
# Actually need to be careful with argument order for quad
sig_10, _ = quad(lambda Ee: dsigma_dEe(10, Ee), 0, Tmax_10)
print(f"  sigma(nu_e, 10 MeV) = {sig_10:.4e} cm^2  (expect ~9.4e-44)")
print()

# ============================================================
# Kinematics
# ============================================================

def Tmax(Enu):
    """Maximum electron recoil energy [MeV]"""
    return 2 * Enu**2 / (me_MeV + 2*Enu)

Ee_threshold = 5.0  # MeV (Borexino analysis cut)

# Minimum Enu that can produce Ee > threshold
# From Tmax(Enu_min) = Ee_threshold:
# 2 Enu^2 / (me + 2 Enu) = Ee_threshold
# 2 Enu^2 - 2 Ee_threshold Enu - me Ee_threshold = 0
a_coeff = 2
b_coeff = -2 * Ee_threshold
c_coeff = -me_MeV * Ee_threshold
Enu_min = (-b_coeff + np.sqrt(b_coeff**2 - 4*a_coeff*c_coeff)) / (2*a_coeff)
Enu_max = 12.28  # upper edge of spectrum data

print("=" * 55)
print("  Kinematics")
print("=" * 55)
print(f"  Ee threshold = {Ee_threshold} MeV")
print(f"  Enu_min      = {Enu_min:.3f} MeV")
print(f"  Enu_max      = {Enu_max} MeV")
print()

# ============================================================
# Target electrons (100 ton pseudocumene C9H12)
# M = 120.19 g/mol, Z_mol = 66 electrons/molecule
# ============================================================

Ne = (1e8 / 120.19) * 6.022e23 * 66
print("=" * 55)
print("  Target")
print("=" * 55)
print(f"  Ne = {Ne:.4e} electrons per 100 ton")
print()

# ============================================================
# Event Rate Integration
# R = Ne * 86400 * integral[ dsigma/dEe * dPhi/dEnu, dEe dEnu ]
# Outer: Enu from Enu_min to Enu_max
# Inner: Ee  from Ee_threshold to Tmax(Enu)
# ============================================================

sec_per_day = 86400

def integrand(Ee, Enu):
    """dblquad convention: f(y, x) with x=outer, y=inner"""
    return dsigma_dEe(Enu, Ee) * flux(Enu)

R, err = dblquad(
    integrand,
    Enu_min, Enu_max,        # outer (Enu)
    lambda Enu: Ee_threshold,  # inner lower (Ee)
    lambda Enu: Tmax(Enu),    # inner upper (Ee)
    epsabs=1e-12, epsrel=1e-10
)

R_cpd = R * Ne * sec_per_day

print("=" * 55)
print("  RESULTS")
print("=" * 55)
print(f"  R_predicted (no oscillations) = {R_cpd:.4f} cpd / 100 ton")
print()

# ============================================================
# Comparison with Borexino
# ============================================================

R_obs = 0.13  # cpd/100t from 0808.2868 Figs 6, 7
ratio = R_obs / R_cpd

print(f"  R_observed  (Borexino)        = {R_obs} cpd / 100 ton")
print(f"  R_obs / R_pred                = {ratio:.4f}")
print()

# ============================================================
# Extract Pee
# R_obs = R_pred * [ Pee + (1-Pee) * sigma_mu/sigma_e ]
# ============================================================

Pee = (ratio - sig_ratio) / (1 - sig_ratio)

print(f"  sigma(nu_mu)/sigma(nu_e)      = {sig_ratio:.4f}")
print(f"  P_ee (extracted)              = {Pee:.4f}")
print(f"  P_ee (MSW-LMA, ~ sin^2 th12) ~ 0.30 - 0.35")
print("=" * 55)
