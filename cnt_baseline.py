"""Homogeneous supercooling of pure gallium from the SI Section S1.1 relations.

Reproduces the calibration baseline from the published equations and Table S1 inputs,
independent of the engine. Used to propagate the corrected solid molar volume:
V_m = 11.81e-6 m^3/mol (CRC solid alpha-Ga) with the Turnbull interfacial energy
(C_T = 0.43) gives gamma = 55 mJ/m^2 and a predicted supercooling of ~68 K, within
1 K of the measured 67.8 +/- 1.1 K (Zhang et al.).

    python cnt_baseline.py
"""
import numpy as np
from scipy.optimize import brentq

dHf = 5590.0        # molar latent heat of fusion, J/mol (NIST/CRC)
Tliq = 302.9        # liquidus, K
kB = 1.380649e-23
NA = 6.02214076e23
I0 = 1e39           # rate prefactor, m^-3 s^-1
I_THRESHOLD = 1e6   # activated-rate threshold, m^-3 s^-1
C_T = 0.43          # Turnbull coefficient (metallic; the value the paper's gamma implies)


def turnbull_gamma(Vm):
    """Solid-liquid interfacial energy from the Turnbull correlation, J/m^2."""
    return C_T * dHf / (NA ** (1 / 3) * Vm ** (2 / 3))


def activated_supercooling(Vm, gamma, f=1.0):
    """First undercooling at which the steady-state nucleation rate crosses the threshold."""
    def rate_gap(dT):
        dGv = (dHf / Vm) * (dT / Tliq)
        barrier = f * 16 * np.pi * gamma ** 3 / (3 * dGv ** 2)
        return I0 * np.exp(-barrier / (kB * (Tliq - dT))) - I_THRESHOLD
    return brentq(rate_gap, 1.0, 150.0)


if __name__ == "__main__":
    Vm = 11.81e-6                 # CRC density of solid alpha-Ga, 5.904 g/cm^3
    gamma = turnbull_gamma(Vm)
    dT = activated_supercooling(Vm, gamma)
    print(f"V_m   = {Vm*1e6:.2f} x 10^-6 m^3/mol (CRC solid alpha-Ga)")
    print(f"gamma = {gamma*1e3:.0f} mJ/m^2 (Turnbull, C_T = {C_T}; {gamma*1e3:.1f} unrounded)")
    print(f"predicted homogeneous supercooling = {dT:.0f} K ({dT:.1f} unrounded)")
    print(f"measured (Zhang et al.)            = 67.8 +/- 1.1 K  ->  agreement within 1 K")
