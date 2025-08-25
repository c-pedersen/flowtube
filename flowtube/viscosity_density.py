'''
Handles the calculation of the viscosity and density of a variety of pure gases.
'''

import numpy as np
import molmass as mm # pyright: ignore[reportMissingImports]
import math
from typing import Protocol

from . import tools

class basic_attrs(Protocol):
    P: float # Pressure in Pa
    T: float # Temperature in K

## van der Waal's Constants - Chemistry LibreTexts
a = {'Ar': 1.355, 'He': 0.0346, 'N2': 1.370, 'O2': 1.382} # bar L2 mol-2
b = {'Ar': 0.03201, 'He': 0.0238, 'N2': 0.0387, 'O2': 0.03186} # L mol-1

## Constants from Appendix A of Reid et al., 1987
# Critical Temperature (K)
T_c = {'Ar': 150.8, 'He': 5.19, 'N2': 126.2, 'O2': 154.6}

# Critical Pressure (bar)
P_c = {'Ar': 48.7, 'He': 2.27, 'N2': 33.9, 'O2': 50.4}

# Dipole Moment (debye)
mu = {'Ar': 0.0, 'He': 0.0, 'N2': 0.0, 'O2': 0.0}

def real_density(obj: basic_attrs, gas: str) -> float:
    ''' Calculates the real density of a gas using van der Waal's equation of state.

    Args:
        obj (basic_attrs): Object with basic attributes (P in Pa, T in K).
        gas (str): Molecular formula of gas (supported: Ar, He, N2, O2).

    Returns:
        float: Real density of gas in kg m-3.
    '''

    # Check if the gas is supported
    if gas not in a.keys():
        raise ValueError(f"Unsupported gas. Supported gases: {', '.join(a.keys())}")
    
    # Molar mass (g mol-1)
    m = float(mm.Formula(gas).mass) # pyright: ignore[reportUnknownArgumentType, reportUnknownMemberType]
    
    # Calculate pressure from Pa to bar
    P_bar = obj.P / tools.P_CF['bar']

    # Coefficients for van der Waal's equation of state solved for the inverse density (cubic function)
    coefficients = [P_bar, -P_bar * b[gas] - tools.UNIVERSAL_GAS_CONSTANT / 100 * obj.T, a[gas], -a[gas] * b[gas]]

    # Find the roots of the cubic function
    roots = np.roots(coefficients) # pyright: ignore[reportUnknownMemberType]

    # Find the real root and convert to kg m-3
    density = 1 / np.real(roots[roots.imag == 0][0]) * m # pyright: ignore[reportUnknownMemberType]

    # Check if the density is negative, which is unphysical
    if density < 0:
        raise ValueError(f"Calculated density for {gas} at T={obj.T} K and P={obj.P} Pa is negative: {density} kg m-3")
    
    return density

def dynamic_viscosity(obj: basic_attrs, gas: str) -> float:
    ''' Estimate absolute/dynamic viscosity of pure gases using the corresponding states method from Reid et al., 1987.

    Args:
        obj (basic_attrs): Object with basic attributes (P in Pa, T in K).
        gas (str): Molecular formula of gas (supported: Ar, He, N2, O2).

    Returns:
        float: Gas viscosity (kg m-1 s-1).
    '''
    # Check if the gas is supported
    if gas not in a.keys():
        raise ValueError(f"Unsupported gas. Supported gases: {', '.join(a.keys())}")

    # Molar mass (g mol-1)
    m = float(mm.Formula(gas).mass) # pyright: ignore[reportUnknownArgumentType, reportUnknownMemberType]

    # Give warning if a polar gas is attempted
    if mu[gas] >= 0.022:
        UserWarning('Polar gases not supported. See eqs. 9-4.16 & 9-4.17 from Reid et al., 1987 to implement.')
    
    
    # Reduced temperature 
    T_r = obj.T / T_c[gas]
    
    # Reduced, Inverse Viscosity ((µP)-1) - eq. 9-4.14 from Reid et al., 1987
    xi = 0.176 * (T_c[gas] / (m ** 3 * P_c[gas] ** 4)) ** (1/6)
    
    # Quantum Gas Correction Factor - eq. 9-4.18 from Reid et al., 1987
    Q = {'He': 1.38, 'H2': 0.76, 'D2': 0.52}
    if gas in ['He', 'H2', 'D2']:
        f_Q = 1.22 * Q[gas] ** 0.15 * (1 + 0.00385 * (((T_r - 12) ** 2) ** (1 / m)) * math.copysign(1, T_r - 12))
    else:
        f_Q = 1
    
    # nu * xi (unitless) - eq. 9-4.15 from Reid et al., 1987
    nu_xi = (0.807 * T_r ** 0.618 - 0.357 * math.exp(-0.449 * T_r) + 0.34 * math.exp(-4.058 * T_r) + 0.018) * f_Q
    
    # Viscosity (kg m-1 s-1)
    nu = nu_xi / xi / 1e7

    return nu

'''
Citations

“A8: Van Der Waal’s Constants for Real Gases.” Chemistry LibreTexts, November 14, 2024. Accessed August 6, 2025. 
https://chem.libretexts.org/Ancillary_Materials/Reference/Reference_Tables/Atomic_and_Molecular_Properties/A8%3A_van_der_Waal’s_Constants_for_Real_Gases. 

“NIST Chemistry Webbook, SRD 69.” Thermophysical Properties of Fluid Systems. Accessed August 6, 2025. 
https://webbook.nist.gov/chemistry/fluid/. 

Reid, R.C., Prausnitz, J.M., Poling, B.E., 1987. The Properties of Gases and Liquids, 4th ed. McGraw-Hill, New York.
''' 
