'''
Handles the calculation of the viscosity and density of a variety of pure gases.
'''

import numpy as np
import molmass as mm
import math

from flowTube import tools

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

def real_density(gas, T, P):
    ''' Calculates the real density of a gas using van der Waal's equation of state.

    Args:
        gas (str): Molecular formula of gas.
        T (float): Temperature (K).
        P (float): Pressure (Pa).

    Returns:
        float: Real density of gas in kg m-3.
    '''
    
    # Calculate pressure from Pa to bar
    P_bar = P / tools.P_CF['bar']

    # Coefficients for van der Waal's equation of state solved for the inverse density (cubic function)
    coefficients = [P_bar, -P_bar * b[gas] - tools.R / 100 * T, a[gas], -a[gas] * b[gas]]

    # Find the roots of the cubic function
    roots = np.roots(coefficients)

    # Find the real root and convert to kg m-3
    density = 1/np.real(roots[roots.imag == 0][0]) * mm.Formula(gas).mass
    
    return density

def dynamic_viscosity(gas, T, P):
    ''' Estimates the absolute/dynamic viscosity of pure gases using the corresponding states method from Reid et al., 1987.
        Estimation of He resulted in a maximum of 2.3 % difference at 1 bar and 20 ˚C with NIST values. Lower temperatures
        result in a better percent difference and pressure only makes a very small difference.

    Args:
        gas (str): Molecular formula of gas.
        T (float): Temperature (K).
        P (float): Pressure (Pa).

    Returns:
        float: Gas viscosity (kg m-1 s-1).
    '''

    # Give warning if a polar gas is attempted
    if mu[gas] >= 0.022:
        warnings.warn('Polar gases not currently supported. See eqs. 9-4.16 & 9-4.17 from Reid et al., 1987 to implement.')
    
    # Reduced temperature 
    T_r = T / T_c[gas]
    
    # Reduced, Inverse Viscosity ((µP)-1) - eq. 9-4.14 from Reid et al., 1987
    xi = 0.176 * (T_c[gas] / (mm.Formula(gas).mass ** 3 * P_c[gas] ** 4)) ** (1/6)
    
    # Quantum Gas Correction Factor - eq. 9-4.18 from Reid et al., 1987
    Q = {'He': 1.38, 'H2': 0.76, 'D2': 0.52}
    if gas in ['He', 'H2', 'D2']:
        F_Q = 1.22 * Q[gas] ** 0.15 * (1 + 0.00385 * (((T_r - 12) ** 2) ** (1 / mm.Formula(gas).mass)) * math.copysign(1, T_r - 12))
    else:
        F_Q = 1
    
    # nu * xi (unitless) - eq. 9-4.15 from Reid et al., 1987
    nu_xi = (0.807 * T_r ** 0.618 - 0.357 * math.exp(-0.449 * T_r) + 0.34 * math.exp(-4.058 * T_r) + 0.018) * F_Q
    
    # Viscosity (kg m-1 s-1)
    nu = nu_xi / xi / 1e7

    return nu

'''
Citations

“A8: Van Der Waal’s Constants for Real Gases.” Chemistry LibreTexts, November 14, 2024. Accessed August 6, 2025. https://chem.libretexts.org/Ancillary_Materials/Reference/Reference_Tables/Atomic_and_Molecular_Properties/A8%3A_van_der_Waal’s_Constants_for_Real_Gases. 

“NIST Chemistry Webbook, SRD 69.” Thermophysical Properties of Fluid Systems. Accessed August 6, 2025. https://webbook.nist.gov/chemistry/fluid/. 

Reid, R.C., Prausnitz, J.M., Poling, B.E., 1987. The Properties of Gases and Liquids, 4th ed. McGraw-Hill, New York.
''' 
