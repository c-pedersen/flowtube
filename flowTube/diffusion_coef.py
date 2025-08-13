'''
Handles the calculation of the diffusion coefficient of a binary gas mixture.
'''

import math
import molmass as mm

from flowTube import tools

## Physical Constants (Appendix B from Reid et al., 1987)
# Characteristic Lennard-Jones Lengths (Å)
sigmas = {'Ar': 3.542, 'He': 2.551, 'Air': 3.711, 'Br2': 4.296, 'Cl2': 4.217, 'HBr': 3.353, 'HCl': 3.339, 'HI': 4.211, 'H2O': 2.641, 'I2': 5.160, 'NO': 3.492, 'N2': 3.798, 'O2': 3.467}

# Characteristic Lennard-Jones Energies (K)
e_ks = {'Ar': 93.3, 'He': 10.22, 'Air': 78.6, 'Br2': 507.9, 'Cl2': 316.0, 'HBr': 449, 'HCl': 344.7, 'HI': 288.7, 'H2O': 809.1, 'I2': 474.2, 'NO': 116.7, 'N2': 71.4, 'O2': 106.7}

def non_polar_Lennard_Jones_potential(e_k, T):
    ''' Calculation of non-polar Lennard-Jones Potential for a binary gas mixture.
        Formulas 11-3.4 to 11-3.6 in Reid et al., 1987

    Args:
        e_k (float): Lennard-Jones Energy (K).
        T (float): Temperature in K.

    Returns:
        float: Non-polar Lennard Jones potential.
    '''

    T_star = T / e_k
    
    return (1.06036 / T_star ** 0.1561 + 0.193 / math.exp(0.47635 * T_star) + 1.03587 / 
            math.exp(1.52996 * T_star) + 1.76474 / math.exp(3.89411 * T_star))

def binary_diffusion_coefficent(gas1, gas2, T, P):
    ''' Calculation of non-polar diffusion coefficient for a low pressure binary gas mixture.
        Formulas 11-3.1 to 11-3.2 in Reid et al., 1987

    Args:
        gas1 (str): Molecular formula of gas 1.
        gas2 (str): Molecular formula of gas 2.
        T (float): Temperature in K.
        P (float): Pressure in Pa.

    Returns:
        float: Diffusion coefficient for binary gas mixture (cm2 s-1)
    '''
    
    # Combined molar mass (g mol-1)
    M = 2 / (1 / mm.Formula(gas1).mass + 1 / mm.Formula(gas2).mass)
    
    # Mean Lennard-Jones Length (Å)
    mean_sigma = (sigmas[gas1] + sigmas[gas2]) / 2

    # Mean Lennard-Jones Energy (K)
    mean_e_k = (e_ks[gas1] * e_ks[gas2]) ** 0.5
    
    # Diffusion Collision Integral (unitless)
    Omega_D = non_polar_Lennard_Jones_potential(mean_e_k, T)
    
    return (0.00266 * T ** 1.5 / ((P / tools.standard_P) * M ** 0.5 * mean_sigma ** 2 * Omega_D))

'''
Citations

Reid, R.C., Prausnitz, J.M., Poling, B.E., 1987. The Properties of Gases and Liquids, 4th ed. McGraw-Hill, New York.
'''
