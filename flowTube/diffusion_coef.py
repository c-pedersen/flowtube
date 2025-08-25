'''
Handles the calculation of the diffusion coefficient of a binary gas mixture.
'''

import numpy as np
import molmass as mm # pyright: ignore[reportMissingImports]
from typing import Protocol

from . import tools

class required_attrs(Protocol):
    P: float # Pressure in Pa
    T: float # Temperature in K
    reactant_gas: str # Molecular formula of reactant gas
    carrier_gas: str # Molecular formula of carrier gas

## Physical Constants (Appendix B from Reid et al., 1987)
# Characteristic Lennard-Jones Lengths (Å)
sigmas: dict[str, float] = {'Ar': 3.542, 'He': 2.551, 'Air': 3.711, 'Br2': 4.296, 'Cl2': 4.217, 'HBr': 3.353, 
                            'HCl': 3.339, 'HI': 4.211, 'H2O': 2.641, 'I2': 5.160, 'NO': 3.492, 'N2': 3.798, 'O2': 3.467}

# Characteristic Lennard-Jones Energies (K)
e_ks: dict[str, float] = {'Ar': 93.3, 'He': 10.22, 'Air': 78.6, 'Br2': 507.9, 'Cl2': 316.0, 'HBr': 449, 'HCl': 344.7, 
                          'HI': 288.7, 'H2O': 809.1, 'I2': 474.2, 'NO': 116.7, 'N2': 71.4, 'O2': 106.7}

def non_polar_Lennard_Jones_potential(e_k: float, T: float) -> float:
    ''' Calculation of non-polar Lennard-Jones Potential for a binary gas mixture.
        Formulas 11-3.4 to 11-3.6 in Reid et al., 1987

    Args:
        e_k (float): Lennard-Jones Energy (K).
        T (float): Temperature in K.

    Returns:
        float: Non-polar Lennard Jones potential.
    '''

    T_star = T / e_k
    
    return (1.06036 / T_star ** 0.1561 + 0.193 / np.exp(0.47635 * T_star) + 1.03587 / 
            np.exp(1.52996 * T_star) + 1.76474 / np.exp(3.89411 * T_star))

def binary_diffusion_coefficent(obj: required_attrs) -> float:
    ''' Calculation of non-polar diffusion coefficient for a low pressure binary gas mixture. Supported gases: Ar, He, 
        Air, Br2, Cl2', HBr, HCl, HI, H2O, I2, NO, N2, and O2. Formulas 11-3.1 to 11-3.2 in Reid et al., 1987

    Args:
        object (required_attrs): Object with required attributes (P in Pa, T in K, reactant_gas, carrier_gas).

    Returns:
        float: Diffusion coefficient for binary gas mixture (cm2 s-1)
    '''
    if (obj.reactant_gas not in sigmas.keys()) or (obj.carrier_gas not in sigmas.keys()):
        raise ValueError(f"Unsupported gas. Supported gases: {', '.join(sigmas.keys())}")
    
    # Combined molar mass (g mol-1)
    m1 = float(mm.Formula(obj.reactant_gas).mass) # pyright: ignore[reportUnknownArgumentType, reportUnknownMemberType]
    m2 = float(mm.Formula(obj.carrier_gas).mass) # pyright: ignore[reportUnknownArgumentType, reportUnknownMemberType]
    m: float = 2 / (1 / m1 + 1 / m2)

    # Mean Lennard-Jones Length (Å)
    mean_sigma = (sigmas[obj.reactant_gas] + sigmas[obj.carrier_gas]) / 2

    # Mean Lennard-Jones Energy (K)
    mean_e_k = (e_ks[obj.reactant_gas] * e_ks[obj.carrier_gas]) ** 0.5 
    
    # Diffusion Collision Integral (unitless)
    Omega_D = non_polar_Lennard_Jones_potential(mean_e_k, obj.T)
    
    return float(0.00266 * obj.T ** 1.5 / ((obj.P / tools.STANDARD_PRESSURE_Pa) * m ** 0.5 * mean_sigma ** 2 * Omega_D))

'''
Citations

Reid, R.C., Prausnitz, J.M., Poling, B.E., 1987. The Properties of Gases and Liquids, 4th ed. McGraw-Hill, New York.
'''
