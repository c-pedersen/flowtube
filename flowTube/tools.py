'''
Constants and unit conversions.
'''

import math
import pandas as pd
from tabulate import tabulate # pyright: ignore[reportMissingModuleSource]
from requests.structures import CaseInsensitiveDict

### Unit conversions ###
def T_in_K(T: float) -> float:
    ''' Convert celcius to Kelvin.

    Args:
        T (float): Temperature in Celcius.

    Returns:
        float: Temperature in Kelvin.
    '''

    return T + STANDARD_TEMPERATURE_K

def P_in_Pa(P: float, units: str) -> float:
    ''' Convert pressure to Pa.

    Args:
        P (float): Pressure.
        units (float): Supported pressure units: Torr, bar, mbar, or hPa.

    Returns:
        float: Pressure in Pa.
    '''

    return P * P_CF[units] 

### Geometric Calculations ###
def cross_sectional_area(diameter: float) -> float:
    ''' Calculate cross sectional area.

    Args:
        diameter (float): Diameter.

    Returns:
        float: Cross sectional area.
    '''

    return math.pi * (diameter / 2) ** 2

### Display Calculations ###
def table(var_names: list[str], var: list[float], var_fmts: list[str], units: list[str]) -> None:
    ''' Print a formatted table of variables.

    Args:
        var_names (list): List of variable names.
        var (list): Variable values.
        var_fmts (list): List of formats for each variable.
        units (list): List of units for each variable.

    Returns:
        None
    '''

    table = pd.DataFrame([var_names, [format(v, '8'+fmt) for v, fmt in zip(var, var_fmts)], units]).T
    print(tabulate(table, disable_numparse = True, tablefmt = 'fancy_grid', showindex=False)) # type: ignore

### Constants ###
STANDARD_TEMPERATURE_K = 273.15 # K
STANDARD_PRESSURE_Pa = 101325 # Pa
UNIVERSAL_GAS_CONSTANT = 8.3145 # kg m2 s-2 K-1 mol-1
BOLTZMANN_CONSTANT = 1.380649e-23 # kg m2 s-2 K-1
AVOGADROS_NUMBER = 6.0221408e23 # mol-1
P_CF = CaseInsensitiveDict({'Torr': 133.322, 'bar': 1e5, 'mbar': 100, 'hPa': 100, 'Pa': 1}) # conversion factors to Pa