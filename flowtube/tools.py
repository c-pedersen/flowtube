"""
Constants and unit conversions.

Constants:

- STANDARD_TEMPERATURE_K: Standard temperature in Kelvin.
- STANDARD_PRESSURE_Pa: Standard pressure in Pascal.
- UNIVERSAL_GAS_CONSTANT: Universal gas constant in kg m2 s-2 K-1 mol-1.
- BOLTZMANN_CONSTANT: Boltzmann constant in kg m2 s-2 K-1.
- AVOGADROS_NUMBER: Avogadro's number in mol-1.
- P_CF: Conversion factors to Pascal for Torr, bar, mbar, hPa, and Pa pressure units.

"""

from warnings import warn
import numpy as np
import pandas as pd
import molmass as mm
from tabulate import tabulate
from requests.structures import CaseInsensitiveDict


### Constants ###
# Physical constants
STANDARD_TEMPERATURE_K = 273.15  # K
STANDARD_PRESSURE_Pa = 101325  # Pa
UNIVERSAL_GAS_CONSTANT = 8.3145  # kg m2 s-2 K-1 mol-1
BOLTZMANN_CONSTANT = 1.380649e-23  # kg m2 s-2 K-1
AVOGADROS_NUMBER = 6.0221408e23  # mol-1

# Conversion factors and calculated constants
MOLES_PER_SCCM = (
    STANDARD_PRESSURE_Pa * (1e-2) ** 3 / UNIVERSAL_GAS_CONSTANT / STANDARD_TEMPERATURE_K
)

P_CF = CaseInsensitiveDict(
    {
        "Torr": 133.322,
        "bar": 1e5,
        "mbar": 100,
        "hPa": 100,
        "Pa": 1,
    }
)  # conversion factors to Pa


### Concentrations conversions ###
def permeation_rate_to_MR(
    flow_rate: float,
    permeation_rate: float,
    reactant_gas: str,
) -> float:
    """Convert permeation rate to mixing ratio.

    Args:
        flow_rate (float): Flow rate through permeation in sccm.
        permeation_rate (float): Permeation rate in ng/min.
        reactant_gas (str): Molecular formula of reactant gas.

    Returns:
        float: Reactant volumetric mixing ratio (mol mol-1).
    """
    # Molar mass (g mol-1)
    molar_mass = float(mm.Formula(reactant_gas).mass)

    # Convert permeation rate to mol / min
    mol_per_min = (permeation_rate * 1e-9) / molar_mass

    # Calculate mol / min of total flow
    total_mol_per_min = flow_rate * MOLES_PER_SCCM

    return mol_per_min / total_mol_per_min


def vapor_pressure_to_MR(
    vapor_pressure: float,
    P_units: str,
    system_pressure: float,
    P_units_system: str,
) -> float:
    """
    Convert vapor pressure to mixing ratio. Assumes the carrier gas
    becomes fully saturated with the volatile substance.

    Args:
        vapor_pressure (float): Vapor pressure.
        P_units (str): Units of vapor pressure.
        system_pressure (float): System pressure.
        P_units_system (str): Units of system pressure.

    Returns:
        float: Reactant volumetric mixing ratio (mol mol-1).
    """
    # Convert vapor pressure to Pa if necessary
    if P_units != "Pa":
        vapor_pressure_Pa = P_in_Pa(vapor_pressure, P_units)
    else:
        vapor_pressure_Pa = vapor_pressure

    # Convert system pressure to Pa if necessary
    if P_units_system != "Pa":
        system_pressure_Pa = P_in_Pa(system_pressure, P_units_system)
    else:
        system_pressure_Pa = system_pressure

    if vapor_pressure_Pa < 0 or not np.isfinite(vapor_pressure_Pa):
        raise ValueError("Vapor pressure must be non-negative and finite.")
    if system_pressure_Pa <= 0 or not np.isfinite(system_pressure_Pa):
        raise ValueError("System pressure must be positive and finite.")
    if vapor_pressure_Pa > system_pressure_Pa:
        raise ValueError("Vapor pressure cannot exceed system pressure.")
    if vapor_pressure_Pa > system_pressure_Pa * 0.01:
        warn("Vapor pressure is greater than 1% of system pressure. " \
        "This may lead to non-ideal behavior and inaccurate mixing ratio calculation.")

    # Calculate mixing ratio as the ratio of vapor pressure to system pressure
    return vapor_pressure_Pa / system_pressure_Pa


### Unit conversions ###
def T_in_K(
    T: float,
) -> float:
    """Convert Celsius to Kelvin.

    Args:
        T (float): Temperature in Celsius.

    Returns:
        float: Temperature in Kelvin.
    """

    return T + STANDARD_TEMPERATURE_K


def P_in_Pa(
    P: float,
    units: str,
) -> float:
    """Convert pressure to Pa.

    Args:
        P (float): Pressure.
        units (float): Supported pressure units: Torr, bar, mbar, or
            hPa.

    Returns:
        float: Pressure in Pa.
    """

    return P * P_CF[units]


### Geometric Calculations ###
def cross_sectional_area(
    diameter: float,
) -> float:
    """Calculate cross sectional area.

    Args:
        diameter (float): Diameter.

    Returns:
        float: Cross sectional area.
    """

    return np.pi * (diameter / 2) ** 2


def partial_cylinder_area(
    height: float,
    width: float,
) -> tuple[float, float]:
    """
    Calculate the cross-sectional area and perimeter of a partial
    cylinder ("boat") given its height and width. Assumes the partial
    cylinder is a segment of a circle with given width as chord length
    and that the height is less than the diameter of the circle.
    Calculations based on:
    https://mathworld.wolfram.com/CircularSegment.html and
    https://www.vcalc.com/wiki/KurtHeckman/Circle-area-of-an-arc-segment-h.

    Args:
        height (float): Height of the segment (cm).
        width (float): Chord length of the segment (cm).

    Returns:
        float: Outer perimeter (cm).
        float: Cross-sectional area (cm^2).
    """

    # Ensure the proper geometry for the following calculation
    if height > width / 2:
        raise ValueError("Boat height cannot be greater than boat radius.")

    # Calculate radius from height and chord length
    r = (height / 2) + (width**2 / (8 * height))

    # Central angle (theta) in radians
    theta = 2 * np.arccos((r - height) / r)

    # Arc length
    arc_length = r * theta

    # Wetted perimeter
    perimeter = arc_length + width

    # Area of the segment
    area = (r**2 / 2) * (theta - np.sin(theta))

    return perimeter, area


### Display Calculations ###
def table(
    title: str,
    var_names: list[str],
    var: list[float],
    var_fmts: list[str],
    units: list[str],
) -> None:
    """Print a formatted table of variables.

    Args:
        title (str): Title of the table.
        var_names (list): List of variable names.
        var (list): Variable values.
        var_fmts (list): List of formats for each variable.
        units (list): List of units for each variable.

    Returns:
        None
    """

    data = pd.DataFrame(
        [var_names, [format(v, "8" + fmt) for v, fmt in zip(var, var_fmts)], units]
    ).T
    table = tabulate(
        data,  # pyright: ignore[reportArgumentType]
        disable_numparse=True,
        tablefmt="fancy_grid",
        showindex=False,
    )

    # Center the title based on the table width
    width = len(table.splitlines()[0])
    print(f"\033[1m{title}\033[0m".center(width))
    print(table)
