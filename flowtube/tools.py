"""
Constants and unit conversions.
"""

import numpy as np
import pandas as pd
from tabulate import tabulate  # pyright: ignore[reportMissingModuleSource]
from requests.structures import CaseInsensitiveDict


### Unit conversions ###
def T_in_K(
    T: float,
) -> float:
    """Convert celcius to Kelvin.

    Args:
        T (float): Temperature in Celcius.

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
        data,
        disable_numparse=True,
        tablefmt="fancy_grid",
        showindex=False,
    )

    # Center the title based on the table width
    width = len(table.splitlines()[0])
    print(f"\033[1m{title}\033[0m".center(width))
    print(table)


### Constants ###
STANDARD_TEMPERATURE_K = 273.15  # K
STANDARD_PRESSURE_Pa = 101325  # Pa
UNIVERSAL_GAS_CONSTANT = 8.3145  # kg m2 s-2 K-1 mol-1
BOLTZMANN_CONSTANT = 1.380649e-23  # kg m2 s-2 K-1
AVOGADROS_NUMBER = 6.0221408e23  # mol-1
P_CF = CaseInsensitiveDict(
    {
        "Torr": 133.322,
        "bar": 1e5,
        "mbar": 100,
        "hPa": 100,
        "Pa": 1,
    }
)  # conversion factors to Pa
