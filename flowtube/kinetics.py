"""
Kinetics calculations for flow tube experiments.
"""

from .flow_calc import full_attrs, carrier_attrs
import numpy as np
from numpy.typing import NDArray
from scipy.stats import linregress


### KPS Method Calculations ###
def diffusion_limited_rate_constant(
    obj: full_attrs,
    N_eff_Shw: float,
    diameter: float,
) -> float:
    """
    Calculate diffusion limited rate constant (s-1) - eq. 10 from Knopf
    et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        N_eff_Shw (float): Effective Sherwood number.
        diameter (float): Diameter of the cylinder (cm).

    Returns:
        float: Diffusion limited rate constant (s-1).
    """

    return 4 * N_eff_Shw * obj.reactant_diffusion_rate / diameter**2


def diffusion_limited_uptake_coefficient(
    obj: full_attrs,
    diameter: float,
    k_diff: float,
) -> float:
    """
    Calculate diffusion limited effective uptake coefficient - eq. 19
    from Knopf et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        diameter (float): Diameter of the cylinder (cm).
        k_diff (float): Diffusion limited rate constant (s-1).

    Returns:
        float: Diffusion limited effective uptake coefficient (cm s-1).
    """

    return diameter / obj.reactant_molec_velocity * k_diff


def correction_factor(
    N_eff_Shw: float,
    Kn: float,
    gamma: NDArray[np.float64] | float,
) -> NDArray[np.float64] | float:
    """
    Calculate correction factor for uptake coefficient - eq. 20 from
    Knopf et al., 2015.

    Args:
        N_eff_Shw (float): Effective Sherwood number (unitless).
        Kn (float): Knudsen number (unitless).
        hypothetical_gamma (float): Hypothetical uptake coefficient
            (unitless).

    Returns:
        float: Correction factor (unitless).
    """

    return 1 / (1 + gamma * 3 / (2 * N_eff_Shw * Kn))


def observed_loss_rate(
    obj: full_attrs,
    diameter: float,
    gamma_eff: NDArray[np.float64] | float,
) -> NDArray[np.float64] | float:
    """
    Calculate observed loss rate (s-1) - eq. 19 from Knopf et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        diameter (float): Diameter of the cylinder (cm).
        gamma_eff (float): Effective uptake coefficient (unitless).

    Returns:
        float: Observed loss rate (s-1).
    """

    return gamma_eff * obj.reactant_molec_velocity / diameter


def cylinder_loss(
    obj: full_attrs,
    diameter: float,
    N_eff_Shw: float,
    Kn: float,
    gamma: NDArray[np.float64] | float,
    time: float,
) -> NDArray[np.float64] | float:
    """
    Calculate penetration (unitless) - eq. 21 from Knopf et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        time (float): Residence time in cylinder (s).

    Returns:
        float: Penetration - fraction of initial reactant after passing
            through cylinder (unitless).
    """

    return 1 - np.exp(
        -gamma
        / (1 + gamma * 3 / (2 * N_eff_Shw * Kn))
        * obj.reactant_molec_velocity
        / diameter
        * time
    )


### Uptake Kinetics Calculations ###
def gamma_from_k(
    obj: full_attrs,
    k: float,
    diameter: float,
) -> float:
    """
    Calculate effective uptake coefficient from rate constant.

    Args:

        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        k (float): Rate constant (s-1).
        diameter (float): Diameter of the cylinder (cm).

    Returns:
        float: Effective uptake coefficient (unitless).
    """

    return k * diameter / obj.reactant_molec_velocity


def fit_first_order_kinetics(
    obj: carrier_attrs,
    concentrations: NDArray[np.float64],
    exposure: NDArray[np.float64],
    exposure_units: str,
) -> tuple[float, float, float, float, float]:
    """
    Fits the observed loss to a first order kinetic model to extract the 
    uptake coefficient.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        concentrations (NDArray[np.float64]): Array of observed
            concentrations (unitless).
        exposure (NDArray[np.float64]): Array of exposures (s or cm).
        exposure_units (str): Units of exposure ("s", "sec", "second", 
            "seconds", "cm", "centimeter", "centimeters").

    Returns:
        tuple[float, float, float, float, float]: Tuple containing the
            slope, intercept, r-value, p-value, and standard error of 
            the regression.
    """

    ### Check for valid inputs ###
    if not isinstance(concentrations, np.ndarray):
        raise TypeError("Concentrations must be a numpy array")
    if not isinstance(exposure, np.ndarray):
        raise TypeError("Exposure must be a numpy array")
    if (concentrations < 0).any():
        raise ValueError("Concentrations must be non-negative")
    
    # Find which flow velocity to use
    if hasattr(obj, "insert_flow_velocity"):
        flow_velocity = obj.insert_flow_velocity # pyright: ignore[reportAttributeAccessIssue]
    elif hasattr(obj, "flow_velocity"):
        flow_velocity = obj.flow_velocity # pyright: ignore[reportAttributeAccessIssue]
    elif hasattr(obj, "FT_flow_velocity"):
        flow_velocity = obj.FT_flow_velocity # pyright: ignore[reportAttributeAccessIssue]
    else:
        raise RuntimeError("Must call initialize() prior to fitting")
    
    # Convert exposure to time
    if exposure_units in ["s", "sec", "second", "seconds"]:
        exposure_time = exposure
    elif exposure_units in ["cm", "centimeter", "centimeters"]:
        exposure_time = exposure / flow_velocity
    else:
        raise ValueError(
            "Unsupported exposure units. "
            "Supported units: s, sec, second, seconds, cm, centimeter, centimeters"
        )

    # Fit data with linear regression
    log_concentrations = np.log(concentrations)

    return linregress(exposure_time, log_concentrations)


""" 
Citations:

Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and 
Penetration of Gas Molecules and Aerosol Particles through Laminar Flow 
Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87, 3746–3754. 
https://doi.org/10.1021/ac5042395
"""