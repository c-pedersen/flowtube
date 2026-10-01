from __future__ import annotations

from typing import TYPE_CHECKING
from warnings import warn

import molmass as mm
import numpy as np

from flowtube import diffusion_coef, tools, viscosity_density

if TYPE_CHECKING:
    from .aerosol_flow_reactor import AerosolFlowReactor
    from .boat_reactor import BoatReactor
    from .coated_wall_reactor import CoatedWallReactor


def validate_init(obj: AerosolFlowReactor | CoatedWallReactor | BoatReactor) -> None:
    """Validate inputs for init of all reactor types.

    Args:
        obj: An instance of a reactor class (AerosolFlowReactor,
            CoatedWallReactor, or BoatReactor).

    Returns:
        None. Raises errors if any validation checks fail.
    """
    from .boat_reactor import BoatReactor
    from .coated_wall_reactor import CoatedWallReactor

    # Check if the gases are supported
    if obj.reactant_gas not in diffusion_coef.sigmas:
        # Validate molecular formulas using molarmass
        try:
            mm.Formula(  # noqa: B018
                obj.reactant_gas
            ).mass  # raises on invalid formula
        except Exception as e:
            raise ValueError(
                f"Invalid reactant gas molecular formula: {obj.reactant_gas}. "
                f"Supported gases: {', '.join(diffusion_coef.sigmas.keys())}, or "
                f"other if manually inputting diffusion coefficient"
            ) from e
    if obj.carrier_gas not in viscosity_density.a:
        raise ValueError(
            f"Unsupported carrier gas. "
            f"Supported gases: {', '.join(viscosity_density.a.keys())}"
        )

    # Check physicality of injector dimensions
    if obj.injector_ID < 0 or obj.injector_OD < 0:
        raise ValueError("Injector ID and OD must be positive")
    elif obj.injector_ID > obj.FT_ID:
        raise ValueError("Injector ID cannot be larger than flow tube ID")
    elif obj.injector_OD > obj.FT_ID:
        raise ValueError("Injector OD cannot be larger than flow tube ID")
    elif obj.injector_ID > obj.injector_OD:
        raise ValueError("Injector ID cannot be larger than injector OD")
    elif obj.injector_ID == 0 or obj.injector_OD == 0:
        raise ValueError("Injector dimensions must be non-zero")

    # Check reactant concentration inputs
    if obj.reactant_conc < 0:
        raise ValueError("Reactant concentration must be non-negative")
    if obj.reactant_conc_type not in [
        "ppm",
        "ppb",
        "ng/min",
        "Pa",
        "hPa",
        "Torr",
        "bar",
        "mbar",
    ]:
        raise ValueError(
            "Unsupported reactant concentration type. "
            "Supported types: 'ppm', 'ppb', 'ng/min', 'Pa', 'hPa', 'Torr', 'bar', 'mbar'"
        )

    ### Check Reactor specific inputs ###
    if isinstance(obj, BoatReactor):
        # Check physicality of boat dimensions
        if (
            obj.boat_length < 0
            or obj.boat_liquid_width < 0
            or obj.boat_cross_section < 0
        ):
            raise ValueError("Boat dimensions must be positive")
        elif (
            obj.boat_liquid_width > obj.FT_ID
            or obj.boat_cross_section > np.pi * (obj.FT_ID / 2) ** 2
        ):
            raise ValueError(
                "Boat liquid width cannot be larger than the flow tube ID, and "
                "boat cross-sectional area cannot be larger than the flow tube "
                "cross-sectional area"
            )
        elif obj.boat_length > obj.FT_length:
            raise ValueError("Boat length cannot be larger than flow tube length")
        if obj.boat_perimeter is not None and obj.boat_perimeter < 0:
            raise ValueError("Boat perimeter must be positive")
    elif isinstance(obj, CoatedWallReactor):
        # Check physicality of insert dimensions
        if np.isnan(obj.insert_ID) != np.isnan(obj.insert_OD):
            raise ValueError(
                "Insert dimensions must all be specified or all be unspecified"
            )
        elif obj.insert_ID <= 0 or obj.insert_OD <= 0:
            raise ValueError("Insert ID and OD must be non-zero and positive")
        elif obj.insert_ID > obj.FT_ID or obj.insert_OD > obj.FT_ID:
            raise ValueError("Insert cannot be larger than flow tube ID")
        elif not np.isnan(obj.insert_ID) and obj.insert_ID > obj.insert_OD:
            raise ValueError("Insert ID cannot be larger than insert OD")


def validate_initialize(
    obj: AerosolFlowReactor | CoatedWallReactor | BoatReactor,
) -> None:
    """Validate inputs for the initialize method.

    Args:
        obj: An instance of a reactor class (AerosolFlowReactor,
            CoatedWallReactor, or BoatReactor).

    Returns:
        None. Raises errors if any validation checks fail.
    """
    from .aerosol_flow_reactor import AerosolFlowReactor

    # Check if flow rates are positive
    if obj.reactant_FR < 0 or obj.reactant_carrier_FR < 0 or obj.carrier_FR < 0:
        raise ValueError("Flow rates must be positive")

    # Check for non-zero flow
    if obj.reactant_FR <= 0:
        raise ValueError("Reactant flow rate must be positive and non-zero")
    if obj.reactant_carrier_FR < 0 or obj.carrier_FR < 0:
        raise ValueError("Flow rates must be positive or zero")

    # Check if the pressure units are supported
    if obj.P_units not in tools.P_CF:
        raise ValueError(
            f"Unsupported pressure units. "
            f"Supported units: {', '.join(tools.P_CF.keys())}"
        )
    elif obj.P < 0:
        raise ValueError("Pressure must be positive")

    # Check if the temperature & temperature gradients are valid numbers
    if obj.T < -273.15:
        raise ValueError("Temperature must be above absolute zero (-273.15 C)")
    if obj.radial_delta_T < 0:
        raise ValueError("Temperature gradients must be positive")

    if obj.axial_distance < 0:
        raise ValueError("Axial distance must be positive")

    # Check if axial distance is positive
    if obj.axial_distance < 0:
        raise ValueError("Axial distance must be positive")

    if isinstance(obj, AerosolFlowReactor):
        # Check physicality of aerosol inputs
        if obj.aerosol_distribution not in ["monodisperse", "lognormal"]:
            raise ValueError(
                "Unsupported aerosol distribution. "
                "Supported distributions: 'monodisperse', 'lognormal'"
            )
        if obj.aerosol_diameter <= 0:
            raise ValueError("Aerosol diameter must be positive")
        if obj.aerosol_diameter < 1:
            warn("Aerosol diameter is <1 nm. Verify that aerosol diameter is in nm.")
        if obj.aerosol_diameter > 1e5:
            warn(
                "Aerosol diameter is >10,000 nm. Verify that aerosol diameter is in nm."
            )
        if obj.aerosol_number_conc < 0:
            raise ValueError("Aerosol number concentration must be non-negative")
        if obj.aerosol_density <= 0:
            raise ValueError("Aerosol density must be positive")
        if obj.aerosol_density > 3 or obj.aerosol_density < 0.5:
            warn(
                "Aerosol density is outside the typical range of 0.5-3 g/cm^3. "
                "Verify that aerosol density is in g/cm^3."
            )
        if obj.aerosol_distribution == "lognormal":
            if np.isnan(obj.aerosol_sigma):
                raise ValueError(
                    "Aerosol geometric standard deviation must be specified for lognormal distribution"
                )
            if obj.aerosol_sigma <= 0:
                raise ValueError(
                    "Aerosol geometric standard deviation must be positive"
                )


def validate_reactant_uptake(
    obj,
    hypothetical_gamma,
    exposure_length,
    exposure_time=None,
    gamma_wall=None,
) -> np.ndarray | float:
    """Validate inputs for the reactant_uptake method.

    Args:
        obj: An instance of a reactor class (AerosolFlowReactor,
            CoatedWallReactor, or BoatReactor).
        hypothetical_gamma (ArrayLike or float): Hypothetical
            uptake coefficient to calculate diffusion correction
            factor.
        exposure_length (float): Length of the exposed surface in
            cm. Default is 1 cm.
        exposure_time (float): Time in minutes over which the
            surface is exposed to the reactant. Default is 10
            minutes (optional).
        gamma_wall (float): Uptake coefficient for the wall
            (optional).

    Returns:
        hypothetical_gamma (np.ndarray or float): Validated hypothetical uptake coefficient.
    """
    if not isinstance(hypothetical_gamma, (int, float)):
        try:
            hypothetical_gamma = np.asarray(hypothetical_gamma, dtype=np.float64)
        except Exception as e:
            raise TypeError(
                "Gamma input must be float or Array-like of float; "
                f"got {type(hypothetical_gamma)}"
            ) from e

        if hypothetical_gamma.ndim != 1:
            raise ValueError("Gamma input must be 1-dimensional.")

    # Verify that the exposure length is a positive number and that
    # it is less than the axial distance of the flow tube or insert
    if exposure_length <= 0:
        raise ValueError("Exposure length must be a positive number.")
    if exposure_length > obj.axial_distance:
        raise ValueError(
            "Exposure length must be less than the axial distance. Set "
            "object.axial_distance = ... with a larger axial_distance."
        )

    # Check exposure time
    if exposure_time is not None and exposure_time <= 0:
        raise ValueError("Exposure time must be a positive number.")

    # Check wall gamma
    if gamma_wall is not None and (gamma_wall < 0 or gamma_wall > 1):
        raise ValueError("Wall gamma must be between 0 and 1")

    # Check if hypothetical_gamma is between 0 and 1
    if np.min(hypothetical_gamma) < 0 or np.max(hypothetical_gamma) > 1:
        raise ValueError("Hypothetical gamma must be between 0 and 1")

    return hypothetical_gamma
