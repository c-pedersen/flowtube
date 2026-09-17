"""
flowtube - A Python package for transport and diffusion calculations in
cylindrical flow reactors.

This package provides tools and utilities for flow reactor analysis
including aerosol flow reactor (AFR), boat reactor, and coated wall
reactor (CWR), viscosity/density, and binary diffusion coefficients
calculations for atmospheric chemistry research.
"""

from typing import TYPE_CHECKING

__version__ = "1.4.0"
__author__ = "Corey Pedersen"
__email__ = "coreyped@gmail.com"

# Import main modules and classes for easy access
from . import diffusion_coef, flow_calc, kinetics, tools, viscosity_density
from .aerosol_flow_reactor import AerosolFlowReactor
from .boat_reactor import BoatReactor
from .coated_wall_reactor import CoatedWallReactor

# Explicit type hint for Pylance
if TYPE_CHECKING:
    from .aerosol_flow_reactor import AerosolFlowReactor as _AerosolFlowReactor

    AerosolFlowReactor: type[_AerosolFlowReactor]

    from .boat_reactor import BoatReactor as _BoatReactor

    BoatReactor: type[_BoatReactor]

    from .coated_wall_reactor import CoatedWallReactor as _CoatedWallReactor

    CoatedWallReactor: type[_CoatedWallReactor]


__all__ = [
    "AerosolFlowReactor",
    "BoatReactor",
    "CoatedWallReactor",
    "diffusion_coef",
    "flow_calc",
    "kinetics",
    "tools",
    "viscosity_density",
]
