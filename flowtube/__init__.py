"""
flowtube - A Python package for transport and diffusion calculations in
cylindrical flow reactors.

This package provides tools and utilities for flow reactor analysis
including coated wall reactor (CWR), boat reactor, viscosity/density,
and binary diffusion coefficients calculations for atmospheric chemistry
research.
"""

from typing import TYPE_CHECKING

__version__ = "1.3.1"
__author__ = "Corey Pedersen"
__email__ = "coreyped@gmail.com"

# Import main modules and classes for easy access
from . import diffusion_coef, flow_calc, kinetics, tools, viscosity_density
from .boat_reactor import BoatReactor
from .coated_wall_reactor import CoatedWallReactor

# Explicit type hint for Pylance
if TYPE_CHECKING:
    from .coated_wall_reactor import CoatedWallReactor as _CoatedWallReactor

    CoatedWallReactor: type[_CoatedWallReactor]
    from .boat_reactor import BoatReactor as _BoatReactor

    BoatReactor: type[_BoatReactor]

__all__ = [
    "BoatReactor",
    "CoatedWallReactor",
    "diffusion_coef",
    "flow_calc",
    "kinetics",
    "tools",
    "viscosity_density",
]
