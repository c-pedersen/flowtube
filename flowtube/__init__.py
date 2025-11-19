"""
flowtube - A Python package for transport and diffusion calculations in
cylindrical flow reactors.

This package provides tools and utilities for flow reactor analysis
including coated wall reactor (CWR), boat reactor, viscosity/density, 
and binary diffusion coefficients calculations for atmospheric chemistry 
research.
"""

from typing import TYPE_CHECKING

__version__ = "1.2.0"
__author__ = "Corey Pedersen"
__email__ = "coreyped@gmail.com"

# Import main modules for easy access
from . import tools
from . import viscosity_density
from . import diffusion_coef
from . import flow_calc
from . import kinetics

# Import the main classes for direct access
from .coated_wall_reactor import CoatedWallReactor
from .boat_reactor import BoatReactor

# Explicit type hint for Pylance
if TYPE_CHECKING:
    from .coated_wall_reactor import CoatedWallReactor as _CoatedWallReactor

    CoatedWallReactor: type[_CoatedWallReactor]
    from .boat_reactor import BoatReactor as _BoatReactor
    BoatReactor: type[_BoatReactor]

__all__ = [
    "CoatedWallReactor",
    "BoatReactor",
    "tools", 
    "viscosity_density",
    "diffusion_coef",
    "flow_calc",
    "kinetics",
]
