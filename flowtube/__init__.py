"""
flowtube - A Python package for transport and diffusion calculations in
cylindrical flow reactors.

This package provides tools and utilities for flow reactor analysis
including coated wall reactor (CWR) calculations, viscosity/density
calculations, and binary diffusion coefficients for atmospheric
chemistry research.
"""

from typing import TYPE_CHECKING

__version__ = "1.0.0"
__author__ = "Corey Pedersen"
__email__ = "coreyped@gmail.com"

# Import main modules for easy access
from . import tools
from . import viscosity_density
from . import diffusion_coef
from . import flow_calc

# Import the main CWR class for direct access
from .coated_wall_reactor import CoatedWallReactor

# Explicit type hint for Pylance
if TYPE_CHECKING:
    from .coated_wall_reactor import CoatedWallReactor as _CoatedWallReactor

    CoatedWallReactor: type[_CoatedWallReactor]

__all__ = [
    "CoatedWallReactor",
    "tools",
    "viscosity_density",
    "diffusion_coef",
    "flow_calc",
]
