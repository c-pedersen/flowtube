"""
flowTube - A Python package for transport and diffusion calculations in cylindrical flow reactors.

This package provides tools and utilities for flow reactor analysis including
coated wall reactor (CWR) calculations, viscosity/density calculations, 
and binary diffusion coefficients for atmospheric chemistry research.
"""

__version__ = "0.1.0"
__author__ = "Corey Pedersen"
__email__ = "coreyped@gmail.com"

# Import main modules for easy access
from . import tools
from . import viscosity_density
from . import diffusion_coef

# Import the main FT class for direct access
from .flowTube import FT

# Define what gets imported with "from flowTube import *"
__all__ = [
    "FT",
    "tools", 
    "viscosity_density",
    "diffusion_coef",
]