"""Type stubs for flowtube package."""

from typing import Any
from . import tools as tools
from . import viscosity_density as viscosity_density
from . import diffusion_coef as diffusion_coef
from . import flow_calc as flow_calc

class CoatedWallReactor:
    """Handles calculations relevant to flow rate, flow diagnostics, transport, and uptake for a coated wall reactor."""
    
    FT_ID: float
    FT_length: float
    injector_ID: float
    injector_OD: float
    reactant_gas: str
    carrier_gas: str
    reactant_MR: float
    insert_ID: float
    insert_length: float
    
    def __init__(
        self,
        FT_ID: float,
        FT_length: float,
        injector_ID: float,
        injector_OD: float,
        reactant_gas: str,
        carrier_gas: str,
        reactant_MR: float,
        insert_ID: float = ...,
        insert_length: float = ...,
    ) -> None: ...
    
    def initialize(
        self,
        reactant_FR: float,
        reactant_carrier_FR: float,
        carrier_FR: float,
        P: float,
        P_units: str,
        T: float,
        radial_delta_T: float = ...,
        axial_delta_T: float = ...,
        disp: bool = ...,
    ) -> None: ...
    
    def flows(
        self,
        reactant_FR: float,
        reactant_carrier_FR: float,
        carrier_FR: float,
        disp: bool = ...,
    ) -> None: ...
    
    def carrier_flow(
        self,
        radial_delta_T: float = ...,
        axial_delta_T: float = ...,
        disp: bool = ...,
    ) -> None: ...
    
    def reactant_diffusion(self, disp: bool = ...) -> None: ...
    
    def reactant_uptake(
        self,
        hypothetical_gamma: Any,
        gamma_wall: float = ...,
        disp: bool = ...,
    ) -> tuple[Any, Any]: ...

__version__: str
__author__: str
__email__: str
__all__: list[str]