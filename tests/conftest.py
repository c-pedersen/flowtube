# tests/conftest.py
import pytest
from flowtube.coated_wall_reactor import CoatedWallReactor
from flowtube.boat_reactor import BoatReactor

@pytest.fixture
def reactor_classes():
    return [CoatedWallReactor, BoatReactor]

@pytest.fixture
def make_constructor_kwargs():
    def _make(Reactor, **overrides):
        # Common ctor args
        base = dict(
            FT_ID=2.6,
            FT_length=100,
            injector_ID=0.60,
            injector_OD=1.20,
            reactant_gas="HCl",
            carrier_gas="N2",
            reactant_MR=1e-6,
        )
        # Class-specific ctor args
        if Reactor is CoatedWallReactor:
            base.update(dict(insert_ID=0.0, insert_length=0.0))
        elif Reactor is BoatReactor:
            base.update(dict(
                boat_width=2.0,
                boat_height=0.963,
                boat_length=53.8,
                boat_wall_thickness=0.11,
            ))
        base.update(overrides)
        return base
    return _make

@pytest.fixture
def make_init_kwargs():
    def _make(Reactor, **overrides):
        base = dict(
            reactant_FR=1.0, 
            reactant_carrier_FR=50.0,
            carrier_FR=5000.0,
            P=2000.0, 
            P_units="Pa",
            T=25.0,
            reactant_diffusion_rate=None,
            radial_delta_T=1.0,
            axial_delta_T=1.0,
            disp=False,
        )
        base.update(overrides)
        return base
    return _make

@pytest.fixture
def build_reactor(make_constructor_kwargs, make_init_kwargs):
    """
    Build helper: constructs Reactor with ctor kwargs, then calls initialize
    with init kwargs (unless call_initialize=False).
    """
    def _build(
            Reactor, 
            constructor_overrides=None,
            init_overrides=None,
            call_initialize=True,
        ):
        constructor_kwargs = make_constructor_kwargs(
            Reactor, 
            **(constructor_overrides or {}),
        )
        init_kwargs = make_init_kwargs(
            Reactor, 
            **(init_overrides or {}),
        )
        obj = Reactor(**constructor_kwargs)
        if call_initialize:
            obj.initialize(**init_kwargs)
        return obj, constructor_kwargs, init_kwargs
    return _build