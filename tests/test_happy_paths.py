# tests/test_happy_paths.py
import pytest
import warnings
from flowtube.coated_wall_reactor import CoatedWallReactor
from flowtube.boat_reactor import BoatReactor

BOTH = [CoatedWallReactor, BoatReactor]


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_valid_constructor_and_init_no_errors(
    Reactor,
    build_reactor,
):
    # Passes if neither __init__ nor initialize raises
    build_reactor(Reactor)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
@pytest.mark.parametrize("P_units", ["Pa", "torr", "Torr", "mbar", "hPa", "bar"])
def test_supported_pressure_units_do_not_raise(
    Reactor,
    build_reactor,
    make_init_kwargs,
    P_units,
):
    warnings.filterwarnings("ignore")
    obj, _, _ = build_reactor(Reactor, call_initialize=False)
    ok = make_init_kwargs(Reactor, P_units=P_units)
    obj.initialize(**ok)  # Should not raise


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
@pytest.mark.parametrize("carrier", ["Ar", "He", "N2", "O2"])
def test_supported_carrier_gases_ctor(
    Reactor,
    make_constructor_kwargs,
    carrier,
):
    kwargs = make_constructor_kwargs(Reactor, carrier_gas=carrier)
    Reactor(**kwargs)  # Should not raise
