# tests/test_common_errors.py
import pytest
from flowtube import CoatedWallReactor
from flowtube import BoatReactor
import numpy as np

""" Tests for common errors across Reactor classes. """

BOTH = [CoatedWallReactor, BoatReactor]


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_valid_constructor_and_init_no_errors(Reactor, build_reactor):
    # Passes if neither __init__ nor initialize raises
    build_reactor(Reactor)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_injector_dimensions(Reactor, make_constructor_kwargs):
    # Test injector ID > OD
    kwargs = make_constructor_kwargs(Reactor, injector_ID=0.5, injector_OD=0.4)
    with pytest.raises(
        ValueError, match=r"Injector ID cannot be larger than injector OD"
    ):
        Reactor(**kwargs)

    # Test non-zero injector dimensions
    kwargs = make_constructor_kwargs(Reactor, injector_ID=0.0)
    with pytest.raises(ValueError, match=r"Injector dimensions must be non-zero"):
        Reactor(**kwargs)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_unsupported_carrier_gas(Reactor, make_constructor_kwargs):
    kwargs = make_constructor_kwargs(Reactor, carrier_gas="Xe")
    with pytest.raises(ValueError, match=r"Unsupported carrier gas"):
        Reactor(**kwargs)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_negative_flows_initialize(Reactor, build_reactor, make_init_kwargs):
    # init error: negative flow in initialize
    obj, _, _ = build_reactor(Reactor, call_initialize=False)
    bad_init = make_init_kwargs(Reactor, reactant_FR=-1.0)
    with pytest.raises(ValueError, match=r"Flow rates.*positive"):
        obj.initialize(**bad_init)

    bad_init = make_init_kwargs(Reactor, carrier_FR=-1.0)
    with pytest.raises(ValueError, match=r"Flow rates.*positive"):
        obj.initialize(**bad_init)

    # init error: zero reactant flow in initialize
    obj, _, _ = build_reactor(Reactor, call_initialize=False)
    bad_init = make_init_kwargs(Reactor, reactant_FR=0.0)
    with pytest.raises(ValueError, match=r"flow rate.*non-zero"):
        obj.initialize(**bad_init)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_unsupported_pressure_units_initialize(
    Reactor,
    build_reactor,
    make_init_kwargs,
):
    obj, _, _ = build_reactor(Reactor, call_initialize=False)
    bad_init = make_init_kwargs(Reactor, P_units="atmz")
    with pytest.raises(ValueError, match=r"Unsupported pressure units"):
        obj.initialize(**bad_init)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_manual_reactant_diffusion_coef(
    Reactor,
    build_reactor,
    make_init_kwargs,
    make_constructor_kwargs,
):
    # Invalid gas but valid formula
    kwargs = make_constructor_kwargs(Reactor, reactant_gas="O8")
    with pytest.raises(
        ValueError,
        match=r"Must input reactant diffusion rate",
    ):
        build_reactor(Reactor, constructor_overrides=kwargs)

    # Invalid gas and invalid formula
    kwargs = make_constructor_kwargs(Reactor, reactant_gas="InvalidGas")
    with pytest.raises(
        ValueError,
        match=r"Invalid reactant gas molecular formula",
    ):
        build_reactor(Reactor, constructor_overrides=kwargs)

    # Valid gas but manual diffusion coefficient
    kwargs = make_constructor_kwargs(Reactor, reactant_gas="O8")
    init_kwargs = make_init_kwargs(Reactor, reactant_diffusion_rate=1.0)
    # Passes if neither __init__ nor initialize raises
    build_reactor(
        Reactor,
        constructor_overrides=kwargs,
        init_overrides=init_kwargs,
    )


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_temperature_below_physical_limit(Reactor, build_reactor, make_init_kwargs):
    obj, _, _ = build_reactor(Reactor, call_initialize=False)
    bad_init = make_init_kwargs(Reactor, T=-273.16)
    with pytest.raises(ValueError, match=r"temperature|absolute zero|below"):
        obj.initialize(**bad_init)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_mixing_ratio_bounds(Reactor, make_constructor_kwargs, build_reactor):
    # >1 is invalid
    kwargs = make_constructor_kwargs(Reactor, reactant_MR=1.01)
    with pytest.raises(ValueError, match=r"(mixing|MR).*0.*1"):
        build_reactor(Reactor, constructor_overrides=kwargs)

    # <0 is invalid
    kwargs = make_constructor_kwargs(Reactor, reactant_MR=-0.01)
    with pytest.raises(ValueError, match=r"(mixing|MR).*0.*1"):
        build_reactor(Reactor, constructor_overrides=kwargs)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_fitting(Reactor, build_reactor):
    # init error: negative flow in initialize
    obj, _, _ = build_reactor(Reactor, call_initialize=False)
    with pytest.raises(RuntimeError, match=r"Must call*"):
        obj.calculate_gamma(
            concentrations=np.array([0.01, 0.02, 0.03]),
            exposure=np.array([0.1, 0.2, 0.3]),
            exposure_units="s",
        )


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_fitting_invalid_exposure_units(Reactor, build_reactor):
    obj, _, _ = build_reactor(Reactor)
    obj.reactant_uptake(hypothetical_gamma=1e-7, disp=False)
    with pytest.raises(ValueError, match=r"Unsupported exposure units"):
        obj.calculate_gamma(
            concentrations=np.array([0.01, 0.02, 0.03]),
            exposure=np.array([0.1, 0.2, 0.3]),
            exposure_units="invalid_unit",
        )


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_fitting_negative_concentrations(Reactor, build_reactor):
    obj, _, _ = build_reactor(Reactor)
    obj.reactant_uptake(hypothetical_gamma=1e-7, disp=False)
    with pytest.raises(ValueError, match=r"Concentrations must be non-negative"):
        obj.calculate_gamma(
            concentrations=np.array([-0.01, 0.02, 0.03]),
            exposure=np.array([0.1, 0.2, 0.3]),
            exposure_units="s",
        )


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_fitting_non_numpy_inputs(Reactor, build_reactor):
    obj, _, _ = build_reactor(Reactor)
    obj.reactant_uptake(hypothetical_gamma=1e-7, disp=False)
    with pytest.raises(TypeError, match=r"Concentrations must be a numpy array"):
        obj.calculate_gamma(
            concentrations=[0.01, 0.02, 0.03],
            exposure=np.array([0.1, 0.2, 0.3]),
            exposure_units="s",
        )
    with pytest.raises(TypeError, match=r"Exposure must be a numpy array"):
        obj.calculate_gamma(
            concentrations=np.array([0.01, 0.02, 0.03]),
            exposure=[0.1, 0.2, 0.3],
            exposure_units="s",
        )
