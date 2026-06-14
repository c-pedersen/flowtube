# tests/test_auto_reinit.py
"""
Tests to verify that modifying a constructor attribute on an already-initialized
CoatedWallReactor or BoatReactor produces the same result as building a fresh
reactor with that attribute set from the start.
"""

import pytest
import numpy as np
from flowtube.coated_wall_reactor import CoatedWallReactor
from flowtube.boat_reactor import BoatReactor


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def build_fresh(Reactor, make_constructor_kwargs, make_init_kwargs, **ctor_overrides):
    """Build and initialize a reactor from scratch with given overrides."""
    ctor = make_constructor_kwargs(Reactor, **ctor_overrides)
    obj = Reactor(**ctor)
    init = make_init_kwargs(Reactor)
    obj.initialize(**init)
    return obj


def build_and_mutate(Reactor, make_constructor_kwargs, make_init_kwargs, attr, value):
    """Build a default reactor, then mutate one constructor attribute."""
    ctor = make_constructor_kwargs(Reactor)
    obj = Reactor(**ctor)
    init = make_init_kwargs(Reactor)
    obj.initialize(**init)
    setattr(obj, attr, value)  # should trigger auto-reinit
    return obj


# ---------------------------------------------------------------------------
# Computed attributes to compare between fresh and mutated reactors
# ---------------------------------------------------------------------------

CWR_COMPUTED_ATTRS = [
    "total_FR",
    "FT_flow_velocity",
    "FT_residence_time",
    "Re_FT",
    "Pe_FT",
    "carrier_dynamic_viscosity",
    "carrier_density",
    "reactant_diffusion_rate",
    "reactant_molec_velocity",
    "reactant_mean_free_path",
    "N_eff_Shw_FT",
    "Kn_FT",
    "FT_conc",
    "FT_conc_molec",
]

BOAT_COMPUTED_ATTRS = [
    "total_FR",
    "flow_velocity",
    "residence_time",
    "Re",
    "carrier_dynamic_viscosity",
    "carrier_density",
    "reactant_diffusion_rate",
    "reactant_molec_velocity",
    "N_eff_Shw_FT",
    "Kn_FT",
    "FT_conc",
    "FT_conc_molec",
]

COMPUTED_ATTRS = {
    CoatedWallReactor: CWR_COMPUTED_ATTRS,
    BoatReactor: BOAT_COMPUTED_ATTRS,
}


def assert_reactors_equal(Reactor, fresh, mutated):
    for attr in COMPUTED_ATTRS[Reactor]:
        fresh_val = getattr(fresh, attr)
        mutated_val = getattr(mutated, attr)
        assert np.isclose(fresh_val, mutated_val), (
            f"Mismatch on '{attr}': fresh={fresh_val}, mutated={mutated_val}"
        )


# ---------------------------------------------------------------------------
# Tests: shared constructor attributes
# ---------------------------------------------------------------------------

BOTH = [CoatedWallReactor, BoatReactor]


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_FT_ID(Reactor, make_constructor_kwargs, make_init_kwargs):
    new_val = 3.0
    fresh = build_fresh(
        Reactor, make_constructor_kwargs, make_init_kwargs, FT_ID=new_val
    )
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "FT_ID", new_val
    )
    assert_reactors_equal(Reactor, fresh, mutated)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_FT_length(Reactor, make_constructor_kwargs, make_init_kwargs):
    new_val = 150.0
    fresh = build_fresh(
        Reactor, make_constructor_kwargs, make_init_kwargs, FT_length=new_val
    )
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "FT_length", new_val
    )
    assert_reactors_equal(Reactor, fresh, mutated)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_injector_ID(Reactor, make_constructor_kwargs, make_init_kwargs):
    new_val = 0.4
    fresh = build_fresh(
        Reactor, make_constructor_kwargs, make_init_kwargs, injector_ID=new_val
    )
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "injector_ID", new_val
    )
    assert_reactors_equal(Reactor, fresh, mutated)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_injector_OD(Reactor, make_constructor_kwargs, make_init_kwargs):
    new_val = 1.5
    fresh = build_fresh(
        Reactor, make_constructor_kwargs, make_init_kwargs, injector_OD=new_val
    )
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "injector_OD", new_val
    )
    assert_reactors_equal(Reactor, fresh, mutated)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_reactant_gas(Reactor, make_constructor_kwargs, make_init_kwargs):
    new_val = "Cl2"
    fresh = build_fresh(
        Reactor, make_constructor_kwargs, make_init_kwargs, reactant_gas=new_val
    )
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "reactant_gas", new_val
    )
    assert_reactors_equal(Reactor, fresh, mutated)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_carrier_gas(Reactor, make_constructor_kwargs, make_init_kwargs):
    new_val = "Ar"
    fresh = build_fresh(
        Reactor, make_constructor_kwargs, make_init_kwargs, carrier_gas=new_val
    )
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "carrier_gas", new_val
    )
    assert_reactors_equal(Reactor, fresh, mutated)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_reactant_conc(Reactor, make_constructor_kwargs, make_init_kwargs):
    new_val = 100.0
    fresh = build_fresh(
        Reactor, make_constructor_kwargs, make_init_kwargs, reactant_conc=new_val
    )
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "reactant_conc", new_val
    )
    assert_reactors_equal(Reactor, fresh, mutated)


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_change_reactant_conc_type(Reactor, make_constructor_kwargs, make_init_kwargs):
    # Switch to ppb and set a matching concentration (30 ppm == 30000 ppb)
    fresh = build_fresh(
        Reactor,
        make_constructor_kwargs,
        make_init_kwargs,
        reactant_conc_type="ppb",
        reactant_conc=30000.0,
    )
    mutated = build_and_mutate(
        Reactor,
        make_constructor_kwargs,
        make_init_kwargs,
        "reactant_conc_type",
        "ppb",
    )
    mutated.reactant_conc = 30000.0
    assert_reactors_equal(Reactor, fresh, mutated)


# ---------------------------------------------------------------------------
# Tests: CWR-only constructor attributes (insert)
# ---------------------------------------------------------------------------


def test_cwr_change_insert_ID(make_constructor_kwargs, make_init_kwargs):
    new_val = 1.0
    ctor_overrides = {"insert_ID": new_val, "insert_OD": 1.5, "insert_length": 50.0}
    fresh = build_fresh(
        CoatedWallReactor, make_constructor_kwargs, make_init_kwargs, **ctor_overrides
    )
    # Build with insert, then mutate insert_ID
    ctor = make_constructor_kwargs(CoatedWallReactor, **ctor_overrides)
    obj = CoatedWallReactor(**ctor)
    obj.initialize(**make_init_kwargs(CoatedWallReactor))
    obj.insert_ID = new_val
    assert_reactors_equal(CoatedWallReactor, fresh, obj)


def test_cwr_change_insert_length(make_constructor_kwargs, make_init_kwargs):
    new_val = 30.0
    ctor_overrides = {"insert_ID": 1.0, "insert_OD": 1.5, "insert_length": new_val}
    fresh = build_fresh(
        CoatedWallReactor, make_constructor_kwargs, make_init_kwargs, **ctor_overrides
    )
    ctor = make_constructor_kwargs(CoatedWallReactor, **ctor_overrides)
    obj = CoatedWallReactor(**ctor)
    obj.initialize(**make_init_kwargs(CoatedWallReactor))
    obj.insert_length = new_val
    assert_reactors_equal(CoatedWallReactor, fresh, obj)


# ---------------------------------------------------------------------------
# Tests: BoatReactor-only constructor attributes
# ---------------------------------------------------------------------------


def test_boat_change_liquid_width(make_constructor_kwargs, make_init_kwargs):
    new_val = 1.5
    fresh = build_fresh(
        BoatReactor,
        make_constructor_kwargs,
        make_init_kwargs,
        boat_liquid_width=new_val,
    )
    mutated = build_and_mutate(
        BoatReactor,
        make_constructor_kwargs,
        make_init_kwargs,
        "boat_liquid_width",
        new_val,
    )
    assert_reactors_equal(BoatReactor, fresh, mutated)


def test_boat_change_boat_length(make_constructor_kwargs, make_init_kwargs):
    new_val = 30.0
    fresh = build_fresh(
        BoatReactor, make_constructor_kwargs, make_init_kwargs, boat_length=new_val
    )
    mutated = build_and_mutate(
        BoatReactor, make_constructor_kwargs, make_init_kwargs, "boat_length", new_val
    )
    assert_reactors_equal(BoatReactor, fresh, mutated)


def test_boat_change_cross_section(make_constructor_kwargs, make_init_kwargs):
    new_val = 2.0
    fresh = build_fresh(
        BoatReactor,
        make_constructor_kwargs,
        make_init_kwargs,
        boat_cross_section=new_val,
    )
    mutated = build_and_mutate(
        BoatReactor,
        make_constructor_kwargs,
        make_init_kwargs,
        "boat_cross_section",
        new_val,
    )
    assert_reactors_equal(BoatReactor, fresh, mutated)


# ---------------------------------------------------------------------------
# Tests: validation still fires on bad values
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_invalid_injector_ID_raises(Reactor, make_constructor_kwargs, make_init_kwargs):
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "injector_OD", 1.20
    )
    with pytest.raises(
        ValueError, match=r"Injector ID cannot be larger than injector OD"
    ):
        mutated.injector_ID = 1.5


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_invalid_carrier_gas_raises(Reactor, make_constructor_kwargs, make_init_kwargs):
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "FT_ID", 2.6
    )
    with pytest.raises(ValueError, match=r"Unsupported carrier gas"):
        mutated.carrier_gas = "Xe"


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_invalid_reactant_gas_raises(
    Reactor, make_constructor_kwargs, make_init_kwargs
):
    mutated = build_and_mutate(
        Reactor, make_constructor_kwargs, make_init_kwargs, "FT_ID", 2.6
    )
    with pytest.raises(ValueError, match=r"Invalid reactant gas molecular formula"):
        mutated.reactant_gas = "NotAGas!!"


def test_boat_invalid_dimensions_raises(make_constructor_kwargs, make_init_kwargs):
    mutated = build_and_mutate(
        BoatReactor, make_constructor_kwargs, make_init_kwargs, "FT_ID", 2.6
    )
    with pytest.raises(ValueError, match=r"Boat.*width.*larger"):
        mutated.boat_liquid_width = 999.0


# ---------------------------------------------------------------------------
# Tests: warning when initialize() has not been called yet
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
def test_no_reinit_before_initialize(Reactor, make_constructor_kwargs, capsys):
    ctor = make_constructor_kwargs(Reactor)
    obj = Reactor(**ctor)

    with pytest.raises(RuntimeError, match=r"initialize\(\) has not been called yet"):
        obj.FT_ID = 3.0  # initialize() not called yet
