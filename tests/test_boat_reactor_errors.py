import pytest
import numpy as np
from flowtube import BoatReactor

""" Tests for errors specific to BoatReactor. """


# Boat-specific geometry
def test_boat_dimensions_must_fit_inside_tube(make_constructor_kwargs):
    kwargs = make_constructor_kwargs(BoatReactor, boat_width=3.0, FT_ID=2.54)
    with pytest.raises(ValueError, match=r"Boat.*width.*height.*larger"):
        BoatReactor(**kwargs)


def test_boat_wall_thickness_positive(make_constructor_kwargs):
    kwargs = make_constructor_kwargs(
        BoatReactor,
        boat_wall_thickness=-0.1,
        boat_width=-0.1,
        boat_height=-0.1,
    )
    with pytest.raises(ValueError, match=r"Boat dimensions must be positive"):
        BoatReactor(**kwargs)


def test_fitting_before_reactant_uptake(build_reactor):
    # init error: negative flow in initialize
    obj, _, _ = build_reactor(BoatReactor)
    with pytest.raises(RuntimeError, match=r"Must call reactant_uptake*"):
        obj.calculate_gamma(
            concentrations=np.array([0.01, 0.02, 0.03]),
            exposure=np.array([0.1, 0.2, 0.3]),
            exposure_units="s",
        )


def test_boat_uptake(build_reactor):
    boat, _, _ = build_reactor(BoatReactor)
    gamma = 1e-6
    boat.reactant_uptake(hypothetical_gamma=gamma)

    assert (
        boat.k * boat.FT_ID / boat.reactant_molec_velocity * boat.geometric_correction
        == pytest.approx(gamma)
    )
    assert 1 - np.exp(
        -gamma
        / boat.geometric_correction
        * boat.reactant_molec_velocity
        / boat.FT_ID
        * boat.residence_time
        / 4,
    ) == pytest.approx(boat.uptake)
