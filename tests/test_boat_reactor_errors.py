import pytest
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
