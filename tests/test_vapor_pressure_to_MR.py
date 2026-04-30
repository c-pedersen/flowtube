# tests/test_vapor_pressure_to_MR.py
"""Unit tests for tools.vapor_pressure_to_MR and its use in both reactor classes."""

import pytest
import warnings
from flowtube.tools import vapor_pressure_to_MR
from flowtube import CoatedWallReactor, BoatReactor

BOTH = [CoatedWallReactor, BoatReactor]


# ---------------------------------------------------------------------------
# Direct tests of vapor_pressure_to_MR
# ---------------------------------------------------------------------------


class TestVaporPressureToMR:
    """Direct tests of the vapor_pressure_to_MR helper."""

    def test_mixing_ratio_uses_system_pressure(self):
        """MR = vapor_pressure_Pa / system_pressure_Pa."""
        result = vapor_pressure_to_MR(
            vapor_pressure=100.0,
            P_units="Pa",
            system_pressure=10000.0,
            P_units_system="Pa",
        )
        assert abs(result - 0.01) < 1e-12

    def test_mixing_ratio_system_pressure_dependence(self):
        """Higher system pressure yields a lower mixing ratio for the same vapor pressure."""
        mr_low = vapor_pressure_to_MR(
            vapor_pressure=50.0,
            P_units="Pa",
            system_pressure=5000.0,
            P_units_system="Pa",
        )
        mr_high = vapor_pressure_to_MR(
            vapor_pressure=50.0,
            P_units="Pa",
            system_pressure=10000.0,
            P_units_system="Pa",
        )
        assert mr_low > mr_high

    @pytest.mark.parametrize(
        "vp, vp_units, sp, sp_units",
        [
            (1.0, "Torr", 760.0, "Torr"),
            (100.0, "Pa", 1.0, "bar"),
            (0.5, "mbar", 1013.25, "mbar"),
            (0.5, "hPa", 1013.25, "hPa"),
        ],
    )
    def test_unit_conversion_consistent(self, vp, vp_units, sp, sp_units):
        """Results are the same regardless of the pressure unit used."""
        from flowtube.tools import P_in_Pa

        expected = P_in_Pa(vp, vp_units) / P_in_Pa(sp, sp_units)
        result = vapor_pressure_to_MR(vp, vp_units, sp, sp_units)
        assert abs(result - expected) < 1e-12

    # --- Validation: vapor_pressure > system_pressure ---

    def test_vapor_exceeds_system_raises(self):
        with pytest.raises(ValueError, match=r"cannot exceed"):
            vapor_pressure_to_MR(
                vapor_pressure=200.0,
                P_units="Pa",
                system_pressure=100.0,
                P_units_system="Pa",
            )

    # --- Validation: system_pressure <= 0 raises BEFORE comparison ---

    def test_zero_system_pressure_raises_positivity_error(self):
        """system_pressure = 0 must raise the positivity error, not the comparison error."""
        with pytest.raises(ValueError, match=r"[Ss]ystem pressure must be positive"):
            vapor_pressure_to_MR(
                vapor_pressure=1.0,
                P_units="Pa",
                system_pressure=0.0,
                P_units_system="Pa",
            )

    def test_negative_system_pressure_raises_positivity_error(self):
        """Negative system_pressure must raise the positivity error, not comparison error."""
        with pytest.raises(ValueError, match=r"[Ss]ystem pressure must be positive"):
            vapor_pressure_to_MR(
                vapor_pressure=1.0,
                P_units="Pa",
                system_pressure=-100.0,
                P_units_system="Pa",
            )

    # --- Validation: vapor_pressure < 0 ---

    def test_negative_vapor_pressure_raises(self):
        with pytest.raises(ValueError, match=r"[Vv]apor pressure must be non-negative"):
            vapor_pressure_to_MR(
                vapor_pressure=-1.0,
                P_units="Pa",
                system_pressure=1000.0,
                P_units_system="Pa",
            )

    # --- Warning: vapor_pressure > 1 % of system_pressure ---

    def test_high_vapor_pressure_warns(self):
        """Vapor pressure > 1 % of system pressure should trigger a UserWarning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            vapor_pressure_to_MR(
                vapor_pressure=200.0,
                P_units="Pa",
                system_pressure=1000.0,
                P_units_system="Pa",
            )
        assert any(issubclass(w.category, UserWarning) for w in caught)

    def test_low_vapor_pressure_no_warn(self):
        """Vapor pressure <= 1 % of system pressure should not trigger a UserWarning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            vapor_pressure_to_MR(
                vapor_pressure=10.0,
                P_units="Pa",
                system_pressure=1000.0,
                P_units_system="Pa",
            )
        user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user_warnings) == 0

    def test_zero_vapor_pressure_returns_zero(self):
        """Zero vapor pressure is valid and should return 0.0."""
        result = vapor_pressure_to_MR(
            vapor_pressure=0.0,
            P_units="Pa",
            system_pressure=1000.0,
            P_units_system="Pa",
        )
        assert result == 0.0


# ---------------------------------------------------------------------------
# Integration tests: vapor_pressure_to_MR through both reactor interfaces
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("Reactor", BOTH, ids=["CoatedWall", "Boat"])
class TestVaporPressureReactorIntegration:
    """Verify vapor_pressure_to_MR behavior when used via reactor initialize()."""

    def test_vapor_pressure_conc_type_succeeds(
        self, Reactor, make_constructor_kwargs, make_init_kwargs
    ):
        """Reactor initializes without error when reactant_conc_type='Pa' and vapor < system."""
        ctor = make_constructor_kwargs(Reactor, reactant_conc_type="Pa", reactant_conc=10.0)
        obj = Reactor(**ctor)
        init = make_init_kwargs(Reactor, P=2000.0, P_units="Pa")
        obj.initialize(**init)  # should not raise
        assert obj.reactant_MR == pytest.approx(10.0 / 2000.0)

    def test_vapor_pressure_exceeds_system_raises_via_reactor(
        self, Reactor, make_constructor_kwargs, make_init_kwargs
    ):
        """Reactor.initialize raises when vapor pressure > system pressure."""
        ctor = make_constructor_kwargs(Reactor, reactant_conc_type="Pa", reactant_conc=5000.0)
        obj = Reactor(**ctor)
        init = make_init_kwargs(Reactor, P=2000.0, P_units="Pa")
        with pytest.raises(ValueError, match=r"cannot exceed"):
            obj.initialize(**init)

    @pytest.mark.parametrize(
        "vp_units",
        ["Pa", "Torr", "bar", "mbar"],
    )
    def test_all_supported_vapor_pressure_units(
        self, Reactor, vp_units, make_constructor_kwargs, make_init_kwargs
    ):
        """Reactor.initialize works for every supported vapor-pressure unit."""
        from flowtube.tools import P_in_Pa

        warnings.filterwarnings("ignore")
        ctor = make_constructor_kwargs(Reactor, reactant_conc=0.001, reactant_conc_type=vp_units)
        obj = Reactor(**ctor)
        init = make_init_kwargs(Reactor, P=101325.0, P_units="Pa")
        obj.initialize(**init)  # should not raise
        expected_mr = P_in_Pa(0.001, vp_units) / 101325.0
        assert obj.reactant_MR == pytest.approx(expected_mr)
