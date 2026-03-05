# tests/test_diffusion_corrected_uptake.py
import pytest
import numpy as np
from flowtube.coated_wall_reactor import CoatedWallReactor
from flowtube.kinetics import correction_factor_from_effective_gamma


def test_diffusion_corrected_uptake_no_insert(build_reactor):
    """
    Without an insert the method uses N_eff_Shw_FT and Kn_FT to correct
    the observed effective gamma.  The true gamma must satisfy:

        gamma = effective_gamma / Cg
        Cg    = 1 - (effective_gamma * 3 / (2 * N_eff_Shw_FT * Kn_FT))
    """
    obj, _, _ = build_reactor(CoatedWallReactor)

    effective_gamma = 1e-4
    result = obj.diffusion_corrected_uptake_coefficient(effective_gamma)

    # Manual calculation using the reactor's own attributes
    Cg_expected = correction_factor_from_effective_gamma(
        N_eff_Shw=obj.N_eff_Shw_FT,
        Kn=obj.Kn_FT,
        effective_gamma=effective_gamma,
    )
    gamma_expected = effective_gamma / Cg_expected

    assert np.isclose(result, gamma_expected)


def test_diffusion_corrected_uptake_with_insert(build_reactor):
    """
    When an insert is present the reactor exposes both flow-tube and insert
    correction parameters. diffusion_corrected_uptake_coefficient must use
    the insert parameters, not the flowtube ones.
    """
    obj, _, _ = build_reactor(
        CoatedWallReactor,
        constructor_overrides={"insert_ID": 1.0, "insert_length": 50.0},
    )

    # Reactor must expose both sets of correction attributes
    assert hasattr(obj, "N_eff_Shw_insert")
    assert hasattr(obj, "Kn_insert")

    effective_gamma = 1e-4
    result = obj.diffusion_corrected_uptake_coefficient(effective_gamma)

    # Result must differ from the FT-based formula …
    Cg_ft = correction_factor_from_effective_gamma(
        N_eff_Shw=obj.N_eff_Shw_FT,
        Kn=obj.Kn_FT,
        effective_gamma=effective_gamma,
    )
    assert not np.isclose(result, effective_gamma / Cg_ft)

    # … and match what insert-based parameters would give
    Cg_insert = correction_factor_from_effective_gamma(
        N_eff_Shw=obj.N_eff_Shw_insert,
        Kn=obj.Kn_insert,
        effective_gamma=effective_gamma,
    )
    assert np.isclose(result, effective_gamma / Cg_insert)


def test_diffusion_correction_increases_gamma(build_reactor):
    """
    In the typical diffusion-limited regime the correction factor Cg < 1,
    so the true gamma must be *larger* than the observed effective gamma.
    """
    obj, _, _ = build_reactor(CoatedWallReactor)

    effective_gamma = 1e-4
    result = obj.diffusion_corrected_uptake_coefficient(effective_gamma)

    assert result > effective_gamma


def test_diffusion_corrected_uptake_warns_unphysical(build_reactor):
    """
    When the computed gamma is outside [0, 1] a UserWarning is emitted.
    This happens when effective_gamma approaches
    2 * N_eff_Shw_FT * Kn_FT / 3 (denominator of Cg approaches zero).
    """
    obj, _, _ = build_reactor(CoatedWallReactor)

    # Choose an effective_gamma large enough to push the denominator past
    # zero, making the corrected gamma negative / unphysical.
    threshold = 2 * obj.N_eff_Shw_FT * obj.Kn_FT / 3
    unphysical_gamma = threshold * 1.1

    with pytest.warns(UserWarning, match="unphysical"):
        obj.diffusion_corrected_uptake_coefficient(unphysical_gamma)
