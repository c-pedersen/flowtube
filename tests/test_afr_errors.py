"""Validation tests for aerosol-specific initialization inputs.

Citations:
Seinfeld, J.H. and Pandis, S.N. (2016) Atmospheric chemistry and
physics: from air pollution to climate change. Third edition. Hoboken,
New Jersey: John Wiley & Sons.

Wagner, C., Hanisch, F., Holmes, N., De Coninck, H., Schuster, G.,
Crowley, J.N., 2008. The interaction of N2 O5 with mineral dust: aerosol
flow tube and Knudsen reactor studies. Atmos. Chem. Phys. 8, 91–109.
https://doi.org/10.5194/acp-8-91-2008
"""

import warnings

import numpy as np
import pytest

from flowtube.aerosol_flow_reactor import AerosolFlowReactor


@pytest.fixture
def initialize_aerosol(make_constructor_kwargs, make_init_kwargs):
    def _initialize(**overrides):
        kwargs = {
            "aerosol_distribution": "monodisperse",
            "aerosol_diameter": 100.0,
            "aerosol_number_conc": 1e4,
            "aerosol_density": 1.0,
        }
        kwargs.update(overrides)
        reactor = AerosolFlowReactor(**make_constructor_kwargs(AerosolFlowReactor))
        init_kwargs = make_init_kwargs(AerosolFlowReactor, **kwargs)
        # Exercise initialize's optional defaults, not the shared fixture's sigma.
        if "aerosol_sigma" not in overrides:
            init_kwargs.pop("aerosol_sigma", None)
        reactor.initialize(**init_kwargs)
        return reactor

    return _initialize


@pytest.mark.parametrize("distribution", ["normal", "", "Monodisperse", None])
def test_unsupported_aerosol_distribution(initialize_aerosol, distribution):
    with pytest.raises(ValueError, match="Unsupported aerosol distribution"):
        initialize_aerosol(aerosol_distribution=distribution)


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("aerosol_diameter", -1.0, "Aerosol diameter must be positive"),
        ("aerosol_diameter", 0.0, "Aerosol diameter must be positive"),
        (
            "aerosol_number_conc",
            -1.0,
            "Aerosol number concentration must be non-negative",
        ),
        ("aerosol_density", -1.0, "Aerosol density must be positive"),
        ("aerosol_density", 0.0, "Aerosol density must be positive"),
        (
            "aerosol_surface_area",
            -1.0,
            "Aerosol surface area must be non-negative",
        ),
    ],
)
def test_invalid_aerosol_numeric_inputs(initialize_aerosol, field, value, message):
    with pytest.raises(ValueError, match=message):
        initialize_aerosol(**{field: value})


@pytest.mark.parametrize("sigma", [0.0, -1.0])
def test_lognormal_sigma_must_be_positive(initialize_aerosol, sigma):
    with pytest.raises(
        ValueError, match="Aerosol geometric standard deviation must be positive"
    ):
        initialize_aerosol(aerosol_distribution="lognormal", aerosol_sigma=sigma)


@pytest.mark.parametrize("overrides", [{}, {"aerosol_sigma": np.nan}])
def test_lognormal_sigma_is_required(initialize_aerosol, overrides):
    with pytest.raises(
        ValueError,
        match="geometric standard deviation must be specified for lognormal",
    ):
        initialize_aerosol(aerosol_distribution="lognormal", **overrides)


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("aerosol_diameter", 0.5, "Aerosol diameter is <1 nm"),
        ("aerosol_diameter", 100001.0, "Aerosol diameter is >"),
        ("aerosol_density", 0.49, "Aerosol density is outside the typical range"),
        ("aerosol_density", 3.01, "Aerosol density is outside the typical range"),
    ],
)
def test_unusual_aerosol_inputs_warn(initialize_aerosol, field, value, message):
    with pytest.warns(UserWarning, match=message):
        initialize_aerosol(**{field: value})


@pytest.mark.parametrize(
    "overrides",
    [
        {},
        {"aerosol_distribution": "lognormal", "aerosol_sigma": 1.5},
        {"aerosol_diameter": 1.0},
        {"aerosol_diameter": 1e5},
        {"aerosol_number_conc": 0.0},
        {"aerosol_density": 0.5},
        {"aerosol_density": 2.0},
        {"aerosol_surface_area": 0.0},
        {"aerosol_surface_area": 1e-5},
        {"aerosol_distribution": "monodisperse", "aerosol_sigma": np.nan},
    ],
)
def test_valid_aerosol_inputs(initialize_aerosol, overrides):
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message="Aerosol", category=UserWarning)
        reactor = initialize_aerosol(**overrides)

    for field, value in overrides.items():
        actual = getattr(reactor, field)
        if isinstance(value, float) and np.isnan(value):
            assert np.isnan(actual)
        else:
            assert actual == value


# From Table 9.3 of Seinfeld and Pandis, 2016.
# Diameters are in µm, and slip corrections are dimensionless.
@pytest.mark.parametrize(
    "diameter, slip_correction",
    [
        (0.001, 216),
        (0.002, 108),
        (0.005, 43.6),
        (0.01, 22.2),
        (0.02, 11.4),
        (0.05, 4.95),
        (0.1, 2.85),
        (0.2, 1.865),
        (0.5, 1.326),
        (1.0, 1.164),
        (2.0, 1.082),
        (5.0, 1.032),
        (10.0, 1.016),
        (20.0, 1.008),
        (50.0, 1.003),
        (100.0, 1.0016),
    ],
)
def test_aerosol_slip_correction_matches_seinfeld_pandis(
    initialize_aerosol, diameter, slip_correction
):
    """
    Test that the aerosol slip correction matches the values in Table 9.3
    of Seinfeld and Pandis, 2016.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        afr = initialize_aerosol(
            aerosol_distribution="monodisperse",
            aerosol_diameter=diameter * 1e3,
            T=25,
            P=101325.0,
            P_units="Pa",
            aerosol_density=1.0,
            carrier_FR=0,
            reactant_FR=0.1,
            reactant_carrier_FR=0,
        )

    assert np.isclose(afr.slip_correction, slip_correction, rtol=0.05)


# from Table 9.5 of Seinfeld and Pandis, 2016. Diamters are in µm,
# diffusivities in cm2 s-1, thermal velocities in cm s-1, and mean free
# paths in µm.
@pytest.mark.parametrize(
    "diameter, diffusivity, thermal_velocity, mean_free_path",
    [
        (0.002, 1.28e-2, 4965, 6.59e-2),
        (0.004, 3.23e-3, 1760, 4.68e-2),
        (0.01, 5.24e-4, 444, 3.00e-2),
        (0.02, 1.30e-4, 157, 2.20e-2),
        (0.04, 3.59e-5, 55.5, 1.64e-2),
        (0.1, 6.82e-6, 14.0, 1.24e-2),
        (0.2, 2.21e-6, 4.96, 1.13e-2),
        (0.4, 8.32e-7, 1.76, 1.21e-2),
        (1.0, 2.74e-7, 0.444, 1.53e-2),
        (2.0, 1.27e-7, 0.157, 2.06e-2),
        (4.0, 6.1e-8, 5.55e-2, 2.8e-2),
        (10.0, 2.38e-8, 1.40e-2, 4.32e-2),
    ],
)
def test_aerosol_properties_match_seinfeld_pandis(
    initialize_aerosol, diameter, diffusivity, thermal_velocity, mean_free_path
):
    """
    Test that the aerosol diffusivity, thermal velocity, and mean free
    path match the values in Table 9.5 of Seinfeld and Pandis, 2016.

    Note: the tolerances are not so tight because the values in the
    table appear to be slightly inconsistent with the equations in the
    text. For example, if you calculate the diffusivity from the mean
    free path and thermal velocity, you get a value that is slightly
    different from the one in the table. The tolerances are set to allow
    for these small discrepancies.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        afr = initialize_aerosol(
            aerosol_distribution="monodisperse",
            aerosol_diameter=diameter * 1e3,
            T=25,
            P=101325.0,
            P_units="Pa",
            aerosol_density=1.0,
            carrier_FR=0,
            reactant_FR=0.1,
            reactant_carrier_FR=0,
        )

    assert np.isclose(afr.aerosol_diffusion_rate, diffusivity, rtol=0.2)
    assert np.isclose(afr.thermal_particle_velocity, thermal_velocity, rtol=0.05)
    assert np.isclose(afr.aerosol_mean_free_path, mean_free_path * 1e-4, rtol=0.2)


@pytest.mark.parametrize(
    "surface_area, k",
    [
        (4e-4, 0.06),
        (3e-4, 0.04),
        (1e-4, 0.015),
    ],
)
def test_uptake_matches_Wagner_2008(
    make_constructor_kwargs, make_init_kwargs, surface_area, k
):
    """
    Test that the reaction rate matches the values in Figure 5 of
    Wagner et al., 2008 based on the reported flow conditions and uptake
    coefficient.
    """
    construtor_kwargs = make_constructor_kwargs(
        AerosolFlowReactor,
        FT_ID=4.0,
        FT_length=120.0,
        injector_ID=0.6,
        injector_OD=0.6,
        reactant_gas="N2O5",
        carrier_gas="N2",
    )
    afr = AerosolFlowReactor(**construtor_kwargs)

    init_kwargs = make_init_kwargs(
        AerosolFlowReactor,
        T=23.0,
        P=101325.0,
        P_units="Pa",
        carrier_FR=2500.0,
        reactant_FR=0.1,
        reactant_carrier_FR=500.0,
        aerosol_distribution="lognormal",
        aerosol_diameter=850.0,
        aerosol_number_conc=1e4,  # not provided, so arbitrarily chosen
        aerosol_sigma=1.5,  # not provided, so arbitrarily chosen
        aerosol_density=2.7,
        aerosol_surface_area=surface_area,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        afr.initialize(**init_kwargs)

    afr.reactant_uptake(hypothetical_gamma=2.3e-2, disp=False)
    assert np.isclose(k, afr.k_rxn, rtol=0.15)

@pytest.mark.parametrize(
    "exposure_distance, number_conc, relative_signal",
    [
        (45, 3e4, 0.5),
        (45, 1.4e4, 0.7),
        (70, 2e4, 0.4),
        (70, 0.8e4, 0.7),
        (105, 1e4, 0.5),
        (105, 2.4e4, 0.2)
    ],
)
def test_uptake_matches_Wagner_2008_with_exposure_length(
    make_constructor_kwargs, make_init_kwargs, exposure_distance, number_conc, relative_signal
):
    """
    Test that the N2O5 loss matches the values in Figure 6 of Wagner et 
    al., 2008.
    """
    construtor_kwargs = make_constructor_kwargs(
        AerosolFlowReactor,
        FT_ID=4.0,
        FT_length=120.0,
        injector_ID=0.6,
        injector_OD=0.6,
        reactant_gas="N2O5",
        carrier_gas="N2",
    )
    afr = AerosolFlowReactor(**construtor_kwargs)

    init_kwargs = make_init_kwargs(
        AerosolFlowReactor,
        T=23.0,
        P=101325.0,
        P_units="Pa",
        carrier_FR=2500.0,
        reactant_FR=0.1,
        reactant_carrier_FR=500.0,
        aerosol_distribution="lognormal",
        aerosol_diameter=850.0,
        aerosol_number_conc=number_conc,
        aerosol_sigma=1.5,  # not provided, so arbitrarily chosen
        aerosol_density=2.7,
        axial_distance=exposure_distance,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        afr.initialize(**init_kwargs)

    # Note the gamma value is taken from Table 1. It's not clear if this is the same as 
    # the gamma calculated from Figure 6, but it seems to be the case.
    afr.reactant_uptake(hypothetical_gamma=1.3e-2, exposure_length=exposure_distance, gamma_wall=1e-7, disp=False)
    assert np.isclose(afr.aerosol_loss, 1-relative_signal, rtol=0.2)