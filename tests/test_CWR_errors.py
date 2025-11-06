from flowtube import CoatedWallReactor, tools, flow_calc
import numpy as np


def test_Knopf_et_al_2015_parameters_do_not_raise(
    make_constructor_kwargs, make_init_kwargs, build_reactor
):
    """
    Test that using parameters from Knopf et al. (2015) does not raise
    errors.

    Reference:
    Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and
    Penetration of Gas Molecules and Aerosol Particles through Laminar
    Flow Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87,
    3746–3754. https://doi.org/10.1021/ac5042395
    """
    ### Table 1: O3 + BSA ###
    kwargs = make_constructor_kwargs(
        CoatedWallReactor,
        FT_ID=0.8,
        FT_length=50,
        injector_ID=0.2,
        injector_OD=0.5,
        reactant_gas="O3",
        carrier_gas="He",
        reactant_MR=1,
    )
    init_kwargs = make_init_kwargs(
        CoatedWallReactor,
        reactant_FR=0.01,
        reactant_carrier_FR=0,
        carrier_FR=1000,
        P=1013,
        P_units="hPa",
        T=295 - 273.15,
        reactant_diffusion_rate=0.1267,
        axial_delta_T=0,
        radial_delta_T=0,
        disp=False,
    )
    obj, _, _ = build_reactor(
        CoatedWallReactor,
        init_overrides=init_kwargs,
        constructor_overrides=kwargs,
    )

    flow_velocity = 30.17  # cm s-1
    flow_rate_ccm = tools.cross_sectional_area(obj.FT_ID) * flow_velocity * 60
    flow_rate_sccm = flow_calc.ccm_to_sccm(obj, flow_rate_ccm)

    init_kwargs["carrier_FR"] = flow_rate_sccm
    obj, _, _ = build_reactor(
        CoatedWallReactor,
        init_overrides=init_kwargs,
        constructor_overrides=kwargs,
    )

    assert np.isclose(obj.FT_residence_time, 1.657, rtol=0.01)
    assert np.isclose(obj.Re_FT, 17.4, rtol=0.2)
    assert np.isclose(obj.Pe_FT, 190.5, rtol=0.01)
    assert np.isclose(obj.reactant_mean_free_path, 1.1e-5, rtol=0.05)
    assert np.isclose(obj.Kn_FT, 2.8e-5, rtol=0.1)
    assert np.isclose(obj.reactant_molec_velocity, 3.6e4, rtol=0.1)

    ### Table 1: OH + Levoglucosan ###
    kwargs = make_constructor_kwargs(
        CoatedWallReactor,
        FT_ID=1.77,
        FT_length=3.6,
        injector_ID=0.2,
        injector_OD=0.5,
        reactant_gas="OH",
        carrier_gas="He",
        reactant_MR=1,
    )
    init_kwargs = make_init_kwargs(
        CoatedWallReactor,
        reactant_FR=0.01,
        reactant_carrier_FR=0,
        carrier_FR=1000,
        P=5.3,
        P_units="hPa",
        T=293 - 273.15,
        reactant_diffusion_rate=188.22,
        axial_delta_T=0,
        radial_delta_T=0,
        disp=False,
    )
    obj, _, _ = build_reactor(
        CoatedWallReactor,
        init_overrides=init_kwargs,
        constructor_overrides=kwargs,
    )

    flow_velocity = 2874.71  # cm s-1
    flow_rate_ccm = tools.cross_sectional_area(obj.FT_ID) * flow_velocity * 60
    flow_rate_sccm = flow_calc.ccm_to_sccm(obj, flow_rate_ccm)

    init_kwargs["carrier_FR"] = flow_rate_sccm
    obj, _, _ = build_reactor(
        CoatedWallReactor,
        init_overrides=init_kwargs,
        constructor_overrides=kwargs,
    )

    assert np.isclose(obj.FT_residence_time, 0.00126, rtol=0.01)
    assert np.isclose(obj.Re_FT, 22.5, rtol=0.1)
    assert np.isclose(obj.Pe_FT, 28.71, rtol=0.1)
    assert np.isclose(obj.reactant_mean_free_path, 9.35e-3, rtol=0.05)
    assert np.isclose(obj.Kn_FT, 0.01, rtol=0.1)
    assert np.isclose(obj.reactant_molec_velocity, 6.04e4, rtol=0.1)
