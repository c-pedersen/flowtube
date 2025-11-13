"""
Conversion functions for flow rates and concentrations in flowtube
package.
"""

from typing import Protocol
from numpy.typing import NDArray
import numpy as np
from . import tools


class basic_attrs(
    Protocol,
):
    P: float  # Pressure in Pa
    T: float  # Temperature in K


class carrier_attrs(
    basic_attrs,
    Protocol,
):
    carrier_dynamic_viscosity: (
        float  # Dynamic viscosity of the carrier gas (kg m-1 s-1)
    )
    carrier_density: float  # Density of the carrier gas (kg m-3)


class full_attrs(
    carrier_attrs,
    Protocol,
):
    reactant_diffusion_rate: float  # Diffusion rate of the reactant (cm2 s-1)
    reactant_molec_velocity: (
        float  # Thermal molecular velocity of the reactant (cm s-1)
    )


### Flow rate conversion functions ###
def sccm_to_ccm(
    obj: basic_attrs,
    FR: float,
) -> float:
    """Convert sccm to cm3 per min.

    Args:
        obj (basic_attrs): Object with basic attributes
            (P in Pa, T in K).
        FR (float): Flow rate in sccm.

    Returns:
        float: Flow rate in cm3 min-1
    """

    return (
        (tools.STANDARD_PRESSURE_Pa / obj.P) * obj.T / tools.STANDARD_TEMPERATURE_K * FR
    )


def ccm_to_sccm(
    obj: basic_attrs,
    FR: float,
) -> float:
    """Convert cm3 per min to sccm

    Args:
        obj (basic_attrs): Object with basic attributes
            (P in Pa, T in K).
        FR (float): Flow rate in cm3 min-1.

    Returns:
        float: Flow rate to sccm.
    """

    return (
        (obj.P / tools.STANDARD_PRESSURE_Pa) * tools.STANDARD_TEMPERATURE_K / obj.T * FR
    )


def sccm_to_velocity(
    obj: basic_attrs,
    FR: float,
    diameter: float,
) -> float:
    """Calculate flow velocity.

    Args:
        obj (basic_attrs): Object with basic attributes
            (P in Pa, T in K).
        FR (float): Flow rate in sccm.
        diameter (float): Diameter in cm.

    Returns:
        float: Flow velocity in cm s-1.
    """

    return sccm_to_ccm(obj, FR) / tools.cross_sectional_area(diameter) / 60


### Physical calculations ###
def MR_to_molec(
    obj: basic_attrs,
    conc: float,
) -> float:
    """Convert mixing ratio in ppb to molecules cm-3.

    Args:
        obj (basic_attrs): Object with basic attributes
            (P in Pa, T in K).
        conc (float): Mixing ratio (ppb - mol mol-1).

    Returns:
        float: Concentration in molec. cm-3.
    """

    return (
        obj.P
        / (tools.UNIVERSAL_GAS_CONSTANT * obj.T)
        * tools.AVOGADROS_NUMBER
        / 100**3
        * conc
        * 1e-9
    )


def molec_velocity(
    obj: basic_attrs,
    molar_mass: float,
) -> float:
    """
    Calculate thermal molecular velocity. Formula matched to values from
    Knopf et al., Anal. Chem., 2015

    Args:
        obj (basic_attrs): Object with basic attributes
            (P in Pa, T in K).
        molar_mass (float): Molar mass of the gas (g mol-1).

    Returns:
        float: Thermal molecular velocity in cm s-1.
    """

    return 100 * np.sqrt(
        8 / np.pi * tools.UNIVERSAL_GAS_CONSTANT * obj.T / molar_mass * 1000
    )


### Flow diagnostics calculations ###
def reynolds_number(
    obj: carrier_attrs,
    FR: float,
    diameter: float,
) -> float:
    """Calculate Reynolds number for a gas flowing through a cylinder.

    Args:
        obj (carrier_attrs): Object with full attributes (P in Pa,
            T in K, carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        FR (float): Total flow rate in sccm.
        diameter (float): Diameter of the cylinder (cm).

    Returns:
        float: Reynolds number.
    """

    return (
        (obj.carrier_density / 100**3)
        * sccm_to_velocity(obj, FR, diameter)
        * diameter
        / (obj.carrier_dynamic_viscosity / 100)
    )


def reynolds_number_irregular(
    obj: carrier_attrs,
    cross_sectional_area: float,
    wetted_perimeter: float,
    FR: float,
) -> float:
    """
    Calculate Reynolds number for a gas flowing through a cylinder.
    Formula 6-14 from Holman and Bhattacharyya 2011.

    Args:
        obj (carrier_attrs): Object with full attributes (P in Pa,
            T in K, carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        cross_sectional_area (float): Cross sectional area of the flow
            passage (cm2).
        wetted_perimeter (float): Wetted perimeter of the flow passage
            (cm).
        FR (float): Total flow rate in sccm.

    Returns:
        float: Reynolds number.
    """
    hydraulic_diameter = 4 * cross_sectional_area / wetted_perimeter

    return (
        (obj.carrier_density / 100**3)
        * sccm_to_ccm(obj, FR)
        / cross_sectional_area
        / 60
        * hydraulic_diameter
        / (obj.carrier_dynamic_viscosity / 100)
    )


def conductance(
    obj: carrier_attrs,
    diameter: float,
    length: float,
) -> float:
    """
    Calculate conductance through cylinder - eq. 3.17 from Moore et al.,
    2009.

    Args:
        obj (carrier_attrs): Object with full attributes (P in Pa,
            T in K, carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        diameter (float): Inner diameter (cm).
        length (float): Length (cm).

    Returns:
        float: Conductance in L s-1.
    """

    return (
        32600
        * diameter**4
        / (obj.carrier_dynamic_viscosity * 1e7 * length)
        * obj.P
        / tools.P_CF["Torr"]
    )


def pressure_gradient(
    obj: basic_attrs,
    conductance: float,
    FR: float,
) -> float:
    """
    Calculate pressure gradient - eqs. 3.9 & 3.10 from Moore et al.,
    2009.

    Args:
        obj (basic_attrs): Object with basic attributes
            (P in Pa, T in K).
        conductance (float): Conductance in L s-1.
        flow_rate (float): Flow rate in sccm.

    Returns:
        float: Pressure gradient ratio.
    """

    return sccm_to_ccm(obj, FR) / 60 / 1000 / conductance


def buoyancy_parameters(
    obj: carrier_attrs,
    delta_T: float,
    distance: float,
    Re: float,
) -> float:
    """Calculate buoyancy parameters.

    Args:
        obj (carrier_attrs): Object with full attributes (P in Pa,
            T in K, carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        delta_T (float): Temperature difference (K).
        distance (float): Distance over which the temperature difference
            is measured (cm) (typically axial or radial).
        Re (float): Reynolds number of the flow tube.

    Returns:
        float: Buoyancy parameter (>1 indicates the flow being driven by
            buoyancy).
    """

    # Grashof Number - eq. 9.12 from Incropera, et al., 2007
    grashof_number = (
        9.81
        / obj.T
        * delta_T
        * distance**3
        / (obj.carrier_dynamic_viscosity / 100 / (obj.carrier_density / 100**3)) ** 2
    )

    return grashof_number / Re**2


def length_to_laminar(
    diameter: float,
    Re: float,
) -> float:
    """
    Entrance length (cm) - length to achieve laminar profile - Bird et
    al., 2002, page 52. Note that the scalar term can vary depending on
    the source. Keyser, 1994 gives a value of 0.0565 for 99% attainment
    of parabolic profile while Hanson and Kosciuch, 2003 provide a value
    of 0.05 for 95% attainment. Using an exponential fit of the values,
    it seems that the 0.035 figure that Bird et al., 2002 gives is for a
    85% attainment. Choose whichever value is most appropriate for your
    purposes.

    Args:
        diameter (float): Diameter of the cylinder (cm).
        Re (float): Reynolds number of the flow tube.

    Returns:
        float: Length to laminar profile (cm).
    """

    return 0.05 * diameter * Re


def mixing_time(
    obj: full_attrs,
    diameter: float,
) -> float:
    """
    Calculate mixing time (s) - Hanson and Lovejoy, Geophys. Res. Lett.,
    1994.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        diameter (float): Diameter of the cylinder (cm).

    Returns:
        float: Mixing time in seconds.
    """

    return (diameter / 2) ** 2 / (5 * obj.reactant_diffusion_rate)


### KPS Method Calculations ###
def N_eff_Shw(
    obj: full_attrs,
    length: float,
    FR: float,
) -> float:
    """
    Calculate the effective Sherwood number - eq. 11 from Knopf et al.,
    2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        length (float): Length of the flow tube (cm).
        FR (float): Total flow rate in cm3 min-1.

    Returns:
        float: Effective Sherwood number.
    """

    # Axial Distance (unitless)
    # - eq. 2 from Knopf et al., Anal. Chem., 2015
    z_star = (
        length * np.pi / 2 * obj.reactant_diffusion_rate / (sccm_to_ccm(obj, FR) / 60)
    )

    return 3.6568 + 0.0978 / (z_star + 0.0154)


def Kn(
    mean_free_path: float,
    char_length: float,
) -> float:
    """Calculate Knudsen number - eq. 8 from Knopf et al., 2015.

    Args:
        mean_free_path (float): Mean free path of the reactant (cm).
        char_length (float): Characteristic length: diameter of cylinder
            for a coated wall reactor (cm).

    Returns:
        float: Knudsen number.
    """

    return 2 * mean_free_path / char_length


def diffusion_limited_rate_constant(
    obj: full_attrs,
    N_eff_Shw: float,
    diameter: float,
) -> float:
    """
    Calculate diffusion limited rate constant (s-1) - eq. 10 from Knopf
    et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        N_eff_Shw (float): Effective Sherwood number.
        diameter (float): Diameter of the cylinder (cm).

    Returns:
        float: Diffusion limited rate constant (s-1).
    """

    return 4 * N_eff_Shw * obj.reactant_diffusion_rate / diameter**2


def diffusion_limited_uptake_coefficient(
    obj: full_attrs,
    diameter: float,
    k_diff: float,
) -> float:
    """
    Calculate diffusion limited effective uptake coefficient - eq. 19
    from Knopf et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        diameter (float): Diameter of the cylinder (cm).
        k_diff (float): Diffusion limited rate constant (s-1).

    Returns:
        float: Diffusion limited effective uptake coefficient (cm s-1).
    """

    return diameter / obj.reactant_molec_velocity * k_diff


def correction_factor(
    N_eff_Shw: float,
    Kn: float,
    gamma: NDArray[np.float64] | float,
) -> NDArray[np.float64] | float:
    """
    Calculate correction factor (gamma_eff / gamma) for uptake 
    coefficient - eq. 15 from Knopf et al., 2015.

    Args:
        N_eff_Shw (float): Effective Sherwood number (unitless).
        Kn (float): Knudsen number (unitless).
        hypothetical_gamma (float): Hypothetical uptake coefficient
            (unitless).

    Returns:
        float: Correction factor (unitless: gamma_eff / gamma).
    """

    return 1 / (1 + gamma * 3 / (2 * N_eff_Shw * Kn))


def observed_loss_rate(
    obj: full_attrs,
    diameter: float,
    gamma_eff: NDArray[np.float64] | float,
) -> NDArray[np.float64] | float:
    """
    Calculate observed first order loss rate (s-1) - eq. 19 from Knopf 
    et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        diameter (float): Diameter of the cylinder (cm).
        gamma_eff (float): Effective uptake coefficient (unitless).

    Returns:
        float: Observed first order loss rate (s-1).
    """

    return gamma_eff * obj.reactant_molec_velocity / diameter


def cylinder_loss(
    obj: full_attrs,
    diameter: float,
    N_eff_Shw: float,
    Kn: float,
    gamma: NDArray[np.float64] | float,
    time: float,
) -> NDArray[np.float64] | float:
    """
    Calculate amount of gas lost to cylinder walls according to a 
    provided gamma - eq. 21 from Knopf et al., 2015.

    Args:
        obj (full_attrs): Object with full attributes (P in Pa, T in K,
            reactant_diffusion_rate in cm2 s-1,
            carrier_dynamic_viscosity in kg m-1 s-1,
            carrier_density in kg m-3).
        time (float): Residence time in cylinder (s).

    Returns:
        float: Penetration - fraction of initial reactant after passing
            through cylinder (unitless).
    """

    return 1 - np.exp(
        -gamma
        / (1 + gamma * 3 / (2 * N_eff_Shw * Kn))
        * obj.reactant_molec_velocity
        / diameter
        * time
    )


""" 
Citations:
    
Bird, R.B., Stewart, W.E., Lightfoot, E.N., 2002. Transport phenomena, 
2nd, Wiley international ed ed. J. Wiley, New York.

Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and 
Penetration of Gas Molecules and Aerosol Particles through Laminar Flow 
Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87, 3746–3754. 
https://doi.org/10.1021/ac5042395

Keyser, L.F., 1984. High-pressure flow kinetics. A study of the hydroxyl
+ hydrogen chloride reaction from 2 to 100 torr. J. Phys. Chem. 88, 
4750–4758. https://doi.org/10.1021/j150664a061

Hanson, D.R., Lovejoy, E.R., 1994. The uptake of N2O5 onto small 
sulfuric acid particles. Geophys. Res. Lett. 21, 2401–2404. 
https://doi.org/10.1029/94GL02288

Hanson, D., Kosciuch, E., 2003. The NH3 Mass Accommodation Coefficient 
for Uptake onto Sulfuric Acid Solutions. J. Phys. Chem. A 107, 
2199–2208. https://doi.org/10.1021/jp021570j

Holman, J. P., & Bhattacharyya, S. (2011). Heat transfer in SI units 
(10th ed.). McGraw-Hill. p. 284

Moore, J.H., Davis, C.C., Coplan, M.A., 2009. Building Scientific 
Apparatus, 4th ed. ed. Cambridge University Press, Leiden.

Incropera, F.P., DeWitt, D.P., Bergman, T.L., Lavine, A.S. (Eds.), 2007.
Fundamentals of heat and mass transfer, 6. ed. ed. Wiley, Hoboken, NJ.

"""
