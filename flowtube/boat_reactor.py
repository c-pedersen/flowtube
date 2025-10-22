import numpy as np
import molmass as mm
from numpy.typing import NDArray
import warnings

from . import tools, diffusion_coef, viscosity_density, flow_calc


class BoatReactor:
    def __init__(
        self,
        FT_ID: float,
        FT_length: float,
        injector_ID: float,
        injector_OD: float,
        reactant_gas: str,
        carrier_gas: str,
        reactant_MR: float,
        boat_width: float,
        boat_height: float,
        boat_length: float,
        boat_wall_thickness: float,
    ) -> None:
        """
        Handles calculations relevant to flow rate, flow diagnostics,
        transport, and uptake for a boat reactor (partial cylinder (less
        than half) inside of a cylinder). Since there is no diffusion
        correction for the boat reactor geometry, the assumption is that
        gas-phase diffusion is negligible. Assumes that the boat reactor
        is filled to the brim with a liquid. The injector is assumed to
        be above the boat reactor at all times. All calculations are for
        over the boat.

        Args:
            FT_ID (float): Inner diameter (cm) of flow tube.
            FT_length (float): Length (cm) of flow tube.
            injector_ID (float): Inner diameter (cm) of reactant
                injector.
            injector_OD (float): Outer diameter (cm) of reactant
                injector.
            reactant_gas (str): Molecular formula of reactant gas
                (supported Ar, He, Air, Br2, Cl2, HBr, HCl, HI, H2O, I2,
                NO, N2, and O2).
            carrier_gas (str): Molecular formula of carrier gas
                (supported: Ar, He, N2, O2).
            reactant_MR (float): Reactant mixing ratio (mol mol-1).
            boat_width (float): Width (cm) of boat reactor.
            boat_height (float): Height (cm) of boat reactor.
            boat_length (float): Length (cm) of boat reactor.
            boat_wall_thickness (float): Wall thickness (cm).

        Returns:
            None
        """

        ### Check for valid inputs ###
        # Check if the gases are supported
        if reactant_gas not in diffusion_coef.sigmas.keys():
            raise ValueError(
                f"Unsupported reactant gas. "
                f"Supported gases: {', '.join(diffusion_coef.sigmas.keys())}"
            )
        if carrier_gas not in viscosity_density.a.keys():
            raise ValueError(
                f"Unsupported carrier gas. "
                f"Supported gases: {', '.join(viscosity_density.a.keys())}"
            )

        # Check physicality of insert dimensions
        if (
            boat_length < 0
            or boat_height < 0
            or boat_width < 0
            or boat_wall_thickness < 0
        ):
            raise ValueError("Boat dimensions must be positive")
        elif boat_height > FT_ID or boat_width > FT_ID:
            raise ValueError("Boat width or height cannot be larger than flow tube ID")
        elif boat_length > FT_length:
            raise ValueError("Boat length cannot be larger than flow tube length")

        # Check physicality of injector dimensions
        if injector_ID < 0 or injector_OD < 0:
            raise ValueError("Injector ID and OD must be positive")
        elif injector_ID > FT_ID:
            raise ValueError("Injector ID cannot be larger than flow tube ID")
        elif injector_OD > FT_ID:
            raise ValueError("Injector OD cannot be larger than flow tube ID")
        elif injector_ID > injector_OD:
            raise ValueError("Injector ID cannot be larger than injector OD")

        # Check mixing ratio
        if reactant_MR < 0 or reactant_MR > 1:
            raise ValueError("Reactant mixing ratio must be between 0 and 1")

        # Initialize variables
        self.FT_ID = FT_ID
        self.FT_length = FT_length
        self.injector_ID = injector_ID
        self.injector_OD = injector_OD
        self.reactant_gas = reactant_gas
        self.carrier_gas = carrier_gas
        self.reactant_MR = reactant_MR
        self.boat_height = boat_height
        self.boat_width = boat_width
        self.boat_length = boat_length
        self.boat_wall_thickness = boat_wall_thickness

    def initialize(
        self,
        reactant_FR: float,
        reactant_carrier_FR: float,
        carrier_FR: float,
        P: float,
        P_units: str,
        T: float,
        radial_delta_T: float = 1,
        axial_delta_T: float = 1,
        disp: bool = True,
    ) -> None:
        """
        Sets experimental conditions and calls calculation functions for
        numerous flow and diffusion parameters.

        Args:
            reactant_FR (float): Reactant flow rate (sccm).
            reactant_carrier_FR (float): Carrier flow rate (sccm) used
                to dilute the reactant.
            carrier_FR (float): Carrier flow rate (sccm) typically
                injected near the start of the flow tube.
            P (float): Pressure.
            P_units (str): Pressure units.
            T (float): Temperature (C).
            radial_delta_T (float): Radial temperature gradient (K)
                (default = 1 K).
            axial_delta_T (float): Axial temperature gradient (K)
                (default = 1 K).
            disp (bool): Display calculated calculated values.

        Returns:
            None
        """
        ### Check for valid inputs ###
        # Check if flow rates are positive
        if reactant_FR < 0 or reactant_carrier_FR < 0 or carrier_FR < 0:
            raise ValueError("Flow rates must be positive")

        # Check if the pressure units are supported
        if P_units not in tools.P_CF.keys():
            raise ValueError(
                f"Unsupported pressure units. "
                f"Supported units: {', '.join(tools.P_CF.keys())}"
            )
        elif P < 0:
            raise ValueError("Pressure must be positive")

        # Check if the temperature & temperature gradients are valid numbers
        if T < -273.15:
            raise ValueError("Temperature must be above absolute zero (-273.15 C)")
        if radial_delta_T < 0 or axial_delta_T < 0:
            raise ValueError("Temperature gradients must be positive")

        self.P = tools.P_in_Pa(P, P_units)
        self.T = tools.T_in_K(T)

        self.flows(reactant_FR, reactant_carrier_FR, carrier_FR, disp=disp)
        self.carrier_flow(
            radial_delta_T=radial_delta_T, axial_delta_T=axial_delta_T, disp=disp
        )
        self.reactant_diffusion(disp=disp)

    def flows(
        self,
        reactant_FR: float,
        reactant_carrier_FR: float,
        carrier_FR: float,
        disp: bool = True,
    ) -> None:
        """Calculates Flow Tube flows.

        Args:
            reactant_FR (float): Reactant flow rate (sccm).
            reactant_carrier_FR (float): Carrier flow rate (sccm) used
                to dilute the reactant.
            carrier_FR (float): Carrier flow rate (sccm) typically
                injected near the start of the flow tube.
            disp (bool): Display calculated calculated values.

        Returns:
            None
        """
        # Check if the flow rates are positive
        if reactant_FR < 0 or reactant_carrier_FR < 0 or carrier_FR < 0:
            raise ValueError("Flow rates must be positive")

        # Lists for displaying values
        var_names: list[str] = []
        var: list[float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        # Flow Rate Setpoints
        var_names += ["Reactant Flow Rate"]
        var += [reactant_FR]
        var_fmts += [".2f"]
        units += ["sccm"]
        var_names += ["Reactant Carrier Flow Rate"]
        var += [reactant_carrier_FR]
        var_fmts += [".1f"]
        units += ["sccm"]

        # Total Flow Rates
        total_reactant_FR = reactant_FR + reactant_carrier_FR
        self.total_FR = reactant_FR + reactant_carrier_FR + carrier_FR
        var_names += ["Total Reactant Flow Rate"]
        var += [total_reactant_FR]
        var_fmts += [".1f"]
        units += ["sccm"]

        # Total Reactant Flow Velocity
        total_reactant_flow_velocity = flow_calc.sccm_to_velocity(
            self, total_reactant_FR, self.injector_ID
        )

        # Calculate the cross-sectional area of the boat reactor and flow tube
        boat_cross_section = tools.partial_cylinder_area(
            self.boat_height, self.boat_width
        )
        self.net_cross_section = (
            tools.cross_sectional_area(self.FT_ID) - boat_cross_section
        )

        # Minimum Carrier Flow Velocity & Rate
        # - to prevent effect mentioned in Li et al., ACP, 2020
        min_carrier_flow_velocity = total_reactant_flow_velocity * 1.33
        min_carrier_FR = flow_calc.ccm_to_sccm(
            self,
            min_carrier_flow_velocity
            * (self.net_cross_section - tools.cross_sectional_area(self.injector_OD))
            * 60,
        )
        var_names += ["Minimum Carrier Flow Rate"]
        var += [min_carrier_FR]
        var_fmts += [".1f"]
        units += ["sccm"]
        if carrier_FR < min_carrier_FR:
            warnings.warn(
                "Carrier flow rate is below the minimum. "
                "This may affect the flow profile in the flow tube."
            )

        # More Flow Rates
        var_names += ["Carrier Flow Rate"]
        var += [carrier_FR]
        var_fmts += [".1f"]
        units += ["sccm"]
        var_names += ["Total Flow Rate"]
        var += [self.total_FR]
        var_fmts += [".1f"]
        units += ["sccm"]

        # Reactant Concentrations (ppb)
        injector_conc = reactant_FR / total_reactant_FR * self.reactant_MR * 1e9
        FT_conc = reactant_FR / self.total_FR * self.reactant_MR * 1e9
        FT_conc_molec = flow_calc.MR_to_molec(self, FT_conc)
        var_names += [f"Injector {self.reactant_gas} Concentration"]
        var += [injector_conc]
        var_fmts += [".3g"]
        units += ["ppb"]
        var_names += [f"Flow Tube {self.reactant_gas} Concentration"]
        var += [FT_conc]
        var_fmts += [".3g"]
        units += ["ppb"]
        var_names += [f"Flow Tube {self.reactant_gas} Concentration"]
        var += [FT_conc_molec]
        var_fmts += [".2e"]
        units += ["molec. cm-3"]

        # Total Flow Velocities
        self.flow_velocity = (
            flow_calc.sccm_to_ccm(self, self.total_FR) / self.net_cross_section / 60
        )
        var_names += ["Flow Tube Velocity"]
        var += [self.flow_velocity]
        var_fmts += [".3g"]
        units += ["cm s-1"]

        # Residence Time over the boat
        self.residence_time = (
            self.boat_length - self.boat_wall_thickness * 2
        ) / self.flow_velocity
        var_names += ["Boat Residence Time"]
        var += [self.residence_time]
        var_fmts += [".3g"]
        units += ["s"]

        ### Display Values ###
        if disp:
            tools.table(
                "Flow Setpoints and Conditions",
                var_names,
                var,
                var_fmts,
                units,
            )

    def carrier_flow(
        self,
        radial_delta_T: float = 1,
        axial_delta_T: float = 1,
        disp: bool = True,
    ):
        """Performs and displays carrier gas transport calculations.

        Args:
            delta_T_radial (float): Radial temperature gradient (K).
            delta_T_axial (float): Axial temperature gradient (K).
            disp (bool): Display calculated values.

        Returns:
            None
        """

        # Lists for displaying values
        var_names: list[str] = []
        var: list[float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        # Carrier Gas Dynamic Viscosity (kg m-1 s-1)
        self.carrier_dynamic_viscosity = viscosity_density.dynamic_viscosity(
            self, self.carrier_gas
        )
        var_names += ["Carrier Gas Dynamic Viscosity"]
        var += [self.carrier_dynamic_viscosity]
        var_fmts += [".2e"]
        units += ["kg m-1 s-1"]

        # Carrier Gas Density (kg m-3)
        self.carrier_density = viscosity_density.real_density(self, self.carrier_gas)
        var_names += ["Carrier Gas Density"]
        var += [self.carrier_density]
        var_fmts += [".3g"]
        units += ["kg m-3"]

        # Reynolds Number - laminar flow if Re < 1800
        Re_FT = flow_calc.reynolds_number(self, self.total_FR, self.FT_ID)
        var_names += ["!!!!!!!!!Flow Tube Reynolds Number"]
        var += [Re_FT]
        var_fmts += [".0f"]
        units += ["unitless"]
        if Re_FT > 1800:
            warnings.warn("Re > 1800. Flow in flow tube may not be laminar")

        # Entrance length (cm) - see flow_calc.py for details
        length_to_laminar = flow_calc.length_to_laminar(self.FT_ID, Re_FT)
        var_names += ["Flow Tube Entrance length"]
        var += [length_to_laminar]
        var_fmts += [".1f"]
        units += ["cm"]

        # Pressure Gradient (%) - see flow_calc.py for details
        FT_pressure_gradient = flow_calc.pressure_gradient(
            self, flow_calc.conductance(self, self.FT_ID, self.FT_length), self.total_FR
        )
        var_names += ["!!!!!!!!Flow Tube Pressure Gradient"]
        var += [FT_pressure_gradient * 100]
        var_fmts += [".2f"]
        units += ["%"]

        # Buoyancy Parameters - see flow_calc.py for details
        radial_buoyancy = flow_calc.buoyancy_parameters(
            self, radial_delta_T, self.FT_ID, Re_FT
        )
        axial_buoyancy = flow_calc.buoyancy_parameters(
            self, axial_delta_T, self.FT_length, Re_FT
        )
        var_names += [f"Radial Buoyancy Parameter (ΔT={radial_delta_T:.1f} C)"]
        var += [radial_buoyancy]
        var_fmts += [".2f"]
        units += ["unitless"]
        var_names += [f"Axial Buoyancy Parameter (ΔT={axial_delta_T:.1f} C)"]
        var += [axial_buoyancy]
        var_fmts += [".2f"]
        units += ["unitless"]
        if radial_buoyancy > 1:
            warnings.warn(
                "Radial buoyancy parameter > 1. "
                "Flow may be affected by buoyancy effects"
            )
        if axial_buoyancy > 1:
            warnings.warn(
                "Axial buoyancy parameter > 1. Flow may be affected by buoyancy effects"
            )

        ### Display Values ###
        if disp:
            tools.table(
                "Fluid Dynamics of Carrier Gas",
                var_names,
                var,
                var_fmts,
                units,
            )

    def reactant_diffusion(self, disp: bool = True):
        """Performs and displays reactant diffusion calculations.

        Args:
            disp (bool): Display calculated calculated values.

        Returns:
            None
        """

        # Lists for displaying values
        var_names: list[str] = []
        var: list[float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        # Reactant Diffusion Rate (cm2 s-1)
        self.reactant_diffusion_rate = diffusion_coef.binary_diffusion_coefficient(self)
        var_names += ["Reactant Diffusion Rate"]
        var += [self.reactant_diffusion_rate]
        var_fmts += [".3g"]
        units += ["cm2 s-1"]

        # Thermal Molecular Velocity (cm s-1)
        # - formula matched to values from Knopf et al., Anal. Chem., 2015
        self.reactant_molec_velocity = flow_calc.molec_velocity(
            self, float(mm.Formula(self.reactant_gas).mass)
        )  # pyright: ignore[reportUnknownArgumentType, reportUnknownMemberType]

        # Reactant Mean Free Path (cm) - Fuchs and Sutugin, 1971
        reactant_mean_free_path = (
            3 * self.reactant_diffusion_rate / self.reactant_molec_velocity
        )

        # Advection Rate (cm2 s-1) - eq. 1 from Knopf et al., Anal. Chem., 2015
        advection_rate = self.flow_velocity * self.FT_ID
        var_names += ["Flow Tube Advection Rate"]
        var += [advection_rate]
        var_fmts += [".3g"]
        units += ["cm2 s-1"]

        # Peclet Number - if > 10 then axial diffusion is negligible
        # - eq. 1 from Knopf et al., Anal. Chem., 2015
        Pe = advection_rate / self.reactant_diffusion_rate
        var_names += ["Peclet Number"]
        var += [Pe]
        var_fmts += [".4g"]
        units += ["unitless"]
        if Pe < 10:
            warnings.warn("Pe < 10. Axial diffusion is non-negligible")

        # Mixing Time (s) - see flow_calc.py for details
        mixing_time = flow_calc.mixing_time(self, self.FT_ID)
        var_names += ["Flow Tube Mixing Time (upper limit)"]
        var += [mixing_time]
        var_fmts += [".2g"]
        units += ["s"]

        # Mixing Length (cm)
        mixing_length = self.flow_velocity * mixing_time
        var_names += ["Flow Tube Mixing Length"]
        var += [mixing_length]
        var_fmts += [".2g"]
        units += ["cm"]

        # Effective Sherwood Number (unitless)
        # - eq. 11 from Knopf et al., Anal. Chem., 2015
        # Note: the boat geometry is not considered and thus this value
        # should be used as a limiting case
        self.N_eff_Shw_FT = flow_calc.N_eff_Shw(self, self.FT_length, self.total_FR)

        # Knudsen Number for reactant-wall/insert interaction
        # - eq. 8 from Knopf et al., Anal. Chem., 2015
        # Note: the boat geometry is not considered and thus this value
        # should be used as a limiting case
        self.Kn_FT = flow_calc.Kn(reactant_mean_free_path, self.FT_ID)

        ### Display Values ###
        if disp:
            tools.table(
                "Reactant Diffusion Parameters",
                var_names,
                var,
                var_fmts,
                units,
            )

    def reactant_uptake(
        self,
        hypothetical_gamma: NDArray[np.float64] | float,
        gamma_wall: float = 5e-6,
        disp: bool = True,
    ) -> NDArray[np.float64] | float:
        """
        Calculates reactant uptake to the boat and loss to flow tube
        walls.

        Args:
            hypothetical_gamma (float or numpy.ndarray): Hypothetical
                uptake coefficient to calculate diffusion correction
                factor.
            gamma_wall (float): Wall uptake coefficient (default: 5e-6
                for halocarbon wax coating - Ivanov et al., 2021).
            disp (bool): Display calculated values.

        Returns:
            float: Loss to boat as ratio to inital amount (fraction).
        """

        ### Check for valid inputs ###
        # Check if hypothetical_gamma is between 0 and 1
        if np.min(hypothetical_gamma) < 0 or np.max(hypothetical_gamma) > 1:  # pyright: ignore[reportUnknownMemberType]
            raise ValueError("Hypothetical gamma must be between 0 and 1")

        # Check if gamma_wall is between 0 and 1
        if gamma_wall < 0 or gamma_wall > 1:
            raise ValueError("Wall uptake coefficient must be between 0 and 1")

        # Lists for displaying values
        var_names: list[str] = []
        var: list[NDArray[np.float64] | float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        # Boat surface area and volume
        liquid_width = self.boat_width - self.boat_wall_thickness * 2
        liquid_length = self.boat_length - self.boat_wall_thickness * 2
        liquid_height = self.boat_height - self.boat_wall_thickness
        liquid_surface_area = liquid_width * liquid_length
        liquid_volume = (
            tools.partial_cylinder_area(liquid_height, liquid_width) * liquid_length
        )
        var_names += ["Boat Surface Area", "Boat Volume"]
        var += [liquid_surface_area, liquid_volume]
        var_fmts += [".1f", ".1f"]
        units += ["cm2", "cm3"]

        # Diffusion Correction - see flow_calc.py for details
        diff_corr = 1 - flow_calc.correction_factor(
            self.N_eff_Shw_FT, self.Kn_FT, hypothetical_gamma
        )
        if diff_corr > 0.05:
            warnings.warn(
                "Diffusion correction is > 5%. "
                "Negligible diffusion may no longer be a valid assumption"
            )
        var_names += [
            "Flow Tube Wall Diffusion Correction "
            "\n(must be small to neglect for boat reactor)"
        ]
        var += [diff_corr * 100]
        var_fmts += [".1f"]
        units += ["%"]

        # Geometric correction for boat geometry – Hanson and Ravishankara, 1993
        cylinder_SA_V_ratio = 4 / self.FT_ID
        actual_SA_V_ratio = liquid_surface_area / (
            self.net_cross_section * liquid_length
        )
        geometric_correction = actual_SA_V_ratio / cylinder_SA_V_ratio
        var_names += ["Boat geometry correction factor"]
        var += [1 / geometric_correction]
        var_fmts += [".2f"]
        units += ["unitless"]

        # Corrected Loss Rate (s-1)
        # - see flow_calc.py and Hanson and Ravishankara, 1993 for details
        k = (
            flow_calc.observed_loss_rate(self, self.FT_ID, hypothetical_gamma)
            / geometric_correction
        )
        var_names += ["Loss Rate"]
        var += [k]
        var_fmts += [".3g"]
        units += ["s-1"]

        # Uptake to boat (fraction) - first order kinetics
        uptake = 1 - np.exp(-k * self.residence_time)
        var_names += ["Loss to Boat"]
        var += [uptake * 100]
        var_fmts += [".1f"]
        units += ["%"]

        # Reactant Wall Loss (fraction)
        # - calculated as if there is no boat, take as an upper limit
        reactant_wall_loss = flow_calc.cylinder_loss(
            self,
            self.FT_ID,
            self.N_eff_Shw_FT,
            self.Kn_FT,
            gamma_wall,
            self.residence_time,
        )
        var_names += ["Estimated Wall Loss (upper limit)"]
        var += [reactant_wall_loss * 100]
        var_fmts += [".2g"]
        units += ["%"]

        ### Display Values ###
        if disp is True or not isinstance(var, np.ndarray):
            tools.table(
                "Reactant Uptake",
                var_names,
                var,  # pyright: ignore[reportArgumentType]
                var_fmts,
                units,
            ) 

        return uptake


""" 
Citations:
    
Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and 
Penetration of Gas Molecules and Aerosol Particles through Laminar Flow 
Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87, 3746–3754. 
https://doi.org/10.1021/ac5042395

Hanson, D., Kosciuch, E., 2003. The NH3 Mass Accommodation Coefficient 
for Uptake onto Sulfuric Acid Solutions. J. Phys. Chem. A 107, 
2199–2208. https://doi.org/10.1021/jp021570j

Hanson, D.R., Ravishankara, A.R., 1993. Uptake of hydrochloric acid and 
hypochlorous acid onto sulfuric acid: solubilities, diffusivities, and 
reaction. J. Phys. Chem. 97, 12309–12319. 
https://doi.org/10.1021/j100149a035

Fuchs, N.A., Sutugin, A.G., 1971. HIGH-DISPERSED AEROSOLS, in: Hidy, 
G.M., Brock, J.R. (Eds.), Topics in Current Aerosol Research, 
International Reviews in Aerosol Physics and Chemistry. Pergamon, p. 1. 
https://doi.org/10.1016/B978-0-08-016674-2.50006-6

Ivanov, A.V., Molina, M.J., Park, J., 2021. Experimental study on HCl 
uptake by MgCl2 and sea salt under humid conditions. J Mass Spectrom 56,
 e4601. https://doi.org/10.1002/jms.4601

Tang, M.J., Cox, R.A., Kalberer, M., 2014. Compilation and evaluation of
gas phase diffusion coefficients of reactive trace gases in the 
atmosphere: volume 1. Inorganic compounds. Atmos. Chem. Phys. 14, 
9233–9247. https://doi.org/10.5194/acp-14-9233-2014
"""
