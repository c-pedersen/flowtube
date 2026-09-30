"""
Main aerosol flow reactor class and associated calculations.

Citations:
Bertram, A.K., Ivanov, A.V., Hunter, M., Molina, L.T., Molina, M.J.,
2001. The Reaction Probability of OH on Organic Surfaces of Tropospheric
Interest. J. Phys. Chem. A 105, 9415-9421.
https://doi.org/10.1021/jp0114034

Hanson, D., Kosciuch, E., 2003. The NH3 Mass Accommodation Coefficient
for Uptake onto Sulfuric Acid Solutions. J. Phys. Chem. A 107,
2199–2208. https://doi.org/10.1021/jp021570j

Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and
Penetration of Gas Molecules and Aerosol Particles through Laminar Flow
Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87, 3746-3754.
https://doi.org/10.1021/ac5042395

Hanson, D.R., Ravishankara, A.R., 1993. Uptake of hydrochloric acid and
hypochlorous acid onto sulfuric acid: solubilities, diffusivities, and
reaction. J. Phys. Chem. 97, 12309-12319.
https://doi.org/10.1021/j100149a035

Fuchs, N.A., Sutugin, A.G., 1971. HIGH-DISPERSED AEROSOLS, in: Hidy,
G.M., Brock, J.R. (Eds.), Topics in Current Aerosol Research,
International Reviews in Aerosol Physics and Chemistry. Pergamon, p. 1.
https://doi.org/10.1016/B978-0-08-016674-2.50006-6

Tang, M.J., Cox, R.A., Kalberer, M., 2014. Compilation and evaluation of
gas phase diffusion coefficients of reactive trace gases in the
atmosphere: volume 1. Inorganic compounds. Atmos. Chem. Phys. 14,
9233-9247. https://doi.org/10.5194/acp-14-9233-2014
"""

from __future__ import annotations

import warnings

import molmass as mm
import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import (
    diffusion_coef,
    flow_calc,
    input_validation,
    kinetics,
    tools,
    viscosity_density,
)

# Attributes than can be updated after initialization and will trigger a re-initialization of the reactor
_CTOR_ATTRS = frozenset(
    {
        "FT_ID",
        "FT_length",
        "injector_ID",
        "injector_OD",
        "reactant_gas",
        "carrier_gas",
        "reactant_conc_type",
        "reactant_conc",
        "reactant_FR",
        "reactant_carrier_FR",
        "carrier_FR",
        "P",
        "P_units",
        "T",
        "manually_inputted_diffusion_rate",
        "reactant_diffusion_rate",
        "radial_delta_T",
        "axial_distance",
        "aerosol_distribution",
        "aerosol_diameter",
        "aerosol_number_conc",
        "aerosol_density",
        "aerosol_sigma",
        "aerosol_surface_area",
        "manually_inputted_surface_area",
    }
)


class AerosolFlowReactor:
    def __init__(
        self,
        FT_ID: float,
        FT_length: float,
        injector_ID: float,
        injector_OD: float,
        reactant_gas: str,
        carrier_gas: str,
        reactant_conc_type: str,
        reactant_conc: float,
    ) -> None:
        """
        Handles calculations relevant to flow rate, flow diagnostics,
        transport, and uptake for a aerosol flow reactor.

        Args:
            FT_ID (float): Inner diameter (cm) of flow tube.
            FT_length (float): Length (cm) of flow tube.
            injector_ID (float): Inner diameter (cm) of reactant
                injector.
            injector_OD (float): Outer diameter (cm) of reactant
                injector.
            reactant_gas (str): Molecular formula of reactant gas
                (supported Ar, He, Air, Br2, Cl2, HBr, HCl, HI, H2O, I2,
                NO, N2, O2, ClONO2/ClNO3, N2O5, O3, and NO2).
            carrier_gas (str): Molecular formula of carrier gas
                (supported: Ar, He, N2, O2).
            reactant_conc_type (str): Type of reactant concentration
                input. Options: "ppm" or "ppb" for mixing ratio,
                "ng/min" for permeation rate, "Pa" for vapor pressure.
            reactant_conc (float): Reactant concentration value.

        Returns:
            None
        """
        # Flag to prevent calling __setattr__ before initialization is complete
        object.__setattr__(self, "_initializing", True)

        ### Initialize variables ###
        self.FT_ID = FT_ID
        self.FT_length = FT_length
        self.injector_ID = injector_ID
        self.injector_OD = injector_OD
        self.reactant_gas = reactant_gas
        self.reactant_conc_type = reactant_conc_type
        self.reactant_conc = reactant_conc
        self.carrier_gas = carrier_gas

        # Validate inputs
        input_validation.validate_init(self)

        # Turn flag off to allow __setattr__ to be used normally
        object.__setattr__(self, "_initializing", False)

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name in _CTOR_ATTRS and not self.__dict__.get("_initializing", True):
            input_validation.validate_init(self)  # validate before re-init

            missing = [a for a in _CTOR_ATTRS if not hasattr(self, a)]
            if missing:
                raise RuntimeError(
                    f"initialize() has not been called yet. Missing: {missing}"
                )

            self.initialize(
                reactant_FR=self.reactant_FR,
                reactant_carrier_FR=self.reactant_carrier_FR,
                carrier_FR=self.carrier_FR,
                P=self.P,
                P_units=self.P_units,
                T=self.T,
                axial_distance=self.axial_distance,
                aerosol_distribution=self.aerosol_distribution,
                aerosol_diameter=self.aerosol_diameter,
                aerosol_number_conc=self.aerosol_number_conc,
                aerosol_density=self.aerosol_density,
                aerosol_sigma=self.aerosol_sigma,
                aerosol_surface_area=self.aerosol_surface_area,
                radial_delta_T=self.radial_delta_T,
                disp=False,
            )

    def initialize(
        self,
        reactant_FR: float,
        reactant_carrier_FR: float,
        carrier_FR: float,
        P: float,
        P_units: str,
        T: float,
        axial_distance: float,
        aerosol_distribution: str,
        aerosol_diameter: float,
        aerosol_number_conc: float,
        aerosol_density: float,
        aerosol_sigma: float = np.nan,
        aerosol_surface_area: float = np.nan,
        reactant_diffusion_rate: float = np.nan,
        radial_delta_T: float = 1,
        disp: bool = True,
    ) -> None:
        """
        Sets experimental conditions and calls calculation functions for
        numerous flow and diffusion parameters. The aerosol is assumed
        to be mixed with the carrier gas so the inputted concentrations
        should reflect that.

        Args:
            reactant_FR (float): Reactant flow rate (sccm).
            reactant_carrier_FR (float): Carrier flow rate (sccm) used
                to dilute the reactant.
            carrier_FR (float): Carrier flow rate (sccm) typically
                injected near the start of the flow tube.
            P (float): Pressure.
            P_units (str): Pressure units.
            T (float): Temperature (C).
            axial_distance (float): Axial distance of exposed reactant
                surface (cm). Also referred to as z.
            aerosol_distribution (str): Aerosol distribution type.
                Options: lognormal or monodisperse.
            aerosol_diameter (float): Aerosol diameter (nm).
            aerosol_number_conc (float): Aerosol number concentration
                (cm-3).
            aerosol_density (float): Aerosol density (g cm-3).
            aerosol_sigma (float): Aerosol geometric standard deviation
                (unitless), required for lognormal distribution.
            aerosol_surface_area (float): Aerosol surface area (cm2
                cm-3), optional input if overriding the calculated value
                from the aerosol distribution, diameter, and number
                concentration.
            reactant_diffusion_rate (float): Reactant diffusion rate
                (cm2 s-1) (optional).
            radial_delta_T (float): Radial temperature gradient (K)
                (default = 1 K).
            disp (bool): Display calculated calculated values.

        Returns:
            None
        """
        # Flag to prevent calling __setattr__ before initialization is complete
        object.__setattr__(self, "_initializing", True)

        try:
            self.P = P
            self.P_units = P_units
            self.T = T
            self.reactant_FR = reactant_FR
            self.reactant_carrier_FR = reactant_carrier_FR
            self.carrier_FR = carrier_FR
            self.radial_delta_T = radial_delta_T
            self.axial_distance = float(axial_distance)
            self.aerosol_distribution = aerosol_distribution
            self.aerosol_diameter = aerosol_diameter
            self.aerosol_number_conc = aerosol_number_conc
            self.aerosol_density = aerosol_density
            self.aerosol_sigma = aerosol_sigma

            # Validate inputs
            input_validation.validate_initialize(self)

            ### Calculated Properties ###
            self.P_Pa = tools.P_in_Pa(self.P, self.P_units)
            self.T_K = tools.T_in_K(self.T)

            # Calculate reactant mixing ratio from input concentration
            if self.reactant_conc_type == "ppm":
                self.reactant_MR = self.reactant_conc * 1e-6
            elif self.reactant_conc_type == "ppb":
                self.reactant_MR = self.reactant_conc * 1e-9
            elif self.reactant_conc_type == "ng/min":
                self.reactant_MR = tools.permeation_rate_to_MR(
                    flow_rate=self.reactant_FR,
                    permeation_rate=self.reactant_conc,
                    reactant_gas=self.reactant_gas,
                )
            elif self.reactant_conc_type in ["Pa", "Torr", "bar", "mbar"]:
                self.reactant_MR = tools.vapor_pressure_to_MR(
                    vapor_pressure=self.reactant_conc,
                    P_units=self.reactant_conc_type,
                    system_pressure=self.P,
                    P_units_system=self.P_units,
                )
            if self.reactant_MR < 0 or self.reactant_MR > 1:
                raise ValueError(
                    "Issue calculating reactant mixing ratio."
                    "Mixing ratio must be between 0 and 1"
                )

            ### Reactant Diffusion Rate ###
            # Verify that the reactant diffusion rate is a number
            try:
                float(reactant_diffusion_rate)
            except ValueError:
                raise TypeError("Reactant diffusion rate must be a number")

            # Check if the user has previously manually inputted a diffusion rate
            try:
                self.manually_inputted_diffusion_rate  # noqa: B018
            except AttributeError:
                if not np.isnan(reactant_diffusion_rate):
                    self.manually_inputted_diffusion_rate = True
                else:
                    self.manually_inputted_diffusion_rate = False
                self.reactant_diffusion_rate = reactant_diffusion_rate
            else:
                if self.manually_inputted_diffusion_rate & ~np.isnan(
                    reactant_diffusion_rate
                ):
                    self.reactant_diffusion_rate = reactant_diffusion_rate

            ### Surface Area Density ###
            # Verify that the aerosol surface area is a number
            try:
                float(aerosol_surface_area)
            except ValueError:
                raise TypeError("Aerosol surface area must be a number")

            # Check if the user has previously manually inputted a surface area
            try:
                self.manually_inputted_surface_area  # noqa: B018
            except AttributeError:
                if not np.isnan(aerosol_surface_area):
                    self.manually_inputted_surface_area = True
                else:
                    self.manually_inputted_surface_area = False
                self.aerosol_surface_area = aerosol_surface_area
            else:
                if self.manually_inputted_surface_area & ~np.isnan(
                    aerosol_surface_area
                ):
                    self.aerosol_surface_area = aerosol_surface_area

            # Perform calculations for flows, carrier gas transport, and reactant diffusion
            self.flows(disp=disp)
            self.carrier_flow(disp=disp)
            self.aerosol_transport(disp=disp)
            self.reactant_diffusion(disp=disp)
        finally:
            # Turn flag off to allow __setattr__ to be used normally
            object.__setattr__(self, "_initializing", False)

    def flows(
        self,
        disp: bool = True,
    ) -> None:
        """Calculates Flow Tube flows.

        Args:
            disp (bool): Display calculated calculated values.

        Returns:
            None
        """
        ### Initialize Lists for displaying values ###
        var_names: list[str] = []
        var: list[float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        ### Calculate Cross Sectional Areas ###
        # Preinjector net cross section is the area between the FT wall
        # and the injector OD
        preinjector_net_cross_section = tools.cross_sectional_area(
            self.FT_ID
        ) - tools.cross_sectional_area(self.injector_OD)

        ### Flow Rates ###
        # Flow Rate Setpoints
        var_names += ["Reactant Flow Rate"]
        var += [self.reactant_FR]
        var_fmts += [".2f"]
        units += ["sccm"]
        var_names += ["Reactant Carrier Flow Rate"]
        var += [self.reactant_carrier_FR]
        var_fmts += [".1f"]
        units += ["sccm"]

        # Total Flow Rates
        total_reactant_FR = self.reactant_FR + self.reactant_carrier_FR
        self.total_FR = self.reactant_FR + self.reactant_carrier_FR + self.carrier_FR
        var_names += ["Total Reactant Flow Rate"]
        var += [total_reactant_FR]
        var_fmts += [".1f"]
        units += ["sccm"]

        ### Total Reactant Flow Velocity ###
        total_reactant_flow_velocity = flow_calc.sccm_to_velocity(
            self, total_reactant_FR, self.injector_ID
        )

        ### Minimum Carrier Flow Velocity & Rate ###
        # to prevent effect mentioned in Li et al., ACP, 2020
        min_carrier_flow_velocity = total_reactant_flow_velocity * 1.33
        min_carrier_FR = flow_calc.ccm_to_sccm(
            self, min_carrier_flow_velocity * preinjector_net_cross_section * 60
        )
        var_names += ["Minimum Carrier Flow Rate"]
        var += [min_carrier_FR]
        var_fmts += [".1f"]
        units += ["sccm"]
        if self.carrier_FR < min_carrier_FR:
            warnings.warn(
                "Carrier flow rate is below the minimum. "
                "This may affect the flow profile in the flow tube."
            )

        ### More Flow Rates ###
        var_names += ["Carrier Flow Rate"]
        var += [self.carrier_FR]
        var_fmts += [".1f"]
        units += ["sccm"]
        var_names += ["Total Flow Rate"]
        var += [self.total_FR]
        var_fmts += [".1f"]
        units += ["sccm"]

        ### Reactant Concentrations ###
        # Concentration inside of the injector (ppb)
        self.injector_conc = (
            self.reactant_FR / total_reactant_FR * self.reactant_MR * 1e9
        )
        var_names += [f"Injector {self.reactant_gas} Concentration"]
        var += [self.injector_conc]
        var_fmts += [".3g"]
        units += ["ppb"]

        # Concentration after the injector (ppb) - FT
        self.FT_conc = self.reactant_FR / self.total_FR * self.reactant_MR * 1e9
        self.FT_conc_molec = flow_calc.MR_to_molec(self, self.FT_conc)
        var_names += 2 * [f"Flow Tube {self.reactant_gas} Concentration"]
        var += [self.FT_conc, self.FT_conc_molec]
        var_fmts += [".3g", ".2e"]
        units += ["ppb", "molec. cm-3"]

        ### Flow Tube Flow Velocity ###
        self.flow_velocity = flow_calc.sccm_to_velocity(self, self.total_FR, self.FT_ID)
        var_names += ["Flow Tube Flow Velocity"]
        var += [self.flow_velocity]
        var_fmts += [".3g"]
        units += ["cm s-1"]

        ### Residence Times ###
        self.residence_time = self.FT_length / self.flow_velocity
        var_names += ["Flow Tube Residence Time"]
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
        disp: bool = True,
    ):
        """Performs and displays carrier gas transport calculations.

        Args:
            delta_T_radial (float): Radial temperature gradient (K).
            disp (bool): Display calculated values.

        Returns:
            None
        """

        ### Initialize Lists for displaying values ###
        var_names: list[str] = []
        var: list[float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        ### Carrier Gas Dynamic Viscosity (kg m-1 s-1) ###
        self.carrier_dynamic_viscosity = viscosity_density.dynamic_viscosity(
            self, self.carrier_gas
        )
        var_names += ["Carrier Gas Dynamic Viscosity"]
        var += [self.carrier_dynamic_viscosity]
        var_fmts += [".2e"]
        units += ["kg m-1 s-1"]

        ### Carrier Gas Density (kg m-3) ###
        self.carrier_density = viscosity_density.real_density(self, self.carrier_gas)
        var_names += ["Carrier Gas Density"]
        var += [self.carrier_density]
        var_fmts += [".3g"]
        units += ["kg m-3"]

        ### Reynolds Number - laminar flow if Re < 1800 ###
        self.Re = flow_calc.reynolds_number(self, self.total_FR, self.FT_ID)
        var_names += ["Flow Tube Reynolds Number"]
        var += [self.Re]
        var_fmts += [".0f"]
        units += ["unitless"]
        if self.Re > 1800:
            warnings.warn("Re > 1800. Flow in flow tube may not be laminar")

        ### Entrance length (cm) - see flow_calc.py for details ###
        length_to_laminar = flow_calc.length_to_laminar(self.FT_ID, self.Re)
        var_names += ["Flow Tube Entrance length"]
        var += [length_to_laminar]
        var_fmts += [".1f"]
        units += ["cm"]

        ### Pressure Gradient (%) - see flow_calc.py for details ###
        conductance = flow_calc.conductance(self, self.FT_ID, self.FT_length)
        pressure_gradient = flow_calc.pressure_gradient(
            self, conductance, self.total_FR
        )
        var_names += ["Flow Tube Pressure Gradient"]
        var += [pressure_gradient * 100]
        var_fmts += [".2f"]
        units += ["%"]

        ### Buoyancy Parameters - see flow_calc.py for details ###
        radial_buoyancy = flow_calc.buoyancy_parameters(
            self, self.radial_delta_T, self.FT_ID, self.Re
        )
        var_names += [f"Radial Buoyancy Parameter (ΔT={self.radial_delta_T:.1f} C)"]
        var += [radial_buoyancy]
        var_fmts += [".2f"]
        units += ["unitless"]
        if radial_buoyancy > 1:
            warnings.warn(
                "Radial buoyancy parameter > 1. "
                "Flow may be affected by buoyancy effects"
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

    def aerosol_transport(
        self,
        disp: bool = True,
    ) -> None:
        """Calculates the transport of aerosols in the flow tube.

        Args:
            disp (bool): Display calculated calculated values.

        Returns:
            None
        """
        ### Initialize Lists for displaying values ###
        var_names: list[str] = []
        var: list[float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        ### Aerosol Surface Area Density (SAD) (cm2 cm-3) ###
        # eq. 6 from Hanson and Kosciuch, 2003
        if self.aerosol_distribution == "monodisperse":
            calculated_aerosol_surface_area = (
                4
                * np.pi
                * (self.aerosol_diameter * 1e-7 / 2) ** 2
                * self.aerosol_number_conc
            )
        elif self.aerosol_distribution == "lognormal":
            calculated_aerosol_surface_area = (
                4
                * np.pi
                * (self.aerosol_diameter * 1e-7 / 2) ** 2
                * self.aerosol_number_conc
                * np.exp(2 * np.log(self.aerosol_sigma) ** 2)
            )
        var_names += ["Calculated Surface Area Density"]
        var += [calculated_aerosol_surface_area]  # pyright: ignore[reportPossiblyUnboundVariable]
        var_fmts += [".3g"]
        units += ["cm2 cm-3"]
        # Check if the user has previously manually inputted an SAD
        if self.manually_inputted_surface_area:
            if self.aerosol_surface_area < 0:
                raise ValueError("Aerosol surface area must be non-negative")
            var_names += [
                "Manually Inputted Surface Area Density \n(used in calculations)"
            ]
            var += [self.aerosol_surface_area]
            var_fmts += [".3g"]
            units += ["cm2 cm-3"]
        else:
            self.aerosol_surface_area = calculated_aerosol_surface_area  # pyright: ignore[reportPossiblyUnboundVariable]

        ### Surface Area Weighted Diameter (nm) ###
        # eq. 8 from Hanson and Kosciuch, 2003
        if self.aerosol_distribution == "lognormal":
            self.surface_area_weighted_diameter = self.aerosol_diameter * np.exp(
                2.5 * np.log(self.aerosol_sigma) ** 2
            )
            var_names += ["Surface Area Weighted Diameter"]
            var += [self.surface_area_weighted_diameter]
            var_fmts += [".3g"]
            units += ["nm"]

        ### Knudsen Number for carrier gas-aerosol interaction ###
        # - eq. 8 from Knopf et al., Anal. Chem., 2015
        Kn_carrier_aerosol = flow_calc.Kn(
            flow_calc.carrier_gas_mean_free_path(self), self.aerosol_diameter * 1e-7
        )
        var_names += ["Knudsen Number (carrier-aerosol)"]
        var += [Kn_carrier_aerosol]
        var_fmts += [".3g"]
        units += ["unitless"]

        ### Slip Correction Factor (unitless) ###
        # eq. 9.34 from Seinfeld and Pandis, 2016
        self.slip_correction = 1 + Kn_carrier_aerosol * (
            1.257 + 0.4 * np.exp(-1.1 / Kn_carrier_aerosol)
        )
        var_names += ["Slip Correction Factor"]
        var += [self.slip_correction]
        var_fmts += [".3g"]
        units += ["unitless"]

        ### Aerosol Diffusion Rate (cm2 s-1) ###
        # - eq. 9.73 from Seinfeld and Pandis, 2016
        self.aerosol_diffusion_rate = (
            tools.BOLTZMANN_CONSTANT
            * self.T_K
            / (
                3
                * np.pi
                * self.carrier_dynamic_viscosity
                * (self.aerosol_diameter * 1e-9)
            )
            * self.slip_correction
            * 100**2
        )
        var_names += ["Aerosol Diffusion Rate"]
        var += [self.aerosol_diffusion_rate]
        var_fmts += [".3g"]
        units += ["cm2 s-1"]

        ### Axial Distance, z* (unitless) ###
        # - eq. 2 from Knopf et al., Anal. Chem., 2015
        z_star = (
            self.axial_distance
            * np.pi
            / 2
            * self.aerosol_diffusion_rate
            / (flow_calc.sccm_to_ccm(self, self.total_FR) / 60)
        )
        var_names += ["Axial Distance, z*"]
        var += [z_star]
        var_fmts += [".3g"]
        units += ["unitless"]

        ### Effective Sherwood Number (unitless) ###
        # see flow_calc.py for details
        N_eff_Shw = flow_calc.N_eff_Shw(z_star=z_star)

        ### Thermal Particle Velocity (cm s-1) ###
        # eq. 9.87 from Seinfeld and Pandis, 2016
        mean_particle_mass = (
            self.aerosol_density
            * (4 / 3)
            * np.pi
            * (self.aerosol_diameter * 1e-7 / 2) ** 3
            / 1000
        )  # kg
        self.thermal_particle_velocity = 100 * np.sqrt(
            8 / np.pi * tools.BOLTZMANN_CONSTANT * self.T_K / mean_particle_mass
        )

        ### Mean Free Path (cm) ###
        # eq. 9.88 from Seinfeld and Pandis, 2016
        self.aerosol_mean_free_path = (
            2 * self.aerosol_diffusion_rate / self.thermal_particle_velocity
        )
        var_names += ["Mean Free Path"]
        var += [self.aerosol_mean_free_path]
        var_fmts += [".3g"]
        units += ["cm"]

        ### Knudsen Number (aerosol-wall) ###
        Kn_aerosol_wall = flow_calc.Kn(self.aerosol_mean_free_path, self.FT_ID)
        var_names += ["Knudsen Number (aerosol-wall)"]
        var += [Kn_aerosol_wall]
        var_fmts += [".3g"]
        units += ["unitless"]

        ### Tube Transmission (%) ###
        # eq. 21 from Knopf et al., Anal. Chem., 2015 - assuming the
        # sticking paramter, gamma is 1
        tube_transmission = np.exp(
            -1
            / (1 + 1 * 3 / (2 * N_eff_Shw * Kn_aerosol_wall))
            * self.thermal_particle_velocity
            / self.FT_ID
            * self.residence_time
        )
        var_names += ["Tube Transmission"]
        var += [tube_transmission * 100]
        var_fmts += [".1f"]
        units += ["%"]

        ### Display Values ###
        if disp:
            tools.table(
                "Aerosol Transport in Flow Tube",
                var_names,
                var,
                var_fmts,
                units,
            )

    def reactant_diffusion(
        self,
        disp: bool = True,
    ) -> None:
        """Performs and displays reactant diffusion calculations.

        Args:
            disp (bool): Display calculated calculated values.

        Returns:
            None
        """

        ### Initialize lists for displaying values ###
        var_names: list[str] = []
        var: list[float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        ### Reactant Diffusion Rate (cm2 s-1) ###
        # Calculate the diffusion rate if it was not manually inputted
        if self.manually_inputted_diffusion_rate:
            if self.reactant_diffusion_rate < 0:
                raise ValueError("Reactant diffusion rate must be non-negative")
            var_names += ["Manually Inputted Reactant Diffusion Rate"]
        else:
            if (
                self.reactant_gas not in diffusion_coef.sigmas
                and self.reactant_gas not in diffusion_coef.e_ks
            ):
                raise ValueError(
                    f"Must input reactant diffusion rate for {self.reactant_gas}"
                )
            else:
                self.reactant_diffusion_rate = (
                    diffusion_coef.binary_diffusion_coefficient(self)
                )
            var_names += ["Calculated Reactant Diffusion Rate \n(Lennard-Jones model)"]
        var += [self.reactant_diffusion_rate]
        var_fmts += [".3g"]
        units += ["cm2 s-1"]

        ### Thermal Molecular Velocity (cm s-1) ###
        # formula matched to values from Knopf et al., Anal. Chem., 2015
        self.reactant_molec_velocity = flow_calc.molec_velocity(
            self, float(mm.Formula(self.reactant_gas).mass)
        )

        ### Reactant Mean Free Path (cm) - Fuchs and Sutugin, 1971 ###
        self.reactant_mean_free_path = (
            3 * self.reactant_diffusion_rate / self.reactant_molec_velocity
        )

        ### Advection Rate (cm2 s-1) - eq. 1 from Knopf et al., Anal. Chem., 2015 ###
        advection_rate = self.flow_velocity * self.FT_ID
        var_names += ["Advection Rate"]
        var += [advection_rate]
        var_fmts += [".3g"]
        units += ["cm2 s-1"]

        ### Peclet Number - if > 10 then axial diffusion is negligible ###
        # - eq. 1 from Knopf et al., Anal. Chem., 2015
        self.Pe = advection_rate / self.reactant_diffusion_rate
        var_names += ["Peclet Number"]
        var += [self.Pe]
        var_fmts += [".4g"]
        units += ["unitless"]
        if self.Pe < 10:
            warnings.warn("Pe < 10. Axial diffusion is non-negligible")

        ### Mixing Time (s) - see flow_calc.py for details ###
        mixing_time = flow_calc.mixing_time(self, self.FT_ID)
        var_names += ["Mixing Time"]
        var += [mixing_time]
        var_fmts += [".2g"]
        units += ["s"]

        ### Mixing Length (cm) ###
        mixing_length = self.flow_velocity * mixing_time
        var_names += ["Mixing Length"]
        var += [mixing_length]
        var_fmts += [".2g"]
        units += ["cm"]

        ### Axial Distance ###
        # - eq. 2 from Knopf et al., Anal. Chem., 2015
        self.z_star = flow_calc.z_star(self, z=self.axial_distance, FR=self.total_FR)

        ### Effective Sherwood Number (unitless) ###
        # - eq. 11 from Knopf et al., Anal. Chem., 2015
        self.N_eff_Shw = flow_calc.N_eff_Shw(z_star=self.z_star)

        ### Knudsen Number for reactant-wall and reactant-aerosol interaction ###
        # - eq. 8 from Knopf et al., Anal. Chem., 2015
        self.Kn_wall = flow_calc.Kn(self.reactant_mean_free_path, self.FT_ID)
        # - and eq. 8 from  Hanson and Kosciuch, 2003
        if self.aerosol_distribution == "monodisperse":
            self.Kn_aerosol = flow_calc.Kn(
                self.reactant_mean_free_path, self.aerosol_diameter * 1e-7
            )
        elif self.aerosol_distribution == "lognormal":
            self.Kn_aerosol = flow_calc.Kn(
                self.reactant_mean_free_path, self.surface_area_weighted_diameter * 1e-7
            )

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
        hypothetical_gamma: ArrayLike | float,
        exposure_length: float = 1,
        gamma_wall: float = np.nan,
        disp: bool = True,
    ) -> None:
        """
        Calculates reactant uptake to aerosol and loss to flow tube
        walls.

        Args:
            hypothetical_gamma (ArrayLike or float): Hypothetical
                uptake coefficient to calculate diffusion correction
                factor.
            exposure_length (float): Length of the exposed surface in
                cm. Default is 1 cm.
            gamma_wall (float): Uptake coefficient for the wall
                (optional).
            disp (bool): Display calculated values.

        Returns:
            None.
        """
        ### Validate inputs ###
        hypothetical_gamma = input_validation.validate_reactant_uptake(
            obj=self,
            hypothetical_gamma=hypothetical_gamma,
            exposure_length=exposure_length,
            gamma_wall=gamma_wall,
        )

        ### Initialize lists for displaying values ###
        var_names: list[str] = []
        var: list[NDArray[np.float64] | float] = []
        var_fmts: list[str] = []
        units: list[str] = []

        ### Diffusion Resistance ###
        # Fuchs and Sutugin, 1971
        Gamma_diff = (self.Kn_aerosol * (1 + self.Kn_aerosol)) / (
            0.75 + 0.283 * self.Kn_aerosol
        )
        var_names += ["Diffusion Resistance (Γ_diff)"]
        var += [Gamma_diff]
        var_fmts += [".3g"]
        units += ["unitless"]

        ### Effective Uptake Coefficient ###
        gamma_eff = 1 / (1 / hypothetical_gamma + 1 / Gamma_diff)
        var_names += ["Effective Uptake Coefficient"]
        var += [gamma_eff]
        var_fmts += [".2e"]
        units += ["unitless"]

        ### Diffusion Correction Mangitude ###
        diffusion_correction = (hypothetical_gamma - gamma_eff) / hypothetical_gamma
        var_names += ["Diffusion Correction"]
        var += [diffusion_correction * 100]
        var_fmts += [".3g"]
        units += ["%"]

        ### Reactant-Particle Collision Rate (cm-3 s-1) ###
        # - eq. 6 from Hanson and Kosciuch, 2003
        k_c = self.aerosol_surface_area * self.reactant_molec_velocity / 4

        ### Observed Rate Contant ###
        # - eq. 7 from Hanson and Kosciuch, 2003
        self.k_obs = hypothetical_gamma * k_c / (1 + hypothetical_gamma / Gamma_diff)
        var_names += ["Observed Rate Constant (k_obs)"]
        var += [self.k_obs]
        var_fmts += [".3g"]
        units += ["s-1"]

        ### Reaction Rate Constant ###
        # standard equation, can be seen in Huynh and McNeill, J. Phys. 
        # Chem. A, 2021, for example
        self.k_rxn = hypothetical_gamma * k_c
        var_names += ["Reaction Rate Constant (k_rxn)"]
        var += [self.k_rxn]
        var_fmts += [".3g"]
        units += ["s-1"]

        ### Diffusion Rate Constant ###
        self.k_diff = 1 / (1 / self.k_obs - 1 / self.k_rxn)
        var_names += ["Diffusion Rate Constant (k_diff)"]
        var += [self.k_diff]
        var_fmts += [".3g"]
        units += ["s-1"]

        ### Reaction Time (s) ###
        reaction_time = exposure_length / self.flow_velocity
        var_names += [f"Reaction Time per {exposure_length:.1f} cm Exposure"]
        var += [reaction_time]
        var_fmts += [".3g"]
        units += ["s"]

        ### Loss to aerosol - see kinetics.py for details ###
        self.aerosol_loss = 1 - np.exp(
            -self.k_obs * exposure_length / self.flow_velocity
        )
        var_names += [f"Loss to Aerosol per {exposure_length:.1f} cm Exposure"]
        var += [self.aerosol_loss * 100]
        var_fmts += [".1f"]
        units += ["%"]

        # If wall loss if included, calculate the wall loss and total loss
        if ~np.isnan(gamma_wall):
            ### Wall Loss per Exposure Length ###
            # see kinetics.py for details
            wall_loss = kinetics.cylinder_loss(
                self,
                self.FT_ID,
                self.N_eff_Shw,
                self.Kn_wall,
                gamma_wall,
                exposure_length / self.flow_velocity,
            )
            self.k_wall = -np.log(1 - wall_loss) / (
                exposure_length / self.flow_velocity
            )
            var_names += [f"Wall Loss per {exposure_length:.1f} cm Exposure"]
            var += [wall_loss * 100]
            var_fmts += [".1f"]
            units += ["%"]

            ### Total Observed Rate Constant ###
            self.k_total = self.k_obs + self.k_wall
            var_names += ["Total Observed Rate Constant (k_obs + k_wall)"]
            var += [self.k_total]
            var_fmts += [".3g"]
            units += ["s-1"]

            ### Total Loss per Exposure Length ###
            self.total_loss = 1 - np.exp(
                -self.k_total * exposure_length / self.flow_velocity
            )
            var_names += [f"Total Loss per {exposure_length:.1f} cm Exposure"]
            var += [self.total_loss * 100]
            var_fmts += [".1f"]
            units += ["%"]

        ### Display Values ###
        if disp and not isinstance(hypothetical_gamma, np.ndarray):
            tools.table(
                "Reactant Uptake",
                var_names,
                var,  # pyright: ignore[reportArgumentType]
                var_fmts,
                units,
            )
