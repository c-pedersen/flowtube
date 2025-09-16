import numpy as np
import molmass as mm # pyright: ignore[reportMissingImports]
from numpy.typing import NDArray
import warnings

from . import tools, diffusion_coef, viscosity_density, flow_calc

class CoatedWallReactor:
    def __init__(self, FT_ID: float, FT_length: float, injector_ID: float, injector_OD: float, reactant_gas: str, 
                 carrier_gas: str, reactant_MR: float, insert_ID: float = 0, insert_length: float = 0) -> None:
        '''
        Handles calculations relevant to flow rate, flow diagnostics, transport, and uptake for a coated wall reactor.
        By default assumes no insert and a fully coated flow tube. To calculate terms for a fully-coated cylindrical 
        insert or for a partially coated flow tube, simply pass values for the insert length and ID.
 
        Args:
            FT_ID (float): Inner diameter (cm) of flow tube.
            FT_length (float): Length (cm) of flow tube.
            injector_ID (float): Inner diameter (cm) of reactant injector.
            injector_OD (float): Outer diameter (cm) of reactant injector.
            reactant_gas (str): Molecular formula of reactant gas (supported Ar, He, Air, Br2, Cl2, HBr, HCl, HI, H2O, 
                I2, NO, N2, and O2).
            carrier_gas (str): Molecular formula of carrier gas (supported: Ar, He, N2, O2).
            reactant_MR (float): Reactant mixing ratio (mol mol-1).
            insert_ID (float, optional): Inner diameter (cm) of insert.
            insert_length (float, optional): Length (cm) of insert.
    
        Returns:
            None
        '''

        ### Check for valid inputs ###
        # Check if the gases are supported
        if reactant_gas not in diffusion_coef.sigmas.keys():
            raise ValueError(f"Unsupported reactant gas. Supported gases: {', '.join(diffusion_coef.sigmas.keys())}")
        if carrier_gas not in viscosity_density.a.keys():
            raise ValueError(f"Unsupported carrier gas. Supported gases: {', '.join(viscosity_density.a.keys())}")
        
        # Check physicality of insert dimensions
        if insert_ID < 0 or insert_length < 0:
            raise ValueError('Insert ID and length must be positive')
        elif insert_ID > FT_ID:
            raise ValueError('Insert ID cannot be larger than flow tube ID')
        elif insert_length > FT_length:
            raise ValueError('Insert length cannot be larger than flow tube length')
        
        # Check physicality of injector dimensions
        if injector_ID < 0 or injector_OD < 0:
            raise ValueError('Injector ID and OD must be positive')
        elif injector_ID > FT_ID:
            raise ValueError('Injector ID cannot be larger than flow tube ID')
        elif injector_OD > FT_ID:
            raise ValueError('Injector OD cannot be larger than flow tube ID')
        elif injector_ID > injector_OD:
            raise ValueError('Injector ID cannot be larger than injector OD')
        
        # Check mixing ratio
        if reactant_MR < 0 or reactant_MR > 1:
            raise ValueError('Reactant mixing ratio must be between 0 and 1')

        # Initialize variables
        self.FT_ID = FT_ID
        self.FT_length = FT_length
        self.injector_ID = injector_ID
        self.injector_OD = injector_OD
        self.reactant_gas = reactant_gas
        self.carrier_gas = carrier_gas
        self.reactant_MR = reactant_MR
        self.insert_ID = insert_ID
        self.insert_length = insert_length

    def initialize(self, reactant_FR: float, reactant_carrier_FR: float, carrier_FR: float, P: float, P_units: str, 
                   T: float, radial_delta_T: float = 1, axial_delta_T: float = 1, disp: bool = True) -> None:
        ''' Sets experimental conditions and calls calculation functions for numerous flow and diffusion parameters.
    
        Args:
            reactant_FR (float): Reactant flow rate (sccm).
            reactant_carrier_FR (float): Carrier flow rate (sccm) used to dilute the reactant.
            carrier_FR (float): Carrier flow rate (sccm) typically injected near the start of the flow tube.
            P (float): Pressure.
            P_units (str): Pressure units.
            T (float): Temperature (C).
            radial_delta_T (float): Radial temperature gradient (K) (default = 1 K).
            axial_delta_T (float): Axial temperature gradient (K) (default = 1 K).
            disp (bool): Display calculated calculated values.
    
        Returns:
            None
        '''
        ### Check for valid inputs ###
        # Check if flow rates are positive
        if reactant_FR < 0 or reactant_carrier_FR < 0 or carrier_FR < 0:
            raise ValueError('Flow rates must be positive')

        # Check if the pressure units are supported
        if P_units not in tools.P_CF.keys():
            raise ValueError(f"Unsupported pressure units. Supported units: {', '.join(tools.P_CF.keys())}")
        elif P < 0:
            raise ValueError('Pressure must be positive')
        
        # Check if the temperature & temperature gradients are valid numbers
        if T < -273.15:
            raise ValueError('Temperature must be above absolute zero (-273.15 C)')
        if radial_delta_T < 0 or axial_delta_T < 0:
            raise ValueError('Temperature gradients must be positive')

        self.P = tools.P_in_Pa(P, P_units)
        self.T = tools.T_in_K(T)

        self.flows(reactant_FR, reactant_carrier_FR, carrier_FR, disp=disp)
        self.carrier_flow(radial_delta_T=radial_delta_T, axial_delta_T=axial_delta_T, disp=disp)
        self.reactant_diffusion(disp=disp)

    def flows(self, reactant_FR: float, reactant_carrier_FR: float, carrier_FR: float, disp: bool = True) -> None:
        ''' Calculates Flow Tube flows.
    
        Args:
            reactant_FR (float): Reactant flow rate (sccm).
            reactant_carrier_FR (float): Carrier flow rate (sccm) used to dilute the reactant.
            carrier_FR (float): Carrier flow rate (sccm) typically injected near the start of the flow tube.
            disp (bool): Display calculated calculated values.
    
        Returns:
            None
        '''
        # Check if the flow rates are positive
        if reactant_FR < 0 or reactant_carrier_FR < 0 or carrier_FR < 0:
            raise ValueError('Flow rates must be positive')
        
        # Lists for displaying values
        var_names: list[str] = []; var: list[float] = []; var_fmts: list[str] = []; units: list[str] = []

        # Flow Rate Setpoints
        var_names += ['Reactant Flow Rate']; var += [reactant_FR]; var_fmts += ['.2f']; units += ['sccm']
        var_names += ['Reactant Carrier Flow Rate']; var += [reactant_carrier_FR]; var_fmts += ['.1f']; 
        units += ['sccm']
        
        # Total Flow Rates
        total_reactant_FR = reactant_FR + reactant_carrier_FR
        self.total_FR = reactant_FR + reactant_carrier_FR + carrier_FR
        var_names += ['Total Reactant Flow Rate']; var += [total_reactant_FR]; var_fmts += ['.1f']; units += ['sccm']
    
        # Total Reactant Flow Velocity
        total_reactant_flow_velocity = flow_calc.sccm_to_velocity(self, total_reactant_FR, self.injector_ID)

        # Minimum Carrier Flow Velocity & Rate - to prevent effect mentioned in Li et al., ACP, 2020
        min_carrier_flow_velocity = total_reactant_flow_velocity * 1.33
        min_carrier_FR = flow_calc.ccm_to_sccm(self, min_carrier_flow_velocity * 
                                               (tools.cross_sectional_area(self.FT_ID) - 
                                                tools.cross_sectional_area(self.injector_OD)) * 60)
        var_names += ['Minimum Carrier Flow Rate']; var += [min_carrier_FR]; var_fmts += ['.1f']; units += ['sccm']
        if carrier_FR < min_carrier_FR:
            warnings.warn('Carrier flow rate is below the minimum. This may affect the flow profile in the flow tube.')

        # More Flow Rates
        var_names += ['Carrier Flow Rate']; var += [carrier_FR]; var_fmts += ['.1f']; units += ['sccm']
        var_names += ['Total Flow Rate']; var += [self.total_FR]; var_fmts += ['.1f']; units += ['sccm']

        # Reactant Concentrations (ppb)
        injector_conc = reactant_FR / total_reactant_FR * self.reactant_MR * 1e9
        FT_conc = reactant_FR / self.total_FR * self.reactant_MR * 1e9
        FT_conc_molec = flow_calc.MR_to_molec(self, FT_conc)
        var_names += [f'Injector {self.reactant_gas} Concentration']; var += [injector_conc]; var_fmts += ['.3g'] 
        units += ['ppb']
        var_names += [f'Flow Tube {self.reactant_gas} Concentration']; var += [FT_conc]; var_fmts += ['.3g']; 
        units += ['ppb']
        var_names += [f'Flow Tube {self.reactant_gas} Concentration']; var += [FT_conc_molec]; var_fmts += ['.2e']
        units += ['molec. cm-3']
        
        # Total Flow Velocities
        self.FT_flow_velocity = flow_calc.sccm_to_velocity(self, self.total_FR, self.FT_ID)
        if self.insert_length > 0:
            self.insert_flow_velocity = flow_calc.sccm_to_velocity(self, self.total_FR, self.insert_ID)
        var_names += ['Flow Tube Velocity']; var += [self.FT_flow_velocity]; var_fmts += ['.3g']; units += ['cm s-1']
        
        # Residence Time
        self.FT_residence_time = self.FT_length / self.FT_flow_velocity
        var_names += ['Flow Tube Residence Time']; var += [self.FT_residence_time]; var_fmts += ['.3g']; units += ['s']
    
        ### Display Values ###
        if disp == True:
            tools.table(var_names, var, var_fmts, units)
    
    def carrier_flow(self, radial_delta_T: float = 1, axial_delta_T: float = 1, disp: bool = True):
        ''' Performs and displays carrier gas transport calculations.
    
        Args:
            delta_T_radial (float): Radial temperature gradient (K).
            delta_T_axial (float): Axial temperature gradient (K).
            disp (bool): Display calculated values.
    
        Returns:
            None
        '''
        
        # Lists for displaying values
        var_names: list[str] = []; var: list[float] = []; var_fmts: list[str] = []; units: list[str] = []
    
        # Carrier Gas Dynamic Viscosity (kg m-1 s-1)
        self.carrier_dynamic_viscosity = viscosity_density.dynamic_viscosity(self, self.carrier_gas)
        var_names += ['Carrier Gas Dynamic Viscosity']; var += [self.carrier_dynamic_viscosity]; var_fmts += ['.2e'] 
        units += ['kg m-1 s-1']
        
        # Carrier Gas Density (kg m-3)
        self.carrier_density = viscosity_density.real_density(self, self.carrier_gas)
        var_names += ['Carrier Gas Density']; var += [self.carrier_density]; var_fmts += ['.3g']; units += ['kg m-3']
        
        # Reynolds Number - laminar flow if Re < 1800
        Re_FT = flow_calc.reynolds_number(self, self.total_FR, self.FT_ID)
        var_names += ['Flow Tube Reynolds Number']; var += [Re_FT]; var_fmts += ['.0f']; units += ['unitless']
        if Re_FT > 1800 :
            warnings.warn('Re > 1800. Flow in flow tube may not be laminar')

        if self.insert_length > 0:
            Re_insert = flow_calc.reynolds_number(self, self.total_FR, self.insert_ID)
            var_names += ['Insert Reynolds Number']; var += [Re_insert]; var_fmts += ['.0f']; units += ['unitless']
            if Re_insert > 1800:
                warnings.warn('Re > 1800. Flow in insert may not be laminar')
        
        # Entrance length (cm) - see flow_calc.py for details
        length_to_laminar = flow_calc.length_to_laminar(self.FT_ID, Re_FT)
        var_names += ['Flow Tube Entrance length']; var += [length_to_laminar]; var_fmts += ['.1f']; units += ['cm']
        
        # Pressure Gradient (%) - see flow_calc.py for details
        FT_conductance = flow_calc.conductance(self, self.FT_ID, self.FT_length)
        if self.insert_length > 0:
            insert_conductance = flow_calc.conductance(self, self.insert_ID, self.FT_length)
            total_conductance = 1 / (1 / insert_conductance + 1 / FT_conductance)

            insert_pressure_gradient = flow_calc.pressure_gradient(self, insert_conductance, self.total_FR)
            var_names += ['Insert Pressure Gradient']; var += [insert_pressure_gradient * 100]; var_fmts += ['.2f'] 
            units += ['%']
            
            total_pressure_gradient = flow_calc.pressure_gradient(self, total_conductance, self.total_FR)
            var_names += ['Total Pressure Gradient']; var += [total_pressure_gradient * 100]; var_fmts += ['.2f'] 
            units += ['%']
        else:
            FT_pressure_gradient = flow_calc.pressure_gradient(self, FT_conductance, self.total_FR)
            var_names += ['Flow Tube Pressure Gradient']; var += [FT_pressure_gradient * 100]; var_fmts += ['.2f'] 
            units += ['%']

        # Buoyancy Parameters - see flow_calc.py for details
        radial_buoyancy = flow_calc.buoyancy_parameters(self, radial_delta_T, self.FT_ID, Re_FT)
        axial_buoyancy = flow_calc.buoyancy_parameters(self, axial_delta_T, self.FT_length, Re_FT)
        var_names += [f'Radial Buoyancy Parameter (ΔT={radial_delta_T:.1f} C)']; var += [radial_buoyancy]
        var_fmts += ['.2f']; units += ['unitless']
        var_names += [f'Axial Buoyancy Parameter (ΔT={axial_delta_T:.1f} C)']; var += [axial_buoyancy]
        var_fmts += ['.2f']; units += ['unitless']
        if radial_buoyancy > 1:
            warnings.warn('Radial buoyancy parameter > 1. Flow may be affected by buoyancy effects')
        if axial_buoyancy > 1:
            warnings.warn('Axial buoyancy parameter > 1. Flow may be affected by buoyancy effects')
    
        ### Display Values ###
        if disp == True:
            tools.table(var_names, var, var_fmts, units)
    
    def reactant_diffusion(self, disp: bool = True):
        ''' Performs and displays reactant diffusion calculations.
    
        Args:
            disp (bool): Display calculated calculated values.
    
        Returns:
            None
        '''

        # Lists for displaying values
        var_names: list[str] = []; var: list[float] = []; var_fmts: list[str] = []; units: list[str] = []

        # Reactant Diffusion Rate (cm2 s-1)
        self.reactant_diffusion_rate = diffusion_coef.binary_diffusion_coefficent(self)
        var_names += ['Reactant Diffusion Rate']; var += [self.reactant_diffusion_rate]; var_fmts += ['.3g'] 
        units += ['cm2 s-1']

        # Thermal Molecular Velocity (cm s-1) - formula matched to values from Knopf et al., Anal. Chem., 2015
        self.reactant_molec_velocity = flow_calc.molec_velocity(
            self, float(mm.Formula(self.reactant_gas).mass)) # pyright: ignore[reportUnknownArgumentType, reportUnknownMemberType]
        
        # Reactant Mean Free Path (cm) - Fuchs and Sutugin, 1971
        reactant_mean_free_path = 3 * self.reactant_diffusion_rate / self.reactant_molec_velocity
        
        # Advection Rate (cm2 s-1) - eq. 1 from Knopf et al., Anal. Chem., 2015
        if self.insert_length > 0:
            advection_rate = self.insert_flow_velocity * self.insert_ID
            var_names += ['Insert Advection Rate']
        else:
            advection_rate = self.FT_flow_velocity * self.FT_ID
            var_names += ['Flow Tube Advection Rate']
        var += [advection_rate]; var_fmts += ['.3g']; units += ['cm2 s-1']
        
        # Peclet Number - if > 10 then axial diffusion is negligible - eq. 1 from Knopf et al., Anal. Chem., 2015
        Pe = advection_rate / self.reactant_diffusion_rate
        var_names += ['Peclet Number']; var += [self.reactant_diffusion_rate]; var_fmts += ['.4g']; units += ['unitless']
        if Pe < 10:
            warnings.warn('Pe < 10. Axial diffusion is non-negligible')
        
        # Mixing Time (s) - see flow_calc.py for details
        if self.insert_length > 0:
            mixing_time = flow_calc.mixing_time(self, self.insert_ID)
            var_names += ['Insert Mixing Time']
        else:
            mixing_time = flow_calc.mixing_time(self, self.FT_ID)
            var_names += ['Flow Tube Mixing Time']
        var += [mixing_time]; var_fmts += ['.2g']; units += ['s']
        
        # Mixing Length (cm)
        if self.insert_length > 0:
            mixing_length = self.insert_flow_velocity * mixing_time
            var_names += ['Insert Mixing Length']
        else:
            mixing_length = self.FT_flow_velocity * mixing_time
            var_names += ['Flow Tube Mixing Length']
        var += [mixing_length]; var_fmts += ['.2g']; units += ['cm']
        
        # Effective Sherwood Number (unitless) - eq. 11 from Knopf et al., Anal. Chem., 2015
        self.N_eff_Shw_FT = flow_calc.N_eff_Shw(self, self.FT_length, self.total_FR)
        if self.insert_length > 0:
            self.N_eff_Shw_insert = flow_calc.N_eff_Shw(self, self.insert_length, self.total_FR)
        
        # Knudsen Number for reactant-wall/insert interaction - eq. 8 from Knopf et al., Anal. Chem., 2015
        self.Kn_FT = flow_calc.Kn(reactant_mean_free_path, self.FT_ID)
        if self.insert_length > 0:
            self.Kn_insert = flow_calc.Kn(reactant_mean_free_path, self.insert_ID)
        
        # Diffusion Limited Rate Constant (s-1) and Uptake Coefficient - see flow_calc.py for details
        if self.insert_length > 0:
            k_diff = flow_calc.diffusion_limited_rate_constant(self, self.N_eff_Shw_insert, self.insert_ID)
            gamma_eff_diff = flow_calc.diffusion_limited_uptake_coefficient(self, self.insert_ID, k_diff)
        else:
            k_diff = flow_calc.diffusion_limited_rate_constant(self, self.N_eff_Shw_FT, self.FT_ID)
            gamma_eff_diff = flow_calc.diffusion_limited_uptake_coefficient(self, self.FT_ID, k_diff)
        var_names += ['Diffusion Limited Rate Constant']; var += [k_diff]; var_fmts += ['.3g']; units += ['s-1']
        var_names += ['Diffusion Limited Effective Uptake Coefficient']; var += [gamma_eff_diff]; var_fmts += ['.2g']
        units += ['unitless']
        
        # Diffusion Limited Uptake Coefficient – diffusion correction limit from Tang et al., Atmos. Chem. Phys., 2014.
        gamma_diff = gamma_eff_diff / 0.1
        var_names += ['Approx. Diffusion Limited Uptake Coefficient']; var += [gamma_diff]; var_fmts += ['.2g']
        units += ['unitless']
        
        ### Display Values ###
        if disp == True:
            tools.table(var_names, var, var_fmts, units)

    def reactant_uptake(self, hypothetical_gamma: NDArray[np.float64] | float, gamma_wall: float = 5e-6, 
                        disp: bool = True) -> tuple[NDArray[np.float64] | float, NDArray[np.float64] | float]:
        ''' Calculates reactant uptake to coated wall or insert and loss to flow tube walls.
    
        Args:
            hypothetical_gamma (float or numpy.ndarray): Hypothetical uptake coefficient to calculate diffusion 
                correction factor.
            gamma_wall (float): Wall uptake coefficient (default: 5e-6 for halocarbon wax coating - 
                Ivanov et al., J. Mass Spectrom., 2021).
            disp (bool): Display calculated values.
    
        Returns:
            float: Gas phase diffusion correction factor (unitless).
            float: Loss to insert or coated wall as ratio of inital amount (fraction).
        '''

        ### Check for valid inputs ###
        # Check if hypothetical_gamma is between 0 and 1
        if np.min(hypothetical_gamma) < 0 or np.max(hypothetical_gamma) > 1: # pyright: ignore[reportUnknownMemberType]
            raise ValueError('Hypothetical gamma must be between 0 and 1')
        
        # Check if gamma_wall is between 0 and 1
        if gamma_wall < 0 or gamma_wall > 1:
            raise ValueError('Wall uptake coefficient must be between 0 and 1')

        # Lists for displaying values
        var_names: list[str] = []; var: list[NDArray[np.float64] | float] = []; var_fmts: list[str] = []; 
        units: list[str] = []

        # Surface area of coated area
        if self.insert_length > 0:
            surface_area = 2 * np.pi * self.insert_ID * self.insert_length / 4
            var_names += ['Insert surface area']
        else:
            surface_area = 2 * np.pi * self.FT_ID * self.FT_length / 4
            var_names += ['Coated wall surface area (1/4 length)']
        var += [surface_area]; var_fmts += ['.1f']; units += ['cm2'] 
    
        # Diffusion Correction Factor - gamma_eff / gamma - eq. 15 from Knopf et al., Anal. Chem., 2015
        if self.insert_length > 0:
            C_g = flow_calc.correction_factor(self.N_eff_Shw_insert, self.Kn_insert, hypothetical_gamma)
            var_names += ['Insert Diffusion Correction Factor (γ_eff/γ)']
        else:
            C_g = flow_calc.correction_factor(self.N_eff_Shw_FT, self.Kn_FT, hypothetical_gamma)
            var_names += ['Flow Tube Diffusion Correction Factor (γ_eff/γ)']
        var += [C_g]; var_fmts += ['.3g']; units += ['unitless']
    
        # Diffusion Correction
        diff_corr = 1 - C_g
        if self.insert_length > 0:
            var_names += ['Insert Diffusion Correction']
        else:
            var_names += ['Flow Tube Diffusion Correction']
        var += [diff_corr * 100]; var_fmts += ['.1f']; units += ['%']
        
        # Effective Uptake Coefficient - eq. 15 from Knopf et al., Anal. Chem., 2015
        gamma_eff = hypothetical_gamma * C_g
        var_names += ['Effective Uptake Coefficient']; var += [gamma_eff]; var_fmts += ['.2e']; units += ['unitless']
        
        # Observed Loss Rate (s-1) - see flow_calc.py for details
        if self.insert_length > 0:
            k_obs = flow_calc.observed_loss_rate(self, self.insert_ID, gamma_eff)
        else:
            k_obs = flow_calc.observed_loss_rate(self, self.FT_ID, gamma_eff)
        var_names += ['Observed Loss Rate']; var += [k_obs]; var_fmts += ['.3g']; units += ['s-1']
        
        # Uptake to coated region - see flow_calc.py for details
        if self.insert_length > 0:
            uptake = flow_calc.cylinder_loss(self, self.insert_ID, self.N_eff_Shw_insert, self.Kn_insert, 
                                             hypothetical_gamma, self.insert_length / self.insert_flow_velocity)
            var_names += ['Insert Loss']
        else:
            uptake = flow_calc.cylinder_loss(self, self.FT_ID, self.N_eff_Shw_FT, self.Kn_FT, hypothetical_gamma, 
                                             self.FT_residence_time / 4)
            var_names += ['Flow Tube Loss - 1/4 Length']
        var += [uptake * 100]; var_fmts += ['.1f']; units += ['%']
        
        # Reactant Wall Loss (entire FT minus insert)
        reactant_wall_loss = flow_calc.cylinder_loss(self, self.FT_ID, self.N_eff_Shw_FT, self.Kn_FT, gamma_wall, 
                                                     self.FT_residence_time * (1 - self.insert_length / self.FT_length))
        var_names += ['Estimated Wall Loss']; var += [reactant_wall_loss * 100]; var_fmts += ['.2g']; units += ['%']
    
        ### Display Values ###
        if disp is True or not isinstance(var, np.ndarray):
            tools.table(var_names, var, var_fmts, units)  # pyright: ignore[reportArgumentType]
    
        return C_g, uptake
    
''' Citations:
    
    Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and Penetration of Gas Molecules and Aerosol Particles 
    through Laminar Flow Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87, 3746–3754. 
    https://doi.org/10.1021/ac5042395
    
    Hanson, D., Kosciuch, E., 2003. The NH3 Mass Accommodation Coefficient for Uptake onto Sulfuric Acid Solutions. 
    J. Phys. Chem. A 107, 2199–2208. https://doi.org/10.1021/jp021570j
    
    Fuchs, N.A., Sutugin, A.G., 1971. HIGH-DISPERSED AEROSOLS, in: Hidy, G.M., Brock, J.R. (Eds.), Topics in Current 
    Aerosol Research, International Reviews in Aerosol Physics and Chemistry. Pergamon, p. 1. 
    https://doi.org/10.1016/B978-0-08-016674-2.50006-6

    Ivanov, A.V., Molina, M.J., Park, J., 2021. Experimental study on HCl uptake by MgCl2 and sea salt under humid 
    conditions. J Mass Spectrom 56, e4601. https://doi.org/10.1002/jms.4601

    Tang, M.J., Cox, R.A., Kalberer, M., 2014. Compilation and evaluation of gas phase diffusion coefficients of 
    reactive trace gases in the atmosphere: volume 1. Inorganic compounds. Atmos. Chem. Phys. 14, 9233–9247. 
    https://doi.org/10.5194/acp-14-9233-2014
'''