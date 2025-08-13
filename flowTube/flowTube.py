import math
import pandas as pd
from tabulate import tabulate
import warnings
import molmass as mm

from . import tools, diffusion_coef, viscosity_density

### Helper Functions ###
## Geometery
# Diameter to cross section
d_to_cs = lambda d: math.pi * (d / 2) ** 2

## Flow Calculations
# Flow rates - sccm to cm3 min-1
sccm_to_ccm = lambda self, FR: (tools.standard_P / self.P) * self.T / tools.standard_T_K * FR

# Flow rates - cm3 min-1 to sccm
ccm_to_sccm = lambda self, FR: (self.P / tools.standard_P) * tools.standard_T_K / self.T * FR

# Flow rate to flow velocity - sccm to cm s-1
sccm_to_velocity = lambda self, FR, cross_section: sccm_to_ccm(self, FR) / cross_section / 60

# Flow velocity to flow rate - cm s-1 to sccm
velocity_to_sccm = lambda self, velocity, cross_section: ccm_to_sccm(self, velocity * cross_section * 60)

# Residence time
residence_time = lambda flow_velocity, distance: distance / flow_velocity

## Display
# Print table of variables
table = lambda var_names, var, var_fmts, units: print(tabulate(pd.DataFrame(
    [var_names, [format(v, '8'+fmt) for v, fmt in zip(var, var_fmts)], units]).T, 
                                                               disable_numparse = True, tablefmt = 'fancy_grid', showindex=False))

class FT:
    def __init__(self, FT_ID, FT_length, injector_ID, injector_OD, reactant_gas, carrier_gas, reactant_MR, 
                 insert_ID=None, insert_length=0):
        '''
        Handles calculations relevant to flow rate, flow diagnostics, transport, and uptake for a coated wall reactor.
        By default assumes no insert. To calculate terms for a fully-coated cylindrical insert, simply pass values for the insert length
        and ID. See bottom of file for citations.
 
        Args:
            FT_ID (float): Inner diameter (cm) of flow tube.
            FT_length (float): Length (cm) of flow tube.
            injector_ID (float): Inner diameter (cm) of reactant injector.
            injector_OD (float): Outer diameter (cm) of reactant injector.
            reactant_gas (str): Molecular formula of reactant gas.
            carrier_gas (str): Molecular formula of carrier gas.
            reactant_MR (float): Reactant mixing ratio (mol mol-1).
            insert_ID (float, optional): Inner diameter (cm) of insert.
            insert_length (float, optional): Length (cm) of insert.
    
        Returns:
            None
        '''
        
        self.FT_ID = FT_ID
        self.FT_length = FT_length
        self.injector_ID = injector_ID
        self.injector_OD = injector_OD
        self.reactant_gas = reactant_gas
        self.carrier_gas = carrier_gas
        self.reactant_MR = reactant_MR
        self.insert_ID = insert_ID
        self.insert_length = insert_length

    def initialize(self, reactant_FR, reactant_carrier_FR, carrier_FR, P, P_units, T, disp=True):
        ''' Sets experimental conditions and calls calculation functions for numerous flow and diffusion parameters.
    
        Args:
            reactant_FR (float): Reactant flow rate (sccm).
            reactant_carrier_FR (float): Carrier flow rate (sccm) used to dilute the reactant.
            carrier_FR (float): Carrier flow rate (sccm) typically injected through a side port near the start of the flow tube.
            P (float): Pressure.
            P_units (str): Pressure units.
            T (float): Temperature (C).
            disp (bool): Display calculated calculated values.
    
        Returns:
            None
        '''

        self.P = tools.P_in_Pa(P, P_units)
        self.T = tools.T_in_K(T)

        self.flows(reactant_FR, reactant_carrier_FR, carrier_FR, disp=disp)
        self.carrier_flow(disp=disp)
        self.reactant_diffusion(disp=disp)
    
    def flows(self, reactant_FR, reactant_carrier_FR, carrier_FR, disp=True):
        ''' Calculates Flow Tube flows.
    
        Args:
            reactant_FR (float): Reactant flow rate (sccm).
            reactant_carrier_FR (float): Carrier flow rate (sccm) used to dilute the reactant.
            carrier_FR (float): Carrier flow rate (sccm) typically injected through a side port near the start of the flow tube.
            disp (bool): Display calculated calculated values.
    
        Returns:
            None
        '''
        
        # Lists for displaying values
        var_names = []; var = []; var_fmts = []; units = []

        # Flow Rate Setpoints
        var_names += ['Reactant Flow Rate']; var += [reactant_FR]; var_fmts += ['.2f']; units += ['sccm']
        var_names += ['Reactant Carrier Flow Rate']; var += [reactant_carrier_FR]; var_fmts += ['.1f']; units += ['sccm']
        
        # Total Flow Rates
        total_reactant_FR = reactant_FR + reactant_carrier_FR
        self.total_FR = reactant_FR + reactant_carrier_FR + carrier_FR
        var_names += ['Total Reactant Flow Rate']; var += [total_reactant_FR]; var_fmts += ['.1f']; units += ['sccm']
    
        # Total Reactant Flow Velocity
        total_reactant_flow_velocity = sccm_to_velocity(self, total_reactant_FR, d_to_cs(self.injector_ID))
        
        # Minimum Carrier Flow Velocity & Rate - to prevent effect mentioned in Li et al., ACP, 2020
        min_carrier_flow_velocity = total_reactant_flow_velocity * 1.33
        min_carrier_FR = velocity_to_sccm(self, min_carrier_flow_velocity, d_to_cs(self.FT_ID) - d_to_cs(self.injector_OD))
        var_names += ['Minimum Carrier Flow Rate']; var += [min_carrier_FR]; var_fmts += ['.1f']; units += ['sccm']
        if carrier_FR < min_carrier_FR:
            warnings.warn('Carrier flow rate is below the minimum')

        # More Flow Rates
        var_names += ['Carrier Flow Rate']; var += [carrier_FR]; var_fmts += ['.1f']; units += ['sccm']
        var_names += ['Total Flow Rate']; var += [self.total_FR]; var_fmts += ['.1f']; units += ['sccm']

        # Reactant Concentrations (ppb)
        injector_conc = reactant_FR / total_reactant_FR * self.reactant_MR * 1e9
        FT_conc = reactant_FR / self.total_FR * self.reactant_MR * 1e9
        var_names += [f'Injector {self.reactant_gas} Concentration']; var += [injector_conc]; var_fmts += ['.3g']; units += ['ppb']
        var_names += [f'Flow Tube {self.reactant_gas} Concentration']; var += [FT_conc]; var_fmts += ['.3g']; units += ['ppb']
        
        # Total Flow Velocity
        total_flow_velocity = sccm_to_velocity(self, self.total_FR, d_to_cs(self.FT_ID))
        var_names += ['Flow Tube Velocity']; var += [total_flow_velocity]; var_fmts += ['.3g']; units += ['cm s-1']
        
        # Residence Time
        self.FT_residence_time = residence_time(total_flow_velocity, self.FT_length)
        var_names += ['Flow Tube Residence Time']; var += [self.FT_residence_time]; var_fmts += ['.3g']; units += ['s']
    
        ### Display Values ###
        if disp == True:
            table(var_names, var, var_fmts, units)
    
    def carrier_flow(self, disp=True):
        ''' Performs and displays carrier gas transport calculations.
    
        Args:
            disp (bool): Display calculated values.
    
        Returns:
            None
        '''
        
        # Lists for displaying values
        var_names = []; var = []; var_fmts = []; units = []
    
        # Carrier Gas Dynamic Viscosity (kg m-1 s-1)
        carrier_dynamic_viscosity = viscosity_density.dynamic_viscosity(self.carrier_gas, self.T, self.P)
        var_names += ['Carrier Gas Dynamic Viscosity']; var += [carrier_dynamic_viscosity]; var_fmts += ['.2e']; units += ['kg m-1 s-1']
        
        # Carrier Gas Density (kg m-3) NOT DYNAMICALLY INTEGRATED
        carrier_density = viscosity_density.real_density(self.carrier_gas, self.T, self.P)
        var_names += ['Carrier Gas Density']; var += [carrier_density]; var_fmts += ['.3g']; units += ['kg m-3']
        
        # Reynolds Number - laminar flow if Re < 1800
        Re_FT = ((carrier_density / 100**3) * sccm_to_velocity(self, self.total_FR, d_to_cs(self.FT_ID)) * 
                 self.FT_ID / (carrier_dynamic_viscosity / 100))
        var_names += ['Flow Tube Reynolds Number']; var += [Re_FT]; var_fmts += ['.0f']; units += ['unitless']
        if Re_FT > 1800:
            warnings.warn('Re > 1800. Flow in flow tube may not be laminar')

        if self.insert_length > 0:
            Re_insert = ((carrier_density / 100**3) * sccm_to_velocity(self, self.total_FR, d_to_cs(self.insert_ID)) * 
                         self.insert_ID / (carrier_dynamic_viscosity / 100))
            var_names += ['Insert Reynolds Number']; var += [Re_insert]; var_fmts += ['.0f']; units += ['unitless']
            if Re_insert > 1800:
                warnings.warn('Re > 1800. Flow in insert may not be laminar')
        
        # Entrance length (cm) - length to achieve laminar profile – Bird et al., 2002, page 52
        length_to_laminar = 0.035 * self.FT_ID * Re_FT
        var_names += ['Flow Tube Entrance length']; var += [length_to_laminar]; var_fmts += ['.1f']; units += ['cm']
        
        # Conductance (L s-1) - eq. 3.17 from Moore et al., 2009
        FT_conductance = 32600 * self.FT_ID ** 4 / (carrier_dynamic_viscosity * 1e7 * 
                                                          (self.FT_length - self.insert_length)) * self.P / tools.P_CF['Torr']
        if self.insert_length > 0:
            insert_conductance = 32600 * self.insert_ID ** 4 / (carrier_dynamic_viscosity * 1e7 * 
                                                                self.insert_length) * self.P / tools.P_CF['Torr'] 
            total_conductance = 1 / (1 / insert_conductance + 1 / FT_conductance)

        
        # Pressure Gradient (%) - eqs. 3.9 & 3.10 from Moore et al., 2009
        if self.insert_length > 0:
            insert_pressure_gradient = sccm_to_ccm(self, self.total_FR) / 60 / 1000 / insert_conductance 
            var_names += ['Insert Pressure Gradient']; var += [insert_pressure_gradient * 100]; var_fmts += ['.2f']; units += ['%']
            
            total_pressure_gradient = sccm_to_ccm(self, self.total_FR) / 60 / 1000 / total_conductance 
            var_names += ['Total Pressure Gradient']; var += [total_pressure_gradient * 100]; var_fmts += ['.2f']; units += ['%']
            
        else:
            FT_pressure_gradient = sccm_to_ccm(self, self.total_FR) / 60 / 1000 / FT_conductance 
            var_names += ['Flow Tube Pressure Gradient']; var += [FT_pressure_gradient * 100]; var_fmts += ['.2f']; units += ['%']
    
        ### Display Values ###
        if disp == True:
            table(var_names, var, var_fmts, units)
    
    def reactant_diffusion(self, disp=True):
        ''' Performs and displays reactant diffusion calculations.
    
        Args:
            disp (bool): Display calculated calculated values.
    
        Returns:
            None
        '''

        # Lists for displaying values
        var_names = []; var = []; var_fmts = []; units = []

        # Reactant Diffusion Rate (cm2 s-1)
        reactant_diffusion_rate = diffusion_coef.binary_diffusion_coefficent(self.reactant_gas, self.carrier_gas, self.T, self.P)
        var_names += ['Reactant Diffusion Rate']; var += [reactant_diffusion_rate]; var_fmts += ['.3g']; units += ['cm2 s-1']
        
        # Advection Rate (cm2 s-1) - eq. 1 from Knopf et al., Anal. Chem., 2015
        if self.insert_length > 0:
            advection_rate = sccm_to_velocity(self, self.total_FR, d_to_cs(self.insert_ID)) * self.insert_ID
            var_names += ['Insert Advection Rate']
        else:
            advection_rate = sccm_to_velocity(self, self.total_FR, d_to_cs(self.FT_ID)) * self.FT_ID
            var_names += ['Flow Tube Advection Rate']
        var += [advection_rate]; var_fmts += ['.3g']; units += ['cm2 s-1']
        
        # Peclet Number - if > 10 then axial diffusion is negligible - eq. 1 from Knopf et al., Anal. Chem., 2015
        Pe = advection_rate / reactant_diffusion_rate
        var_names += ['Peclet Number']; var += [reactant_diffusion_rate]; var_fmts += ['.4g']; units += ['unitless']
        if Pe < 10:
            warnings.warn('Pe < 10. Axial diffusion is non-negligible')
        
        # Mixing Time (s) - Hanson and Lovejoy, Geophys. Res. Lett., 1994
        if self.insert_length > 0:
            mixing_time = (self.insert_ID / 2)**2 / (5 * reactant_diffusion_rate) 
            var_names += ['Insert Mixing Time']
        else:
            mixing_time = (self.FT_ID / 2)**2 / (5 * reactant_diffusion_rate) 
            var_names += ['Flow Tube Mixing Time']
        var += [mixing_time]; var_fmts += ['.2g']; units += ['s']
        
        # Mixing Length (cm) - Hanson and  Kosciuch, J. Phys. Chem. A, 2003
        if self.insert_length > 0:
            mixing_length = sccm_to_velocity(self, self.total_FR, d_to_cs(self.insert_ID)) * mixing_time
            var_names += ['Insert Mixing Length']
        else:
            mixing_length = sccm_to_velocity(self, self.total_FR, d_to_cs(self.FT_ID)) * mixing_time
            var_names += ['Flow Tube Mixing Length']
        var += [mixing_length]; var_fmts += ['.2g']; units += ['cm']
        
        # Axial Distance (unitless) - eq. 2 from Knopf et al., Anal. Chem., 2015
        z_star_FT = self.FT_length * math.pi / 2 * reactant_diffusion_rate / (sccm_to_ccm(self, self.total_FR) / 60)
        if self.insert_length > 0:
            z_star_insert = self.insert_length * math.pi / 2 * reactant_diffusion_rate / (sccm_to_ccm(self, self.total_FR) / 60)
        
        # Effective Sherwood Number (unitless) - eq. 11 from Knopf et al., Anal. Chem., 2015
        self.N_eff_Shw_FT = 3.6568 + 0.0978 / (z_star_FT + 0.0154)
        if self.insert_length > 0:
            self.N_eff_Shw_insert = 3.6568 + 0.0978 / (z_star_insert + 0.0154)
        
        # Thermal Molecular Velocity (cm s-1) - formula matched to values from Knopf et al., Anal. Chem., 2015
        self.reactant_molec_velocity = 100 * math.sqrt(8 / math.pi * tools.R * self.T / mm.Formula(self.reactant_gas).mass * 1000)
        
        # Reactant Mean Free Path (cm) - Fuchs and Sutugin, 1971
        reactant_mean_free_path = 3 * reactant_diffusion_rate / self.reactant_molec_velocity
        
        # Knudsen Number for reactant-wall/insert interaction - eq. 8 from Knopf et al., Anal. Chem., 2015
        self.Kn_FT =  2 * reactant_mean_free_path / self.FT_ID
        if self.insert_length > 0:
            self.Kn_insert = 2 * reactant_mean_free_path / self.insert_ID
        
        # Diffusion Limited Rate Constant (s-1) - eq. 10 from Knopf et al., Anal. Chem., 2015
        if self.insert_length > 0:
            k_diff = 4 * self.N_eff_Shw_insert * reactant_diffusion_rate / self.insert_ID ** 2
        else:
            k_diff = 4 * self.N_eff_Shw_FT * reactant_diffusion_rate / self.FT_ID ** 2
        var_names += ['Diffusion Limited Rate Constant']; var += [k_diff]; var_fmts += ['.3g']; units += ['s-1']
        
        # Diffusion Limited Effective Uptake Coefficient - eq. 19 from Knopf et al., Anal. Chem., 2015
        if self.insert_length > 0:
            gamma_eff_diff = self.insert_ID / self.reactant_molec_velocity * k_diff
        else:
            gamma_eff_diff = self.FT_ID / self.reactant_molec_velocity * k_diff
        var_names += ['Diffusion Limited Effective Uptake Coefficient']; var += [gamma_eff_diff]; var_fmts += ['.2g']
        units += ['unitless']
        
        # Diffusion Limited Uptake Coefficient – diffusion correction limit from Tang et al., Atmos. Chem. Phys., 2014.
        gamma_diff = gamma_eff_diff / 0.1
        var_names += ['Approx. Diffusion Limited Uptake Coefficient']; var += [gamma_diff]; var_fmts += ['.2g']; units += ['unitless']
        
        ### Display Values ###
        if disp == True:
            table(var_names, var, var_fmts, units)
    
    def reactant_uptake(self, hypothetical_gamma, disp=True):
        ''' Calculates reactant uptake to coated wall insert and loss to flow tube walls.
    
        Args:
            hypothetical_gamma (float): Hypothetical uptake coefficient to calculate diffusion correction factor.
            disp (bool): Display calculated values.
    
        Returns:
            float: Gas phase insert diffusion correction factor (unitless).
            float: Percent loss in insert (%).
        '''

        # Flag unphysical gammas
        if hypothetical_gamma > 1:
            warnings.warn('γ must be less than or equal to 1')

        # Lists for displaying values
        var_names = []; var = []; var_fmts = []; units = []
        
        ### HCl Uptake ###
        # Residence Time in Insert (s)
        if self.insert_length > 0:
            insert_residence_time = self.insert_length / sccm_to_velocity(self, self.total_FR, d_to_cs(self.insert_ID))
            var_names += ['Insert Residence Time']; var += [insert_residence_time]; var_fmts += ['.2g']; units += ['s']
    
        # Diffusion Correction Factor - gamma_eff / gamma - eq. 15 from Knopf et al., Anal. Chem., 2015
        if self.insert_length > 0:
            C_g = 1 / (1 + hypothetical_gamma * 3 / (2 * self.N_eff_Shw_insert * self.Kn_insert))
            var_names += ['Insert Diffusion Correction Factor (γ_eff/γ)']
        else:
            C_g = 1 / (1 + hypothetical_gamma * 3 / (2 * self.N_eff_Shw_FT * self.Kn_FT))
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
        
        # Observed Loss Rate (s-1) - eq. 19 from Knopf et al., Anal. Chem., 2015
        if self.insert_length > 0:
            k_obs = gamma_eff * self.reactant_molec_velocity / self.insert_ID
        else:
            k_obs = gamma_eff * self.reactant_molec_velocity / self.FT_ID
        var_names += ['Observed Loss Rate']; var += [k_obs]; var_fmts += ['.3g']; units += ['s-1']
        
        # Percent Loss
        if self.insert_length > 0:
            loss = 1 - math.exp(- gamma_eff * self.reactant_molec_velocity / self.insert_ID * insert_residence_time)
            var_names += ['Insert Loss']
        else:
            loss = 1 - math.exp(- gamma_eff * self.reactant_molec_velocity / self.FT_ID * self.FT_residence_time / 4)
            var_names += ['Flow Tube Loss - 1/4 Length']
        var += [loss * 100]; var_fmts += ['.1f']; units += ['%']
        
        # Reactant Transmission through Flow Tube Walls - calculated over entire FT minus the insert - eq. 21 from Knopf et al., 2015
        gamma_wall = 9e-6 # halocarbon-wax coated wall (760 Torr, 207 K) – 10x upper limit calculated from Huynh and McNeill, 2021.
        reactant_transmission = math.exp(-gamma_wall / (1 + gamma_wall * 3 / (2 * self.N_eff_Shw_FT * self.Kn_FT)) * 
                                         self.reactant_molec_velocity / self.FT_ID * self.FT_residence_time)
        var_names += ['Est. Reactant Transmission in Uncoated, Empty Flow Tube']; var += [reactant_transmission * 100]; 
        var_fmts += ['.2g']; units += ['%']
    
        ### Display Values ###
        if disp == True:
            table(var_names, var, var_fmts, units)
    
        return C_g, loss
    
    
    ''' Citations:
    
        Moore, J.H., Davis, C.C., Coplan, M.A., 2009. Building Scientific Apparatus, 4th ed. ed. Cambridge University Press, Leiden.
        
        Seinfeld, J.H., Pandis, S.N., 2016. Atmospheric chemistry and physics: from air pollution to climate change, Third edition. ed.
        John Wiley & Sons, Hoboken, New Jersey.
        
        Bird, R.B., Stewart, W.E., Lightfoot, E.N., 2002. Transport phenomena, 2nd, Wiley international ed ed. J. Wiley, New York.
        
        Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and Penetration of Gas Molecules and Aerosol Particles through
        Laminar Flow Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87, 3746–3754. https://doi.org/10.1021/ac5042395
        
        Hanson, D.R., Lovejoy, E.R., 1994. The uptake of N2O5 onto small sulfuric acid particles. Geophys. Res. Lett. 21, 2401–2404.
        https://doi.org/10.1029/94GL02288
        
        Hanson, D., Kosciuch, E., 2003. The NH3 Mass Accommodation Coefficient for Uptake onto Sulfuric Acid Solutions. J. Phys. Chem.
        A 107, 2199–2208. https://doi.org/10.1021/jp021570j
        
        Fuchs, N.A., Sutugin, A.G., 1971. HIGH-DISPERSED AEROSOLS, in: Hidy, G.M., Brock, J.R. (Eds.), Topics in Current Aerosol
        Research, International Reviews in Aerosol Physics and Chemistry. Pergamon, p. 1. 
        https://doi.org/10.1016/B978-0-08-016674-2.50006-6
    
        Huynh, H.N., McNeill, V.F., 2021. Heterogeneous Reactivity of HCl on CaCO3 Aerosols at Stratospheric Temperature. ACS Earth
        Space Chem. 5, 1896–1901. https://doi.org/10.1021/acsearthspacechem.1c00151

        Tang, M.J., Cox, R.A., Kalberer, M., 2014. Compilation and evaluation of gas phase diffusion coefficients of reactive trace
        gases in the atmosphere: volume 1. Inorganic compounds. Atmos. Chem. Phys. 14, 9233–9247. 
        https://doi.org/10.5194/acp-14-9233-2014

    '''
