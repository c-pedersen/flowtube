#!/usr/bin/env python3
"""
Example usage of flowtube package for coated wall reactor calculations.

If insert diameter and length are not specified, the entire flow tube will be treated as a coated wall reactor with 
the specified inner diameter and length. If they insert dimensions are specified, the insert will be treated as
the coated wall reactor with the specified inner diameter and length.
"""

import flowtube

def main():
    """Example flow reactor calculation."""
    
    # Create a flow tube reactor instance
    cwr = flowtube.CoatedWallReactor(
        FT_ID=2.6,              # Flow tube inner diameter (cm)
        FT_length=50.0,         # Flow tube length (cm)
        injector_ID=1,        # Injector inner diameter (cm)
        injector_OD=1.2,        # Injector outer diameter (cm)
        reactant_gas='HCl',     # Reactant gas formula (options are Ar, He, Air, Br2, Cl2', HBr, HCl, HI, H2O, I2, NO, 
                                # N2, and O2)
        carrier_gas='N2',       # Carrier gas formula (options are Ar, He, N2, O2)
        reactant_MR=1e-6,       # Reactant mixing ratio (mol/mol)
        insert_ID=1.5,          # OPTIONAL: Insert inner diameter (cm) 
        insert_length=20.0      # OPTIONAL: Insert length (cm)
    )

    # Initialize with experimental conditions
    print("Initializing flow reactor with experimental conditions...")
    cwr.initialize(
        reactant_FR=10.0,       # Reactant flow rate (sccm)
        reactant_carrier_FR=100, # Reactant carrier flow rate (sccm)
        carrier_FR=1000,        # Carrier flow rate (sccm)
        P=50,                   # Pressure
        P_units='Torr',         # Pressure units (options are Torr, bar, mbar, hPa, or Pa))
        T=25                    # Temperature (°C)
    )

    # Calculate uptake with a hypothetical uptake coefficient
    print("\nCalculating uptake with hypothetical gamma = 0.01...")
    _, _ = cwr.reactant_uptake(hypothetical_gamma=0.01)

if __name__ == "__main__":
    main()
