#!/usr/bin/env python3
"""
Example usage of flowtube package for boat reactor calculations.
"""

import flowtube


def main():
    """Example flow reactor calculation."""

    # Create a flow tube reactor instance
    boat = flowtube.BoatReactor(
        FT_ID=2.6,  # Flow tube inner diameter (cm)
        FT_length=100.0,  # Flow tube length (cm)
        injector_ID=1,  # Injector inner diameter (cm)
        injector_OD=1.2,  # Injector outer diameter (cm)
        reactant_gas="HCl",  # Reactant gas formula (options are Ar, He, Air, Br2, Cl2',
        # HBr, HCl, HI, H2O, I2, NO, N2, and O2)
        carrier_gas="N2",  # Carrier gas formula (options are Ar, He, N2, O2)
        reactant_MR=1e-6,  # Reactant mixing ratio (mol/mol)
        boat_width=2.0,  # Boat width (cm)
        boat_height=0.963,  # Boat height (cm)
        boat_length=53.8,  # Boat length (cm)
        boat_wall_thickness=0.11,  # Boat wall thickness (cm)
    )

    # Initialize with experimental conditions
    print("Initializing flow reactor with experimental conditions...")
    boat.initialize(
        reactant_FR=10.0,  # Reactant flow rate (sccm)
        reactant_carrier_FR=100,  # Reactant carrier flow rate (sccm)
        carrier_FR=1000,  # Carrier flow rate (sccm)
        P=50,  # Pressure
        P_units="Torr",  # Pressure units (options are Torr, bar, mbar, hPa, or Pa))
        T=25,  # Temperature (°C)
    )

    # Calculate uptake with a hypothetical uptake coefficient
    print("\nCalculating uptake with hypothetical gamma = 1E-7...")
    boat.reactant_uptake(hypothetical_gamma=1E-7)


if __name__ == "__main__":
    main()
