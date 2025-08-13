#!/usr/bin/env python3
"""
Example usage of flowTube package for flow reactor calculations.
"""

import sys
import os

from flowTube import flowTube

def main():
    """Example flow reactor calculation."""
    
    # Create a flow tube reactor instance
    reactor = flowTube.FT(
        FT_ID=2.5,              # Flow tube inner diameter (cm)
        FT_length=50.0,         # Flow tube length (cm)
        injector_ID=0.3,        # Injector inner diameter (cm)
        injector_OD=0.6,        # Injector outer diameter (cm)
        reactant_gas='HCl',     # Reactant gas formula
        carrier_gas='N2',       # Carrier gas formula
        reactant_MR=1e-6,       # Reactant mixing ratio (mol/mol)
        insert_ID=1.5,          # Insert inner diameter (cm) 
        insert_length=20.0      # Insert length (cm)
    )

    # Initialize with experimental conditions
    print("Initializing flow reactor with experimental conditions...")
    reactor.initialize(
        reactant_FR=10.0,       # Reactant flow rate (sccm)
        reactant_carrier_FR=100, # Reactant carrier flow rate (sccm)
        carrier_FR=1000,        # Carrier flow rate (sccm)
        P=50,                   # Pressure
        P_units='Torr',         # Pressure units
        T=25                    # Temperature (°C)
    )

    # Calculate uptake with a hypothetical uptake coefficient
    print("\nCalculating uptake with hypothetical gamma = 0.01...")
    correction_factor, loss_percent = reactor.reactant_uptake(
        hypothetical_gamma=0.01
    )

if __name__ == "__main__":
    main()
