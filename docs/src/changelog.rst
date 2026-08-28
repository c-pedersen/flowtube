Changelog
=========

1.4.0
-------------------
- Removed functionality for calculating the diffusion corrected uptake coefficient
- Added support for changing individual parameters of the `boat_reactor` and `coated_wall_reactor` modules after initialization by automatically re-running the relevant calculations when a parameter is updated.
- Added exposure time return for kinetics fitting functions in `boat_reactor` and `coated_wall_reactor` modules.
- Added Lennard-Jones parameters for ClONO2, N2O5, O3, and NO2 to `diffusion_coef` module in order to support diffusion coefficient calculations for these species.
- Allowed for the use of a custom diffusion coefficient even if the species is in the database of the `diffusion_coef` module.
- Removed axial temperature gradient parameter and associated calculations.
- Updated minimum carrier flow rate calculation to account for the presence of an insert
- Added calculation of concentration inside the insert and used this concentration to calculate the fraction of unreacted surface sites after a given exposure time in the `coated_wall_reactor` module.
- Fixed issue where a manually inputted diffusion coefficient would be overwritten if any attributes were updated.
- Removed wall loss from `coated_wall_reactor` and moved it to optional for the `boat_reactor` module.
- Miscellaneous improvements and bug fixes.

1.3.1
-------------------
- Fixed critical error in vapor pressure to mixing ratio conversion function that caused incorrect results.
- Fixed critical issue with documentation build
- Added feature to `coated_wall_reactor` module to calculate fraction of unreacted surface sites after a given exposure time
- Added support for accounting for flow around the outside of inserts in the `coated_wall_reactor` module when computing insert flow velocity
- Change calculate_gamma to calculate_gamma_effective.
- Minor improvements

1.3.0
-------------------
- Improved documentation with usage examples and explanations of the API.
- Added a changelog to track changes and updates to the library.
- Added support for non-half-cylinder boats in the `boat_reactor` module.
- Added support for various reactant gas sources including gas cylinders, permeation tubes, and volatile sources.
- Improved testing
- Squashed some bugs and improved error handling in various modules.
