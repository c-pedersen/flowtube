Changelog
=========

1.3.2
-------------------
- Added exposure time return for kinetics fitting functions in `boat_reactor` and `coated_wall_reactor` modules.
- Added Lennard-Jones parameters for ClONO2, N2O5, O3, and NO2 to `diffusion_coef` module in order to support diffusion coefficient calculations for these species.
- Allowed for the use of a custom diffusion coefficient even if the species is in the database of the `diffusion_coef` module.

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
