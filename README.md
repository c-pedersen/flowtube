# flowtube Package
[![PyPI Version](https://img.shields.io/pypi/v/flowtube.svg)](https://img.shields.io/pypi/v/flowtube.svg?cacheSeconds=3600)
[![License: MIT](https://cdn.prod.website-files.com/5e0f1144930a8bc8aace526c/65dd9eb5aaca434fac4f1c34_License-MIT-blue.svg)](/LICENSE)
[![Tests](https://github.com/c-pedersen/flowtube/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/c-pedersen/flowtube/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/c-pedersen/flowtube/branch/main/graph/badge.svg?ts=20251106)](https://codecov.io/gh/c-pedersen/flowtube)

A Python package for transport and diffusion calculations in cylindrical
flow reactors using the KPS method published in Knopf et al., Anal. 
Chem., 2015. The package is currently designed for coated wall reactor 
(CWR) and boat reactor analysis with plans to support aerosol flow 
reactors in future versions.

Author: Corey Pedersen 

## Features
Flow Calculations: Flow rates, velocities, residence times, and Reynolds 
numbers

Transport Properties: Dynamic viscosity and density

Diffusion Analysis: Binary diffusion coefficients, Péclet numbers, and 
mixing parameters

Uptake Modeling: Diffusion correction factors and uptake coefficient 
calculations

Support for Inserts: Handles coated cylindrical inserts within flow 
tubes

Support for Boats: Handles boat reactors which are placed in side of a 
flow tube and filled with an analyte solution.

## Scientific Applications
This package is designed for atmospheric chemistry and aerosol research,
 particularly for:
- Heterogeneous reaction kinetics studies
- Uptake coefficient measurements
- Flow reactor design and optimization
- Transport phenomena analysis in laminar flow systems

## Installation

```bash
pip install flowtube
```

## Acronyms
`sccm = standard cubic centimeter per minute` <br>
`FT = flow tube` <br>
`FR = flow rate` <br>
`MR = mixing ratio` <br>
`CWR = coated wall reactor` <br>

## Assumptions
KPS method assumptions (Knopf et al., 2015):
1. The interacting gas species is a trace gas in the bulk flow
2. Laminar flow is established in the flow reactor
3. The gas temperature and viscosity are homogeneous
4. The axial diffusion velocity is negligible compared to bulk flow 
    velocity
5. The amount of gas species taken up is small compared to its reservoir
6. There is an absence of gas-phase reactions impacting gas species 
    concentration

## Bibliography (more citations within package files)
Knopf, D.A., Pöschl, U., Shiraiwa, M., 2015. Radial Diffusion and 
Penetration of Gas Molecules and Aerosol Particles through Laminar Flow 
Reactors, Denuders, and Sampling Tubes. Anal. Chem. 87, 3746–3754. 
https://doi.org/10.1021/ac5042395
