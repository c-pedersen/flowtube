### Helper Functions ###

## Unit conversions
# temperatures: C to K
T_in_K = lambda t: t + standard_T_K

# pressure: to Pa
P_in_Pa = lambda value, units: value * P_CF[units]


### Constants ###
# standard temperature in K
standard_T_K = 273.15

# standard pressure in Pa
standard_P = 101325

# universal gas constant (kg m2 s-2 K-1 mol-1)
R = 8.3145

# Boltzmann Constant (kg m2 s-2 K-1)
k = 1.380649e-23

# conversion factors to get to units of Pa
P_CF = {'Torr': 133.322, 'bar': 1e5, 'mbar': 100, 'hPa': 100, 'Pa': 1}