import numpy as np


sigma = 72.75e-3 # Converted into SI
dalton = 1.66053906660e-27
m1 = 18.02 * dalton #m1 = 2.99e-23 * 1e-3# Converted into SI
v1 = 2.99e-23 * 1e-6 # converted into SI
kb = 1.38e-23

# Tetens equation
# https://en.wikipedia.org/wiki/Vapour_pressure_of_water
def p_sat(t):
    t -= 273.15
    return 610.78 * np.exp((17.27 * t) / (237.3 + t))

def cluster(t, p_h2o):
    S = p_h2o / p_sat(t)
    return (32.0 * np.pi / 3.0) \
         * ((v1**2.0) * (sigma ** 3.0))/(((kb * t)**3.0) * (np.log(S)**3.0))

def nucleate(t, p_h2o):
    N = p_h2o / (kb * t)
    S = p_h2o / p_sat(t)
    inside = - (16.0 * np.pi / 3.0) \
             * ((v1**2.0) * (sigma**3.0)) \
             / (((kb*t)**3.0) * (np.log(S)**2.0))
    return np.sqrt(2.0 * sigma / (np.pi * m1)) \
         * (v1 * (N**2.0) / S) \
         * np.exp(inside)
