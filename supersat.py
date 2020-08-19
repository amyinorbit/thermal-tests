#!/usr/bin/env python3
"""Calculate the height and radius of a thermal
"""
from matplotlib import pyplot as plt
import numpy as np

alpha = 3400e-6 # Thermal expansion coefficient for air
g = -9.81
m = 2.54
a = 1.90
t_ref = 18

# We use NASA's basic atmosphere model -- more than enough for our purpose
# https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html

dtdz = -0.00649
N2 = alpha * g * dtdz
print('buoyancy period = %.1f s' % (2 * np.pi / np.sqrt(N2)))

def compute_thermal(t_0, t_p, r_0):
    """Computes the main thermal parameters (z_end and R_end)
    """
    V_0 = m * r_0 ** 3.0
    g_p = alpha * -g * t_p
    B_0 = g_p * V_0

    z_end = 3.17 * (B_0 / N2) ** (1/4)
    r_end = ((4*a*B_0) / (3 * m**2 * N2))**(1/4)
    return (z_end, r_end)

def remap(x, old_lo, old_hi, lo, hi):
    """maps a value from a range to another
    """
    old_span = old_hi - old_lo
    span = hi - lo
    return lo + (x - old_lo) * span / old_span

# Tetens equation, converted to Pascals
# https://en.m.wikipedia.org/wiki/Vapour_pressure_of_water
def p_sat(t):
    return 610.78 * np.exp((17.27 * t) / (237.3 + t))


t_p = 10
r_0 = 100
z_end, r_end = compute_thermal(t_ref, t_p, r_0)

print("thermal height: %f" % z_end)

z = np.linspace(0, z_end)
t_0 = t_ref + dtdz * z
t = t_0 + np.linspace(t_p, 0)

# now that we have our temperature anomaly along the termal, we can start doing
# interesting things with humidity

H = 1 # humidity at ground level *everywhere*
p_s = p_sat(t)
p_s_thermal = H * p_sat(t_ref)




plt.plot(t_0, z, label="T0")
plt.plot(t, z, label="T")
plt.xlabel("temperature")
plt.ylabel("height")
plt.tight_layout()
plt.legend()
plt.show()
plt.close()

#plt.plot(p_s, z)
plt.plot(p_s_thermal/p_s, z)
plt.show()
