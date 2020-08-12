#!/usr/bin/env python3
"""Calculate the height and radius of a thermal
"""
from matplotlib import pyplot as plt
import numpy as np

alpha = 3400e-6 # Thermal expansion coefficient for air
g = -9.81
m = 2.54
a = 1.90
t_ref = 15

# We use NASA's basic atmosphere model -- more than enough for our purpose
# https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
z = np.arange(0, 3e3)
t_0 = t_ref - 0.00649 * z
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

tp = 5
r = np.arange(2, 2e3)
z_end, r_end = compute_thermal(t_ref, tp, r)

plt.figure()
plt.title('Impact of thermal radius on height (%.0fC increase)' % tp)
plt.plot(r, z_end, label='z_end')
plt.plot(r, r_end, label=t'r_end')
plt.xlabel('Thermal Ground Radius')
plt.legend()


r = 20
tp = np.arange(0, 15)
z_end, r_end = compute_thermal(t_ref, tp, r)

plt.figure()
plt.title('Impact of thermal radius on height R=%.0fm' % r)
plt.plot(tp, z_end, label='z_end')
plt.plot(tp, r_end, label='r_end')
plt.xlabel('Temperature Anomaly')
plt.legend()


plt.show()
plt.close()
