#!/usr/bin/env python3
"""Calculate the height and radius of a thermal
"""
from matplotlib import pyplot as plt
import numpy as np
import nucleator

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


t_p = 5
r_0 = 100
z_end, r_end = compute_thermal(t_ref, t_p, r_0)

print("thermal height: %f" % z_end)

z = np.linspace(0, z_end)
t_0 = t_ref + dtdz * z
t = t_0 + np.linspace(t_p, 0)

# now that we have our temperature anomaly along the termal, we can start doing
# interesting things with humidity

H = 1 # humidity at ground level *everywhere*
p_s = nucleator.p_sat(t+273.15)
p_h2o = H * nucleator.p_sat(t_ref+273.15)
j_star = nucleator.nucleate(t+273.15, p_h2o)


plt.title("Background and in-thermal temperature")
plt.plot(t_0, z, label="$T_0$")
plt.plot(t, z, label=r"$T$")
plt.xlabel("temperature")
plt.ylabel("height")
plt.tight_layout()
plt.legend()
plt.show()
plt.close()

#plt.plot(p_s, z)
plt.plot(p_h2o/p_s, z)
plt.show()
plt.close()

plt.plot(j_star, z)
plt.plot()
plt.show()
