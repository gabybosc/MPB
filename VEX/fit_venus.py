import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_MAG_pds

import sys

sys.path.append("..")
from funciones import datenum, donde, fechas, tiempos

np.set_printoptions(precision=4)


year, month, day, doy = fechas()
ti, tf = tiempos()

# t, B, pos = importar_MAG_pds(year, doy, ti, tf)
t, B, pos = importar_MAG_pds(2011, 213, 0, 10)
pos_RV = pos / 6050


def altitude(SZA):
    alt = 0.11 * SZA ** 2 - 0.22 * SZA + 389
    return alt / 6050


sza = np.linspace(0, np.pi, 100)
alt = 1 + altitude(sza * 180 / np.pi)

y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

yy = y_alt[x_alt >= 0]
xx = x_alt[x_alt >= 0]

yz = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)

fig, ax = plt.subplots()
ax.plot(pos_RV[:, 0], yz)
ax.plot(xx, yy, color="#5647b4", linestyle="-.")
ax.axis("equal")
ax.set_xlim(0, 3)
ax.set_ylim(0, 2.5)
circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
ax.add_artist(circle)
ax.set_title("VENUS VSO coordinates", fontsize=16)
ax.set_xlabel(r"$X_{VSO}$ ($R_V$)", fontsize=14)
ax.set_ylabel(r"$(Y²_{VSO} + Z²_{VSO} )^{1/2}$ ($R_V$)", fontsize=14)
plt.show()
