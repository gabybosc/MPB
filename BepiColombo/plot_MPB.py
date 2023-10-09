import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from leer_datos import importar_bepi
from funciones_bepi import tiempos_UTC
import sys

sys.path.append("..")
from funciones import Bpara_Bperp

np.set_printoptions(precision=4)

"""
comparar con los datos de 1Hz a ver si están bien calibrados
puedo o hacer un avg o un downsampling
"""

t, B, pos = importar_bepi(13.5, 14.1)
Bnorm = np.linalg.norm(B, axis=1)
pos_RV = pos / 6050


Bpara, Bperp, tpara = Bpara_Bperp(B, t, 13.5, 14.1)

yy = 2021
mm = 8
dd = 10
tiempo_mag = tiempos_UTC(yy, mm, dd, t)
tiempo_paraperp = tiempos_UTC(yy, mm, dd, tpara)

outs = [13.8807, 13.8933, 13.9094, 13.929]
MPB = tiempos_UTC(yy, mm, dd, outs)

fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()
ax1 = plt.subplot2grid((3, 1), (0, 0))
ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
for ax in [ax1, ax2, ax3]:
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()

ax1.plot(tiempo_mag, Bnorm, linewidth=0.5)
ax1.set_ylabel("|B| (nT)")
ax1.set_title(f"Bepi-Colombo MAG 2021-08-10")

ax2.plot(tiempo_mag, B[:, 0], label="Bx VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 1], label="By VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 2], label="Bz VSO", linewidth=0.5)
ax2.set_ylabel("B components (nT)")

ax3.plot(tiempo_paraperp, Bpara, linewidth=0.5, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_paraperp, Bperp, "-.", linewidth=0.5, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Relative variation \n of B")
ax3.set_xlabel("Tiempo (UTC)")
ax3.set_ylim([-0.1, 1])


for ax in [ax2, ax3]:
    ax.legend(loc="upper left")
for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
for ax in [ax1, ax2, ax3]:
    ax.axvspan(xmin=MPB[1], xmax=MPB[2], facecolor="#79B953", alpha=0.5)
    ax.axvspan(xmin=MPB[0], xmax=MPB[1], facecolor="#cdcdcd", alpha=0.7)
    ax.axvspan(xmin=MPB[2], xmax=MPB[3], facecolor="#cdcdcd", alpha=0.7)
    # en un radio de 10 min de la MPB
    ax.set_xlim([MPB[0] - np.timedelta64(10, "m"), MPB[-1] + np.timedelta64(10, "m")])


plt.show()
