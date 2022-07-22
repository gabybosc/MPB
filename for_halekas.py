import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as md
import os as os

sys.path.append("..")

from funciones import (
    angulo,
    ancho_mpb,
    corrientes,
    datenum,
    fechas,
    tiempos,
    donde,
    UTC_to_hdec,
)

ti, tf = 17, 18
year, month, day, doy = 2016, "03", 16, "076"


path = f"../../datos/clweb/{year}-{month}-{day}/18/"

dens = np.loadtxt(path + "densidad_solo_protones.asc")
temp = np.loadtxt(path + "temperatura_solo_protones.asc")

density = dens[:, -1]
temperature = temp[:, -1]

t = dens[:, 3] + dens[:, 4] / 60 + dens[:, 5] / 3600  # hdec

tiempo = np.array([np.datetime64(datenum(year, int(month), day, x)) for x in t])

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")


ax1 = plt.subplot2grid((2, 1), (0, 0))
ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)
for ax in [ax1, ax2]:
    # ax.set_xlim(t[12000], t[-17000])
    ax.xaxis.set_major_formatter(xfmt)
    ax.grid()


ax1.plot(tiempo, density)
# ax1.set_ylim(ymin=0.1, ymax=2e5)
ax1.set_ylabel("H⁺ density (cm⁻³)")
ax1.set_title(f"MAVEN SWIFA 600-2000eV {year}-{month}-{day}")
plt.setp(ax1.get_xticklabels(), visible=False)

ax2.plot(tiempo, temperature)
ax2.set_ylabel("temperatue (eV)")
ax2.set_xlabel("Time (UTC)")

plt.show()
