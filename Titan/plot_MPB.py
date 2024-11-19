import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.widgets import MultiCursor
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from cycler import cycler

sys.path.append("..")
from funciones_plot import regiones
from funciones import (
    donde,
    Bpara_Bperp,
    datenum,
    Mij,
    autovectores,
    error,
    find_nearest,
    angulo,
    ancho_mpb,
    corrientes,
    UTC_to_hdec,
)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)

path = "../../../datos/Titan/t96_tswis_1s.ascii"
datos = np.loadtxt(path)
tiempo = datos[:, 0]
i = donde(tiempo, 23)
f = donde(tiempo, 26)
t1, t2, t3, t4 = 24.54175182, 24.55058123, 24.57651763, 24.58203602

t, B, Bnorm, posicion = datos[i:f, 0], datos[i:f, 1:4], datos[i:f, 4], datos[i:f, 5:8]
B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, t[0] + 0.2, t[-1] - 0.2)

tiempo_mag = np.array([np.datetime64(datenum(2013, 11, 30, x)) for x in t])
tiempo_mag_delta = np.array([np.datetime64(datenum(2013, 11, 30, x)) for x in t_plot])
tm1 = donde(t, t1)
tm2 = donde(t, t2)
tm3 = donde(t, t3)
tm4 = donde(t, t4)
t_bs = donde(t, 24 + UTC_to_hdec("00:24:20"))
t_mpr = donde(t, 24 + UTC_to_hdec("00:40:30"))

# tiempo_lim = [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]

fig = plt.figure(
    1, constrained_layout=True
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)

ax1 = plt.subplot2grid((3, 1), (0, 0))
plt.plot(tiempo_mag, Bnorm)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_title("Cassini MAG 2013-11-30")
ax1.set_ylabel(r"|$\mathbf{B}$| (nT)")

ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax2.plot(tiempo_mag, B)
ax2.set_ylabel(r"Componentes de $\mathbf{B}$ (nT)")
ax2.legend(["Bx", "By", "Bz"])
plt.setp(ax2.get_xticklabels(), visible=False)

ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
plt.plot(tiempo_mag_delta, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
plt.plot(
    tiempo_mag_delta, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / B"
)
ax3.set_ylabel(r"Variación relativa de $\mathbf{B}$")
ax3.set_xlabel("Tiempo (hdec)")
ax3.legend()

for ax in [ax1, ax2, ax3]:
    ax.set_xlim(
        tiempo_mag[donde(t, 24 + UTC_to_hdec("00:20:00"))],
        tiempo_mag[donde(t, 24 + UTC_to_hdec("00:50:00"))],
    )
    ax.grid()
    for xc in [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]:
        ax.axvline(x=xc, color="k", linewidth=1.5)

regiones(ax1, tiempo_mag, tm1, tm4, t_bs, t_mpr, "yes")
regiones(ax2, tiempo_mag, tm1, tm4, t_bs, t_mpr)
regiones(
    ax3,
    tiempo_mag_delta,
    donde(t_plot, t1),
    donde(t_plot, t4),
    donde(t_plot, 24 + UTC_to_hdec("00:24:20")),
    donde(t_plot, 24 + UTC_to_hdec("00:40:30")),
)
ax1.legend()
plt.show()
