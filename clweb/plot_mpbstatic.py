import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mpl
from importar_datos import (
    importar_mag,
    importar_lpw,
    importar_swea,
    importar_swia,
    importar_t1t2t3t4,
    importar_static,
)
from cycler import cycler
import sys

sys.path.append("..")
from funciones import (
    find_nearest,
    donde,
    datenum,
    fechas,
    tiempos,
    UTC_to_hdec,
    proyecciones,
)

np.set_printoptions(precision=4)

"""
Este script plotea mag, swea, swia y lpw en la región de interés para cualquier
cruce.

"""
year, month, day, doy = 2017, "02", "03", "034"  # fechas()
ti, tf = 19, 20  # tiempos()

path = f"../../../datos/clweb/{year}-{month}-{day}/"

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
t1, t2, t3, t4 = (
    19.4571,
    19.4804,
    19.4959,
    19.5088,
)  # importar_t1t2t3t4(year, doy, int(ti))

t_up = t1 - 0.015
t_down = t4 + 0.015

Bnorm = np.linalg.norm(B, axis=1)


# ######### SWEA

swea, t_swea, JE = importar_swea(year, month, day, ti, tf)


# ######################################################################## SWIA

swia, t_swia, density = importar_swia(year, month, day, ti, tf)

# ########################## STATIC
static, t_static, mass, counts = importar_static(year, month, day, ti, tf)

# ############ tiempos UTC
year = int(year)
month = int(month)
day = int(day)

tiempo_mag = np.array(
    [np.datetime64(datenum(year, month, day, x)) for x in t]
)  # datenum es una función mía
tiempo_swea = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swea])
tiempo_swia = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swia])
tiempo_static = np.array(
    [np.datetime64(datenum(year, month, day, x)) for x in t_static]
)

tm1 = np.where(t == find_nearest(t, t1))
tm2 = np.where(t == find_nearest(t, t2))
tm3 = np.where(t == find_nearest(t, t3))
tm4 = np.where(t == find_nearest(t, t4))
# tm_up = np.where(t == find_nearest(t, t_up))
# tm_down = np.where(t == find_nearest(t, t_down))
# tmva = np.where(t == find_nearest(t, np.mean([t1, t2, t3, t4])))
# tbs = np.where(t == find_nearest(t, UTC_to_hdec("18:02:00")))
# tbs_swifa = np.where(t_swifa == find_nearest(t_swifa, UTC_to_hdec("18:02:00")))[0][0]
# tmpr = np.where(t == find_nearest(t, UTC_to_hdec("18:19:00")))

tiempo_lim = [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]

B1, B2, B3 = proyecciones(B)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)
mpl.rcParams.update({"font.size": 12})


#  #######  La fig de la tesis pero sin LPW y con STATIC
fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()

ax1.xaxis.set_major_formatter(xfmt)
ax1 = plt.subplot2grid((5, 1), (1, 0))
ax1.plot(tiempo_mag, B[:, 0], label="Bx MSO")
ax1.plot(tiempo_mag, B[:, 1], label="By MSO")
ax1.plot(tiempo_mag, B[:, 2], label="Bz MSO")
ax1.legend(loc="center right")
ax1.set_ylabel("B components (nT)")

ax2 = plt.subplot2grid((5, 1), (0, 0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, Bnorm)
plt.ylabel("|B| (nT)")
ax2.set_title(f"MAVEN MAG SWEA SWIA {year}-{month}-{day}")

ax3 = plt.subplot2grid((5, 1), (2, 0), sharex=ax1)
ax3.xaxis.set_major_formatter(xfmt)
ax3.set_ylabel("Cuentas de \n masa")
ax3.scatter(tiempo_static, mass, marker="s", c=np.log(counts), cmap="inferno")
ax3.set_yscale("log")


ax4 = plt.subplot2grid((5, 1), (3, 0), sharex=ax1)
ax4.xaxis.set_major_formatter(xfmt)
plt.semilogy(tiempo_swea, JE)
ax4.legend(["50", "100", "150"])
ax4.legend(loc="center right")
ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")


ax5 = plt.subplot2grid((5, 1), (4, 0), sharex=ax1)
ax5.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_swia, density)
ax5.set_ylabel("Proton \n density (cm⁻³)")
ax5.set_xlabel("Time (UTC) \nSZA (º)\nDistance (RM)")
ax5.xaxis.set_label_coords(-0.05, -0.05)


for ax in [ax1, ax2, ax3, ax4]:
    # regiones(ax, tiempo_mag, tm1, tm4, tm_up, tm_down, tbs, tmpr)
    plt.setp(ax.get_xticklabels(), visible=False)


for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()
    # ax.axvline(x=tiempo_mag[tbs][0], color="c", linewidth=1.5)
    for xc in tiempo_lim:
        ax.axvline(x=xc, color="k", linewidth=1.5)


plt.tight_layout()

plt.show()
