import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from importar_datos import (
    importar_mag,
    importar_lpw,
    importar_swea,
    importar_swia,
    importar_t1t2t3t4,
    importar_static,
    importar_swea_heatmap,
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
year, month, day, doy = 2016, "03", 16, "076"
# year, month, day, doy = 2017, "02", "03", "034"  # fechas()
ti, tf = 17.5, 19  # tiempos()
# ti, tf = 19, 20.5

path = f"../../../datos/clweb/{year}-{month}-{day}/"

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
t1, t2, t3, t4 = (
    18.2167,
    18.2204,
    18.235,
    18.2476,
    # 19.4571,
    # 19.4804,
    # 19.4959,
    # 19.5088,
)  # importar_t1t2t3t4(year, doy, int(ti))

t_up = t1 - 0.015
t_down = t4 + 0.015

Bnorm = np.linalg.norm(B, axis=1)

# ######### SWEA

# swea, t_swea, JE = importar_swea(year, month, day, ti, tf)
swea, t_swea, energy, swea_counts = importar_swea_heatmap(year, month, day, ti, tf)
# ######################################################################## SWIA

swia, t_swia, density, vel_swia, vel_norm_swia = importar_swia(year, month, day, ti, tf)

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
tm_up = np.where(t == find_nearest(t, t_up))
tm_down = np.where(t == find_nearest(t, t_down))
tmva = np.where(t == find_nearest(t, np.mean([t1, t2, t3, t4])))
tbs = np.where(t == find_nearest(t, UTC_to_hdec("18:02:00")))
# tbs = np.where(t == find_nearest(t, UTC_to_hdec("19:10:00")))
# tbs_swifa = np.where(t_swifa == find_nearest(t_swifa, UTC_to_hdec("18:02:00")))[0][0]
tmpr = np.where(t == find_nearest(t, UTC_to_hdec("18:19:00")))
# tmpr = np.where(t == find_nearest(t, UTC_to_hdec("20:00:00")))

tiempo_lim = [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]

B1, B2, B3 = proyecciones(B)

# mpl.rcParams["axes.prop_cycle"] = cycler(
#     "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
# )
plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)
mpl.rcParams.update({"font.size": 12})

#  #######  La fig de la tesis pero sin LPW y con STATIC
fig = plt.figure(
    2, figsize=(10, 20)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()

ax1.xaxis.set_major_formatter(xfmt)
ax1 = plt.subplot2grid((4, 1), (0, 0))
ax1.plot(tiempo_mag, B[:, 0], label="Bx MSO")
ax1.plot(tiempo_mag, B[:, 1], label="By MSO")
ax1.plot(tiempo_mag, B[:, 2], label="Bz MSO")
ax1.plot(tiempo_mag, Bnorm, linestyle="--", label="|B|")
ax1.legend(loc="upper right")
ax1.set_ylabel("Campo magnético (nT)")
ax1.set_title(f"MAVEN MAG STATIC SWEA SWIA {year}-{month}-{day}")
divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="3%", pad=0.5)
cax1.axis("off")

ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
ax2.set_ylabel("Espectro de masas\n de iones")
z2 = ax2.scatter(tiempo_static, mass, marker="s", c=np.log(counts), cmap="viridis")
ax2.set_yscale("log")
divider2 = make_axes_locatable(ax2)
cax2 = divider2.append_axes("right", size="3%", pad=0.5)
cbar2 = plt.colorbar(z2, cax=cax2)
cbar2.ax.set_ylabel(r"log$_{10}$ JE", rotation=90)

ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
ax3.xaxis.set_major_formatter(xfmt)
# plt.semilogy(tiempo_swea, JE)
# ax4.legend(["50", "100", "150"])
# ax4.legend(loc="center right")
z3 = ax3.scatter(tiempo_swea, energy, marker="s", c=np.log(swea_counts), cmap="viridis")
ax3.set_yscale("log")
ax3.set_ylabel("Flujo diferencial \n de e- del SW \n (cm⁻² sr⁻¹ s⁻¹)")
divider3 = make_axes_locatable(ax3)
cax3 = divider3.append_axes("right", size="3%", pad=0.5)
cbar3 = plt.colorbar(z3, cax=cax3)
cbar3.ax.set_ylabel(r"log$_{10}$ JE", rotation=90)

ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)
ax4.xaxis.set_major_formatter(xfmt)
ax4.plot(tiempo_swia, density)
ax4.set_ylabel("Densidad de \n protones (cm⁻³)")
divider4 = make_axes_locatable(ax4)
cax4 = divider4.append_axes("right", size="3%", pad=0.5)
cax4.axis("off")
# ax4.xaxis.set_label_coords(-0.05, -0.05)

ax4.set_xlabel("Tiempo (UTC) \nSZA (º)\nDistancia (RM)")
ax4.xaxis.set_label_coords(-0.05, -0.05)
ax4.set_xticklabels(
    [
        "17:40\n40\n1.96",
        "17:50\n27\n1.74",
        "18:00\n10\n1.51",
        "18:10\n15\n1.28",
        "18:20\n46\n1.10",
        "18:30\n85\n1.05",
        "18:40\n123\n1.15",
        "18:50\n152\n1.35",
    ],
    fontdict=None,
    minor=False,
)

for ax in [ax1, ax2, ax3]:
    # regiones(ax, tiempo_mag, tm1, tm4, tm_up, tm_down, tbs, tmpr)
    plt.setp(ax.get_xticklabels(), visible=False)

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()
    ax.axvline(x=tiempo_mag[tbs][0], color="k", linewidth=1.5)
    ax.axvline(x=tiempo_mag[tmpr][0], color="k", linewidth=1.5)
    ax.axvline(x=tiempo_mag[tm2][0], color="k", linewidth=1.5)
    # for xc in tiempo_lim:
    #     ax.axvline(x=xc, color="k", linewidth=1.5)

plt.tight_layout()
plt.show()

# from funciones import SZA
#
# for t_utc in ["17:40", "17:50", "18:20", "18:30", "18:40", "18:50"]:
#     idx = donde(t, UTC_to_hdec(t_utc))
#     print(t_utc, SZA(posicion, idx), np.linalg.norm(posicion[idx]) / 3390)
