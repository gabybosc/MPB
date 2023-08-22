import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mpl
from importar_datos import (
    importar_mag,
    importar_lpw,
    importar_swea,
    importar_swicfa,
    importar_swia,
    importar_static,
    importar_t1t2t3t4,
)
from cycler import cycler
import sys

sys.path.append("..")
from funciones import (
    find_nearest,
    datenum,
    donde,
    Bpara_Bperp,
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
year, month, day, doy = fechas()
ti, tf = tiempos()


mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day, int(ti))

t_up = t1 - 0.015
t_down = t4 + 0.015

# zoom_inicial = donde(t, t1-0.05)
# zoom_final = donde(t, t4 + 0.05)

Bnorm = np.linalg.norm(B, axis=1)

path = f"../../../datos/MAG_1s/{year}/"  # a los datos de 1s PDS
mag_low = np.loadtxt(
    path + f"mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts", skiprows=160
)
tlow = mag_low[:, 6]  # el dia decimal
tlow = (tlow - int(doy)) * 24  # para que me de sobre la cantidad de horas

Mlow = np.size(tlow)  # el numero de datos
# el campo
Blow = np.zeros((Mlow, 3))
for i in range(7, 10):
    Blow[:, i - 7] = mag_low[:, i]


B_para, B_perp_norm, t_plot = Bpara_Bperp(Blow, tlow, t[0], t[-1])

# ######### SWEA

swea, t_swea, JE = importar_swea(year, month, day, ti, tf)

# ######################################################################## SWIA

swia, t_swia, density, vel_cut, vel_norm = importar_swia(year, month, day, ti, tf)

# t_swica, t_swifa, density_swica, density_swifa = importar_swicfa(
#     year, month, day, ti, tf
# )

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
# tiempo_swica = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swica])
# tiempo_swifa = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swifa])
tiempo_static = np.array(
    [np.datetime64(datenum(year, month, day, x)) for x in t_static]
)
tiempo_low = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_plot])

tm1 = np.where(t == find_nearest(t, t1))
tm2 = np.where(t == find_nearest(t, t2))
tm3 = np.where(t == find_nearest(t, t3))
tm4 = np.where(t == find_nearest(t, t4))
tm_up = np.where(t == find_nearest(t, t_up))
tm_down = np.where(t == find_nearest(t, t_down))
tmva = np.where(t == find_nearest(t, np.mean([t1, t2, t3, t4])))
# tbs = np.where(t == find_nearest(t, UTC_to_hdec("18:02:00")))
# tbs_swifa = np.where(t_swifa == find_nearest(t_swifa, UTC_to_hdec("18:02:00")))[0][0]
# tmpr = np.where(t == find_nearest(t, UTC_to_hdec("18:19:00")))

tiempo_lim = [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]

B1, B2, B3 = proyecciones(B)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)

# ##### funciones para el plot que se repiten
# def regiones(ax, tiempo_mag, tm1, tm4, tm_up, tm_down, tbs, tmpr):
#     ax.axvspan(
#         xmin=tiempo_mag[tbs][0], xmax=tiempo_mag[tm1][0], facecolor="#FE6779", alpha=0.6
#     )  # Magnetosheath
#     ax.axvspan(
#         xmin=tiempo_mag[tm1][0], xmax=tiempo_mag[tm4][0], facecolor="#79B953", alpha=0.6
#     )  # MPB
#     ax.axvspan(
#         xmin=tiempo_mag[tm4][0],
#         xmax=tiempo_mag[tmpr][0],
#         facecolor="#428AE0",
#         alpha=0.5,
#     )  # MPR


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
ax1.set_ylabel("B components (nT)")

ax2 = plt.subplot2grid((5, 1), (0, 0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, Bnorm)
plt.ylabel("|B| (nT)")
ax2.set_title(f"MAVEN MAG SWEA SWIA {year}-{month}-{day}")

ax3 = plt.subplot2grid((5, 1), (2, 0), sharex=ax1)
ax3.xaxis.set_major_formatter(xfmt)
ax3.plot(tiempo_low, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_low, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Relative variation \n of B")

ax4 = plt.subplot2grid((5, 1), (3, 0), sharex=ax1)
ax4.xaxis.set_major_formatter(xfmt)
plt.semilogy(tiempo_swea, JE)
ax4.legend(["50", "100", "150"])
ax4.legend(loc="center right")
ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")


ax5 = plt.subplot2grid((5, 1), (4, 0), sharex=ax1)
ax5.xaxis.set_major_formatter(xfmt)
ax5.plot(tiempo_swia, density)
# plt.plot(tiempo_swica, density_swica, label=r"$n_{sw}$ SWICA")
# plt.plot(tiempo_swifa[:tbs_swifa], density_swifa[:tbs_swifa], label=r"$n_{sw}$ SWIFA")
ax5.set_ylabel("Proton \n density (cm⁻³)")

ax6 = plt.subplot2grid((5, 1), (2, 0), sharex=ax1)
ax6.xaxis.set_major_formatter(xfmt)
ax6.scatter(tiempo_static, mass, marker="s", c=np.log(counts), cmap="inferno")
ax6.set_yscale("log")
ax6.set_xlabel("Time (UTC) \nSZA (º)\nDistance (RM)")
ax6.set_ylabel("Cuentas de \n masa")
ax6.xaxis.set_label_coords(-0.05, -0.05)

# ax6.set_xticklabels(
#     [
#         "17:55\n19\n1.62",
#         "18:00\n10\n1.51",
#         "18:05\n5\n1.39",
#         "18:10\n15\n1.28",
#         "18:15\n29\n1.18",
#         "18:20\n46\n1.10",
#     ],
#     fontdict=None,
#     minor=False,
# )


for ax in [ax1, ax2, ax3, ax4, ax5]:
    # regiones(ax, tiempo_mag, tm1, tm4, tm_up, tm_down, tbs, tmpr)
    plt.setp(ax.get_xticklabels(), visible=False)

for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()
    ax.legend()
    # ax.axvline(x=tiempo_mag[tbs][0], color="c", linewidth=1.5)
    for xc in tiempo_lim:
        ax.axvline(x=xc, color="k", linewidth=1.5)

for xc in tiempo_lim:
    ax6.axvline(x=xc, color="white", linewidth=1.5)

# ax5.text(tiempo_mag[tbs][0] - 50000000, 0.75, "BS", fontsize=14, color="c", rotation=90)
# ax5.text(tiempo_mag[tm1][0] - 50000000, 2, "t1", fontsize=13, color="k", rotation=90)
# ax5.text(tiempo_mag[tm2][0] - 40000000, 0.55, "t2", fontsize=13, color="k", rotation=90)
# ax5.text(tiempo_mag[tm3][0] - 50000000, 2, "t3", fontsize=13, color="k", rotation=90)
# ax5.text(tiempo_mag[tm4][0] - 50000000, 0.55, "t4", fontsize=13, color="k", rotation=90)

plt.tight_layout()

plt.show()
