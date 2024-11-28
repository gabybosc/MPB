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
i = donde(tiempo, 24)
f = donde(tiempo, 25)
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
tmva_i = donde(t, UTC_to_hdec("24:33:05"))
tmva_f = donde(t, UTC_to_hdec("24:33:25"))

fig = plt.figure(
    1, constrained_layout=True
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)

ax1 = plt.subplot2grid((3, 1), (0, 0))
plt.plot(tiempo_mag, Bnorm)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_title("Cassini MAG 2013-12-01")
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
ax3.set_ylim([-0.1, 2.1])
ax3.set_ylabel(r"Variación relativa de $\mathbf{B}$")
ax3.set_xlabel("Tiempo (UTC)")
ax3.legend()

for ax in [ax1, ax2, ax3]:
    ax.set_xlim(
        tiempo_mag[donde(t, 24 + UTC_to_hdec("00:20:00"))],
        tiempo_mag[donde(t, 24 + UTC_to_hdec("00:50:00"))],
    )
    ax.grid()
    for xc in [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]:
        ax.axvline(x=xc, color="k", linewidth=1.5)

regiones(ax1, tiempo_mag, tm1, tm4, t_bs, tm1, tm4, t_mpr, "yes")
regiones(ax2, tiempo_mag, tm1, tm4, t_bs, tm1, tm4, t_mpr)
regiones(
    ax3,
    tiempo_mag_delta,
    donde(t_plot, t1),
    donde(t_plot, t4),
    donde(t_plot, 24 + UTC_to_hdec("00:24:20")),
    donde(t_plot, t1),
    donde(t_plot, t4),
    donde(t_plot, 24 + UTC_to_hdec("00:40:30")),
)
ax1.legend()
plt.show()

fig = plt.figure(
    1, constrained_layout=True
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)

ax1 = plt.subplot2grid((2, 1), (0, 0))
plt.plot(tiempo_mag, Bnorm)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_title("Cassini MAG 2013-12-01")
ax1.set_ylabel(r"|$\mathbf{B}$| (nT)")

ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)
ax2.plot(tiempo_mag, B)
ax2.set_ylabel(r"Componentes de $\mathbf{B}$ (nT)")
ax2.legend(["Bx", "By", "Bz"])
ax2.set_xlabel("Tiempo (UTC)")

for ax in [ax1, ax2]:
    ax.set_xlim(
        tiempo_mag[donde(t, 24 + UTC_to_hdec("00:30:00"))],
        tiempo_mag[donde(t, 24 + UTC_to_hdec("00:37:00"))],
    )
    ax.grid()
    for xc in [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]:
        ax.axvline(x=xc, color="k", linewidth=1.5)

regiones(
    ax1,
    tiempo_mag,
    donde(t, t1 - 0.015),
    tm1,
    tmva_i,
    tmva_f,
    tm4,
    donde(t, t4 + 0.015),
    label_yn="yes",
    label_r1="Upstream",
    label_r2="MVA",
    label_r3="Downstream",
    color_r1="#581845",
    color_r2="#FFC300",
    color_r3="#C70039",
)
regiones(
    ax2,
    tiempo_mag,
    donde(t, t1 - 0.015),
    tm1,
    tmva_i,
    tmva_f,
    tm4,
    donde(t, t4 + 0.015),
    label_yn="no",
    label_r1="Upstream",
    label_r2="MVA",
    label_r3="Downstream",
    color_r1="#581845",
    color_r2="#FFC300",
    color_r3="#C70039",
)
ax1.legend()
plt.show()

pos_RT = posicion / 2575
orbita = np.sqrt(pos_RT[:, 1] ** 2 + pos_RT[:, 2] ** 2)
nmva = np.array([0.81616256, 0.47551288, -0.32827759])
nmva_2d = np.array([0.81616256, np.sqrt(0.47551288**2 + -(0.32827759**2))])
R = np.array(
    [
        pos_RT[donde(t_plot, t1), 0],
        np.sqrt(pos_RT[donde(t_plot, t1), 1] ** 2 + pos_RT[donde(t_plot, t1), 2] ** 2),
    ]
)

# fig, ax = plt.subplots()
# ax.plot(pos_RT[:, 0], orbita)
# ax.scatter(
#     pos_RT[0, 0], orbita[0], s=50, zorder=2, marker="o", color="r", label="start"
# )
# ax.scatter(
#     pos_RT[-1, 0], orbita[-1], s=50, zorder=2, marker="x", color="k", label="end"
# )
# ax.quiver(R[0], R[1], nmva[0], nmva[1], color="C2", scale=10, label="MVA")
# ax.scatter(R[0], R[1])
# ax.axis("equal")
# ax.set_xlim(0, 2.5)
# ax.set_ylim(0, 2.5)
# circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
# ax.add_artist(circle)
# ax.set_title("Titan TSWIS coordinates", fontsize=16)
# ax.set_xlabel(r"$X_{TSWIS}$ ($R_T$)", fontsize=14)
# ax.set_ylabel(r"$(Y²_{TSWIS} + Z²_{TSWIS} )^{1/2}$ ($R_T$)", fontsize=14)
# plt.legend()

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)

ax1 = plt.subplot2grid((3, 1), (0, 0))
ax1.plot(pos_RT[:, 0], pos_RT[:, 1])
ax1.add_artist(circle)
ax1.set_title("Titan TSWIS coordinates", fontsize=16)
ax1.set_xlabel(r"$X_{TSWIS}$ ($R_T$)", fontsize=14)
ax1.set_ylabel(r"$Y_{TSWIS}$ ($R_T$)", fontsize=14)

circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
ax2 = plt.subplot2grid((3, 1), (1, 0))
ax2.plot(pos_RT[:, 0], pos_RT[:, 2])
ax2.add_artist(circle)
ax2.set_title("Titan TSWIS coordinates", fontsize=16)
ax2.set_xlabel(r"$X_{TSWIS}$ ($R_T$)", fontsize=14)
ax2.set_ylabel(r"$Z_{TSWIS}$ ($R_T$)", fontsize=14)

circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
ax3 = plt.subplot2grid((3, 1), (2, 0))
ax3.plot(pos_RT[:, 1], pos_RT[:, 2])
ax3.add_artist(circle)
ax3.set_title("Titan TSWIS coordinates", fontsize=16)
ax3.set_xlabel(r"$Y_{TSWIS}$ ($R_T$)", fontsize=14)
ax3.set_ylabel(r"$Z_{TSWIS}$ ($R_T$)", fontsize=14)
for ax in [ax1, ax2, ax3]:
    ax.set_ylim(-5, 5)
    ax.set_xlim(-5, 5)
    ax.axis("equal")
    # ax.legend()

plt.show()
