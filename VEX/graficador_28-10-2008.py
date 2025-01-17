import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from _importar_datos import importar_MAG
import sys
import matplotlib.dates as md

sys.path.append("..")
from funciones import Bpara_Bperp, UTC_to_hdec, datenum, donde, SZA

"""
Plotea la fig del 28-10-2008 - la que uso en la tesis
"""
plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

# plt.rcParams["axes.prop_cycle"] = cycler(
#     "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
# )

lista = np.loadtxt("../outputs/VEX_times.txt")
year, month, day, doy = (2008, 10, 28, 302)
t1, t2, t3, t4 = (
    8.541556527954604,
    8.544405851015947,
    8.551476393427427,
    8.556111111111111,
)
tbs = UTC_to_hdec("08:26:40")
tmpr = UTC_to_hdec("08:34:25")

ti = t1 - 0.2
tf = t4 + 0.2
if ti < 0:
    ti = 0
if tf > 24:
    tf = 24

t, B, pos, cl, tpos = importar_MAG(year, doy, ti, tf)
if cl:
    Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)  # si son datos de clweb 1s
else:
    # para datos de PDS filtrados y diezmados
    Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)
Bnorm = np.linalg.norm(B, axis=1)

# IMA = np.genfromtxt(
#     "dens_20081028-futaana.txt", dtype=str, skip_header=1
# )  # estos son los de futaana pero no me gustan
# t_dens = [
#     (UTC_to_hdec(IMA[i, 0]) + UTC_to_hdec(IMA[i, 1])) / 2 for i in range(len(IMA))
# ]
# dens = [float(IMA[i, 2]) for i in range(len(IMA))]
IMA = np.loadtxt("densidad_20081028-clweb.txt", skiprows=1)  # estos son los del clweb
t_dens = np.array(IMA[:, 0] + IMA[:, 1] / 60 + IMA[:, 2] / 3600)  # hdec
dens = IMA[:, -1]

tiempo_mag = np.array([np.datetime64(datenum(year, month, day, x)) for x in t])
tiempo_para = np.array([np.datetime64(datenum(year, month, day, x)) for x in tpara])
tiempo_dens = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_dens])


def regiones(ax, ti, tf, c, a=0.5):
    ax.axvspan(
        xmin=ti,
        xmax=tf,
        facecolor=c,
        alpha=a,
    )


def lineas_t1t2t3t4(ax, tiempo_mag, t):
    ax.axvline(x=tiempo_mag[donde(t, t1)], color="k", linewidth=1.5)
    ax.axvline(x=tiempo_mag[donde(t, t2)], color="k", linewidth=1.5)
    ax.axvline(x=tiempo_mag[donde(t, t3)], color="k", linewidth=1.5)
    ax.axvline(x=tiempo_mag[donde(t, t4)], color="k", linewidth=1.5)
    ax.axvline(x=tiempo_mag[donde(t, tbs)], color="#FF1493", linewidth=1.5)


"""
fig 1
"""

plt.clf()
fig = plt.figure(1, constrained_layout=True)
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
xfmt = md.DateFormatter("%H:%M")
ax1 = plt.subplot2grid((4, 1), (0, 0))
ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)

for ax in [ax1, ax2, ax3, ax4]:
    regiones(ax, tiempo_mag[donde(t, t1)], tiempo_mag[donde(t, t4)], "#79B953")
    regiones(ax, tiempo_mag[donde(t, tbs)], tiempo_mag[donde(t, t1)], "#FE6779")
    regiones(ax, tiempo_mag[donde(t, t4)], tiempo_mag[donde(t, tmpr)], "#428AE0")
    lineas_t1t2t3t4(ax, tiempo_mag, t)
    ax.grid()

ax1.plot(tiempo_mag, Bnorm, linewidth=1)
ax1.set_ylabel("|B| (nT)")
ax1.legend(["MPB", "MS", "MPR"])
ax1.set_title(f"VEX MAG IMA {year}-{month}-{day}")

ax2.plot(tiempo_mag, B[:, 0], label="Bx VSO", linewidth=1)
ax2.plot(tiempo_mag, B[:, 1], label="By VSO", linewidth=1)
ax2.plot(tiempo_mag, B[:, 2], label="Bz VSO", linewidth=1)
ax2.set_ylabel("Componentes de B (nT)")
ax2.legend(loc="upper right")

ax3.plot(tiempo_para, Bpara, linewidth=1, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_para, Bperp, linewidth=1, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Variación relativa de B")
ax3.set_xlabel("Tiempo (hdec)")
if max(Bpara) > 1 or max(Bperp) > 1:
    ax3.set_ylim([-0.1, 1])
ax3.legend(loc="upper right")

ax4.scatter(tiempo_dens, dens, c="#003f5c")
ax4.set_ylabel("Densidad de \nprotones del SW" + r"(cm$^{-3}$)")
ax4.set_xlim(tiempo_mag[donde(t, 8.4)], tiempo_mag[donde(t, 8.67)])

ax4.set_xlabel("Tiempo (UTC) \nSZA (º)\nDistancia (RM)")
ax4.xaxis.set_label_coords(-0.05, -0.05)
ax4.set_xticklabels(
    [
        "08:25\n33\n1.37",
        "08:30\n49\n1.19",
        "08:35\n70\n1.07",
        "08:40\n95\n1.04",
    ],
    fontdict=None,
    minor=False,
)

for ax in [ax1, ax2, ax3]:
    plt.setp(ax.get_xticklabels(), visible=False)

figure = plt.gcf()  # get current figure
figure.set_size_inches(8, 9)
# when saving, specify the DPI
# plt.savefig(f"{year}-{month}-{day}.png", dpi=150)
plt.show()

# los SZA
tmva = 8.549442222
for tc in [t1, t2, t3, t4, tmva]:
    print(SZA(pos, donde(t, tc)))

"""
fig 2
"""

fig2 = plt.figure(2, constrained_layout=True)
fig2.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
xfmt = md.DateFormatter("%H:%M")
axz1 = plt.subplot2grid((2, 1), (0, 0))
axz2 = plt.subplot2grid((2, 1), (1, 0), sharex=axz1)

for ax in [axz1, axz2]:
    regiones(
        ax, tiempo_mag[donde(t, t2)], tiempo_mag[donde(t, tmva)], "#FFC300", a=0.75
    )  # MVA
    regiones(
        ax, tiempo_mag[donde(t, t1 - 0.015)], tiempo_mag[donde(t, t1)], "#581845"
    )  # up
    regiones(
        ax, tiempo_mag[donde(t, t4)], tiempo_mag[donde(t, t4 + 0.015)], "#C70039"
    )  # down
    lineas_t1t2t3t4(ax, tiempo_mag, t)
    ax.grid()

axz1.plot(tiempo_mag, Bnorm)
axz1.set_ylabel("|B| (nT)")
axz1.set_title(f"VEX MAG IMA {year}-{month}-{day}")
axz1.legend(["MVA", "Upstream", "Downstream"])
plt.setp(axz1.get_xticklabels(), visible=False)

axz2.plot(tiempo_mag, B[:, 0], label="Bx VSO")
axz2.plot(tiempo_mag, B[:, 1], label="By VSO")
axz2.plot(tiempo_mag, B[:, 2], label="Bz VSO")
axz2.set_ylabel("Componentes de B (nT)")
axz2.set_xlabel("Tiempo (UTC)")
axz2.set_xlim(tiempo_mag[donde(t, 8.5)], tiempo_mag[donde(t, 8.6)])
axz2.legend(loc="upper left")
axz2.set_xticklabels(
    [
        "08:31",
        "08:32",
        "08:33",
        "08:34",
        "08:35",
        "08:36",
    ],
    fontdict=None,
    minor=False,
)
# figure = plt.gcf()  # get current figure
# figure.set_size_inches(8, 9)
# # when saving, specify the DPI
# plt.savefig(f"{year}-{month}-{day}.png", dpi=150)
plt.show()

from funciones import SZA, hdec_to_UTC

for t_utc in ["08:25", "08:30", "08:35", "08:40"]:
    idx = donde(t, UTC_to_hdec(t_utc))
    print(t_utc, SZA(pos, idx), np.linalg.norm(pos[idx]) / 6050)

for t1234 in [t1, t2, t3, t4, tmva]:
    idx = donde(t, t1234)
    print(hdec_to_UTC(t1234), SZA(pos, idx), np.linalg.norm(pos[idx]) - 6050)
