import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import sys
import matplotlib.dates as md

sys.path.append("..")
from funciones import Bpara_Bperp, UTC_to_hdec, datenum, donde, SZA, hdec_to_UTC

"""
Plotea la fig de titan inbound
"""
plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

# plt.rcParams["axes.prop_cycle"] = cycler(
#     "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
# )

path = "../../../datos/Titan/t96_tswis_1s.ascii"
datos = np.loadtxt(path)
tiempo = datos[:, 0]
i = donde(tiempo, 24)
f = donde(tiempo, 26)
t1, t2, t3, t4 = (
    24.7371487623411579,
    24.7811189855693299,
    24.9116972242469316,
    24.9383458443852176,
)
t, B, Bnorm, posicion = datos[i:f, 0], datos[i:f, 1:4], datos[i:f, 4], datos[i:f, 5:8]
B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, t[0] + 0.2, t[-1] - 0.2)

tiempo_mag = np.array([np.datetime64(datenum(2013, 11, 30, x)) for x in t])
tiempo_mag_delta = np.array([np.datetime64(datenum(2013, 11, 30, x)) for x in t_plot])
tm1 = donde(t, t1)
tm2 = donde(t, t2)
tm3 = donde(t, t3)
tm4 = donde(t, t4)
tbs = UTC_to_hdec("25:15:35")
tmpr = 24.58203602
tmva_i = UTC_to_hdec("24:47:42")
tmva_f = UTC_to_hdec("24:52:41")

for t_utc in [
    "24:20",
    "24:25",
    "24:30",
    "24:35",
    "24:40",
    "24:45",
    "24:50",
    "24:55",
    "25:00",
    "25:05",
    "25:10",
    "25:15",
    "25:20",
]:
    idx = donde(t, UTC_to_hdec(t_utc))
    print(t_utc, SZA(posicion, idx) - 90, np.linalg.norm(posicion[idx]) / 2575)

for t1234 in [t1, t2, t3, t4, tmva_i]:
    idx = donde(t, t1234)
    print(
        hdec_to_UTC(t1234),
        SZA(posicion, idx) - 90,
        np.linalg.norm(posicion[idx]) - 2575,
    )


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
ax1 = plt.subplot2grid((3, 1), (0, 0))
ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)

for ax in [ax1, ax2, ax3]:
    regiones(ax, tiempo_mag[donde(t, t1)], tiempo_mag[donde(t, t4)], "#79B953")
    regiones(
        ax,
        tiempo_mag[donde(t, t4)],
        tiempo_mag[donde(t, tbs)],
        "#FE6779",
    )
    regiones(
        ax,
        tiempo_mag[donde(t, tmpr)],
        tiempo_mag[donde(t, t1)],
        "#428AE0",
    )
    lineas_t1t2t3t4(ax, tiempo_mag, t)
    ax.grid()
    ax.set_xlim(
        tiempo_mag[donde(t, 24 + UTC_to_hdec("00:30:00"))],
        tiempo_mag[donde(t, 24 + UTC_to_hdec("01:20:00"))],
    )

ax1.plot(tiempo_mag, Bnorm, linewidth=1)
ax1.set_ylabel("|B| (nT)")
ax1.legend(["MPB", "MS", "MPR"])
ax1.set_title("Cassini MAG 2013-12-01")

ax2.plot(tiempo_mag, B[:, 0], label="Bx TSWIS", linewidth=1)
ax2.plot(tiempo_mag, B[:, 1], label="By TSWIS", linewidth=1)
ax2.plot(tiempo_mag, B[:, 2], label="Bz TSWIS", linewidth=1)
ax2.set_ylabel("Componentes de B (nT)")
ax2.legend(loc="upper right")

ax3.plot(tiempo_mag_delta, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_mag_delta, B_perp_norm, linewidth=1, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Variación relativa de B")
ax3.set_xlabel("Tiempo (hdec)")
if max(B_para) > 1 or max(B_perp_norm) > 1:
    ax3.set_ylim([-0.1, 1])
ax3.legend(loc="upper right")

ax3.set_xlabel(r"Tiempo (UTC) " "\n" r" SZA ($^\circ$) " "\n" r" Distancia (R$_T$)")
ax3.xaxis.set_label_coords(-0.1, -0.03)
ax3.set_xticklabels(
    [
        "00:30\n48\n2.14",
        "00:35\n49\n1.75",
        "00:40\n43\n1.55",
        "00:45\n28\n1.61",
        "00:50\n14\n1.91",
        "00:55\n4\n2.36",
        "00:00\n-2\n2.88",
        "00:05\n-7\n3.45",
        "00:10\n-10\n4.04",
        "00:15\n-13\n4.64",
        "00:20\n-15\n5.25",
    ],
    fontdict=None,
    minor=False,
)

for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)

figure = plt.gcf()  # get current figure
figure.set_size_inches(8, 9)
# when saving, specify the DPI
# plt.savefig(f"{year}-{month}-{day}.png", dpi=150)
plt.show()

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
        ax,
        tiempo_mag[donde(t, tmva_i)],
        tiempo_mag[donde(t, tmva_f)],
        "#FFC300",
        a=0.75,
    )  # MVA
    regiones(
        ax, tiempo_mag[donde(t, t1 - 0.06)], tiempo_mag[donde(t, t1)], "#C70039"
    )  # down
    regiones(
        ax, tiempo_mag[donde(t, t4)], tiempo_mag[donde(t, t4 + 0.06)], "#581845"
    )  # up
    lineas_t1t2t3t4(ax, tiempo_mag, t)
    ax.grid()

axz1.plot(tiempo_mag, Bnorm)
axz1.set_ylabel("|B| (nT)")
axz1.set_title("Cassini MAG 2013-12-01")
axz1.legend(["MVA", "Downstream", "Upstream"], loc="upper right")
plt.setp(axz1.get_xticklabels(), visible=False)

axz2.plot(tiempo_mag, B[:, 0], label="Bx TSWIS")
axz2.plot(tiempo_mag, B[:, 1], label="By TSWIS")
axz2.plot(tiempo_mag, B[:, 2], label="Bz TSWIS")
axz2.set_ylabel("Componentes de B (nT)")
axz2.set_xlabel("Tiempo (UTC)")
axz2.set_xlim(
    tiempo_mag[donde(t, UTC_to_hdec("24:40:01"))],
    tiempo_mag[donde(t, UTC_to_hdec("24:59:59"))],
)
axz2.legend(loc="upper right")
axz2.set_xticklabels(
    [
        "00:45",
        "00:50",
        "00:55",
        # "00:34",
        # "00:35",
        # "00:36",
    ],
    fontdict=None,
    minor=False,
)
# figure = plt.gcf()  # get current figure
# figure.set_size_inches(8, 9)
# # when saving, specify the DPI
# plt.savefig(f"{year}-{month}-{day}.png", dpi=150)
plt.show()
