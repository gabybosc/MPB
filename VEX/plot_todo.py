import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import sys
import matplotlib.dates as md
from _importar_datos import importar_MAG, importar_ELS_clweb, importar_t1t2t3t4

sys.path.append("..")
from funciones import find_nearest, Bpara_Bperp, fechas, datenum, UTC_to_hdec

energias = [35, 50, 75, 100]

year, month, day, doy = fechas()
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day)
ti = t1 - 1
tf = t4 + 1

t, B, pos, cl, tpos = importar_MAG(year, doy, ti, tf)
if cl:
    Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)  # si son datos de clweb 1s
else:
    # para datos de PDS filtrados y diezmados
    Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)
Bnorm = np.linalg.norm(B, axis=1)

IMA = np.genfromtxt(f"dens_{year}{month}{day}-futaana.txt", dtype=str, skip_header=1)
t_dens = [
    (UTC_to_hdec(IMA[i, 0]) + UTC_to_hdec(IMA[i, 1])) / 2 for i in range(len(IMA))
]
dens = [float(IMA[i, 2]) for i in range(len(IMA))]

# IMA = np.loadtxt(
#     f"densidad_{year}{month}{day}-clweb.txt", skiprows=1
# )  # estos son los del clweb
# t_dens = np.array(IMA[:, 0] + IMA[:, 1] / 60 + IMA[:, 2] / 3600)  # hdec
# dens = IMA[:, -2]
# vel = IMA[:, -1]

tiempo_mag = np.array(
    [np.datetime64(datenum(int(year), int(month), int(day), x)) for x in t]
)
tiempo_para = np.array(
    [np.datetime64(datenum(int(year), int(month), int(day), x)) for x in tpara]
)
tiempo_dens = np.array(
    [np.datetime64(datenum(int(year), int(month), int(day), x)) for x in t_dens]
)

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
ax4.set_xlabel("Tiempo (UTC)")

for ax in [ax1, ax2, ax3]:
    plt.setp(ax.get_xticklabels(), visible=False)

figure = plt.gcf()  # get current figure
figure.set_size_inches(8, 9)
# when saving, specify the DPI
# plt.savefig(f"{year}-{month}-{day}.png", dpi=150)
plt.show()
#
# fig = plt.figure(
#     1, constrained_layout=True
# )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
# fig.subplots_adjust(
#     top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
# )
#
# ax1 = plt.subplot2grid((4, 1), (0, 0))
# ax1.plot(t, Bnorm, linewidth=0.5)
# ax1.set_ylabel("|B| (nT)")
# ax1.grid()
#
# ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
# ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=0.5)
# ax2.plot(t, B[:, 1], label="By VSO", linewidth=0.5)
# ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=0.5)
# ax2.set_ylabel("B components (nT)")
# ax2.legend(loc="upper left")
# ax2.grid()
#
# ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
# ax3.plot(tpara, Bpara, linewidth=0.5, label="B ||")
# ax3.plot(tpara, Bperp, linewidth=0.5, label="B perp")
# ax3.set_ylabel("variación de Bpara perp")
# ax3.set_xlabel("Tiempo (hdec)")
# ax3.set_ylim([-0.1, 5])
# ax3.legend(loc="upper left")
# ax3.grid()
#
# ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)
# for energia in energias:
#     index = np.where(energy == find_nearest(energy, energia))[
#         0
#     ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
#     plt.semilogy(
#         t_els[index],
#         JE[index],
#         label=f"{energia} eV",
#         linewidth=0.5,
#     )
# ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")
# ax4.legend(loc="center right")
# ax4.grid()
#
# plt.show()
