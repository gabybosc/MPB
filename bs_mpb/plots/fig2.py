import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import sys
from os.path import exists
from loader import importar_tiempos
from importar_datos import importar_mag_1s, importar_swea, importar_swia

sys.path.append("../..")
from funciones import Bpara_Bperp, UTC_to_hdec, donde

plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

path = "../../../../datos/bs_mpb/"

date, t_BS, t_MPB = importar_tiempos(path)

g = 50

year, month, day = date[g].split("-")
t_mpb = []
t_bs = []
for i in range(3):
    t_mpb.append(UTC_to_hdec(t_MPB[i][g]))
    t_bs.append(UTC_to_hdec(t_BS[i][g]))

if t_mpb[0] < t_bs[0]:
    ti = t_mpb[0] - 0.25
    tf = t_bs[2] + 0.25

if t_mpb[0] > t_bs[0]:
    ti = t_bs[0] - 0.25
    tf = t_mpb[2] + 0.25

mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
swia, t_swia, i_density, i_temp, vel_mso = importar_swia(year, month, day, ti, tf)
energias = [50 + i * 25 for i in range(6)]
if type(t_swea) is not int:
    JE_pds = np.zeros((len(t_swea), len(energias)))

    for i, e in enumerate(energias):
        j = donde(energia, e)
        JE_pds[:, i] = flux_plot[j]
else:
    JE_pds = 0

# idx_mpb = donde(t, t_mpb[1])
Bnorm = np.linalg.norm(B, axis=1)
Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)

plt.clf()
fig = plt.figure(1, constrained_layout=True)
fig.subplots_adjust(
    top=0.95,
    bottom=0.1,
    left=0.05,
    right=0.95,
    hspace=0.005,
    wspace=0.15,
)
plt.title("Spacebar when ready to click:")

ax1 = plt.subplot2grid((5, 1), (0, 0))
ax2 = plt.subplot2grid((5, 1), (1, 0), sharex=ax1)
ax3 = plt.subplot2grid((5, 1), (2, 0), sharex=ax1)
ax4 = plt.subplot2grid((5, 1), (3, 0), sharex=ax1)
ax5 = plt.subplot2grid((5, 1), (4, 0), sharex=ax1)

ax1.plot(tpara, Bpara, label=r"|$\Delta B \parallel$| / B", linewidth=0.5)
ax1.plot(tpara, Bperp, "-.", label=r"|$\Delta B \perp$| / B", linewidth=0.5)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel(r"|$\Delta B$|/ B")
ax1.set_xlim([t[0], t[-1]])
if max(Bpara) > 1:
    ax1.set_ylim([-0.1, 1])
ax1.grid()
ax1.legend()
ax1.set_title(f"{year}-{month}-{day}")

ax2.plot(t, B)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("Bx, By, Bz (nT)")
ax2.legend(["Bx", "By", "Bz"])
ax2.grid()

ax3.plot(t, Bnorm)
ax3.grid()
plt.setp(ax3.get_xticklabels(), visible=False)
if max(Bnorm) > 70 and Bnorm[donde(t, t_mpb[1])] < 40:
    ax2.set_ylim([-50, 50])
    ax3.set_ylim([0, 50])
if Bnorm[donde(t, t_mpb[1])] < 20:
    ax2.set_ylim([-20, 20])
    ax3.set_ylim([0, 30])
elif max(Bnorm) > 70 and Bnorm[donde(t, t_mpb[1])] > 40:
    ax2.set_ylim([-100, 100])
    ax3.set_ylim([0, 100])
ax3.set_ylabel("|B| (nT)")

plt.setp(ax4.get_xticklabels(), visible=False)
ax4.semilogy(t_swea, JE_pds)
ax4.legend(energias, loc="upper right")
ax4.grid()
ax4.set_ylabel("Diff. en. flux")

ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
ax5.plot(t_swia, i_density)
if type(i_density) is not int:
    if max(i_density) > 30 and i_density[donde(t_swia, t_mpb[1])] < 20:
        ax5.set_ylim([-0.1, 20])
ax5.grid()
ax5.set_xlabel("Tiempo (hdec)")

for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.axvline(x=t_mpb[1], c="#FF1493", label="MPB")
    ax.axvline(x=t_bs[1], c="#07aec7", label="BS")
    ax.axvspan(xmin=t_mpb[0], xmax=t_mpb[2], facecolor="#FF1493", alpha=0.3)
    ax.axvspan(xmin=t_bs[0], xmax=t_bs[2], facecolor="#07aec7", alpha=0.3)
ax5.legend()

figure = plt.gcf()  # get current figure
figure.set_size_inches(16, 8)
plt.show()
