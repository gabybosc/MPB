import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_MAG
import sys


sys.path.append("..")
from funciones import Bpara_Bperp, doy_to_day

"""
Plotea y guarda automáticamente las figs con datos de MAG y los tiempos t1t2t3t4
"""


lista = np.loadtxt("../outputs/VEX_times.txt")
fig_path = "../outputs/VEX/MPB/"

for l in lista:
    year = int(l[0])
    doy = str(int(l[1])).zfill(3)

    year, month, day = doy_to_day(year, doy)

    ti = l[2] - 0.2
    tf = l[4] + 0.2
    if ti < 0:
        ti = 0
    if tf > 24:
        tf = 24

    t, B, pos, cl = importar_MAG(year, doy, ti, tf)
    if cl == True:
        Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)  # si son datos de clweb 1s
    else:
        # para datos de PDS filtrados y diezmados
        Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)
    Bnorm = np.linalg.norm(B, axis=1)

    plt.clf()
    fig = plt.figure(1, constrained_layout=True)
    fig.subplots_adjust(
        top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
    )

    ax1 = plt.gca()

    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax1.plot(t, Bnorm, linewidth=0.5)
    ax1.set_ylabel("|B| (nT)")
    ax1.set_title(f"VEX MAG {year}-{month}-{day}")
    ax1.grid()

    ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
    ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=0.5)
    ax2.plot(t, B[:, 1], label="By VSO", linewidth=0.5)
    ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=0.5)
    ax2.set_ylabel("B components (nT)")
    ax2.legend(loc="upper left")
    ax2.grid()

    ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
    ax3.plot(tpara, Bpara, linewidth=0.5, label="B ||")
    ax3.plot(tpara, Bperp, linewidth=0.5, label="B perp")
    ax3.set_ylabel("variación de Bpara perp")
    ax3.set_xlabel("Tiempo (hdec)")
    ax3.set_ylim([-0.1, 5])
    ax3.legend(loc="upper left")
    ax3.grid()

    for ax in [ax1, ax2]:
        plt.setp(ax.get_xticklabels(), visible=False)
    for ax in [ax1, ax2, ax3]:
        ax.axvspan(xmin=l[2], xmax=l[-1], facecolor="#79B953", alpha=0.5)
        # ax.axvspan(xmin=MPB[0], xmax=MPB[1], facecolor="#cdcdcd", alpha=0.7)
        # ax.axvspan(xmin=MPB[2], xmax=MPB[3], facecolor="#cdcdcd", alpha=0.7)
    figure = plt.gcf()  # get current figure
    figure.set_size_inches(8, 9)
    # when saving, specify the DPI
    plt.savefig(fig_path + f"{year}-{month}-{day}.png", dpi=150)
    # plt.show()
