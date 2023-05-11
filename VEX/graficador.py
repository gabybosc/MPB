import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_MAG_pds
import sys
from os.path import exists

sys.path.append("..")
from funciones import day_to_doy, Bpara_Bperp

lista = np.loadtxt("../outputs/orbitas_VEX.txt", dtype=str)
fig_path = "../../../../Pictures/VEX/"

for l in lista:
    if float(l[2]) < 60:
        year, month, day = l[0].split("-")
        year, doy = day_to_doy(year, month, day)

        if not exists(fig_path + f"{l[0]}.png"):  # si no está ya la figura
            hh = int(l[1].split(":")[0])
            ti = hh - 2
            tf = hh + 2
            if ti < 0:
                ti = 0
            if tf > 24:
                tf = 24

            t, B, pos = importar_MAG_pds(year, doy, ti, tf)
            Bnorm = np.linalg.norm(B, axis=1)
            Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)

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

            figure = plt.gcf()  # get current figure
            figure.set_size_inches(8, 9)
            # when saving, specify the DPI
            plt.savefig(fig_path + f"{l[0]}.png", dpi=300)
            # plt.show()
