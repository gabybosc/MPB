import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from importar_datos import (
    importar_mag,
    importar_lpw,
    importar_swea,
    importar_swia,
    importar_static,
)
import sys

sys.path.append("..")
from funciones import find_nearest, Bpara_Bperp, fechas, tiempos, donde
from funciones_plot import onpick1


"""
Le paso los datos de clweb y elijo los tiempos t1t2t3t4
Usa los datos de baja resolución para calcular el B_para y B_perp
"""


np.set_printoptions(precision=4)

year, month, day, doy = 2016, "03", 16, 76
ti, tf = 18, 19  #  tiempos()

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
Bnorm = np.linalg.norm(B, axis=1)


# ###################################################################SWEA

swea, t_swea, energias = importar_swea(year, month, day, ti, tf)
energy = swea[:, 7]
JE_total = swea[:, -1]

inicio_swea = donde(t_swea, ti)  # debería ser 0
fin_swea = donde(t_swea, tf)
# ######################################################################## SWIA

swia, t_swia, density = importar_swia(year, month, day, ti, tf)

# ######################################################################### LPW
lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)

# ####################################################################### STATIC
static, t_static, mass, counts = importar_static(year, month, day, ti, tf)

index = np.array((int(year), int(day)))

happy = False
while not happy:
    val = []
    while len(val) < 4:
        plt.clf()  # clear figure
        fig = plt.figure(
            1, constrained_layout=True
        )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
        fig.subplots_adjust(
            top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
        )
        plt.title("Spacebar when ready to click:")

        ax1 = plt.subplot2grid((3, 2), (0, 0))
        plt.plot(t_plot, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
        plt.plot(
            t_plot, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / B"
        )
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r"|$\Delta B$|/ B")

        ax2 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
        ax2.plot(t, B[:, 0], label="Bx")
        ax2.plot(t, B[:, 1], label="By")
        ax2.plot(t, B[:, 2], label="Bz")
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.set_ylabel("Bx, By, Bz (nT)")

        ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
        plt.plot(t, Bnorm)
        ax3.set_ylabel("|B| (nT)")
        ax3.set_xlabel("Tiempo (hdec)")

        ax4 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
        for energia in energias:
            index = np.where(energy == find_nearest(energy, energia))[
                0
            ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
            JE = JE_total[index]
            plt.semilogy(
                t_swea[inicio_swea:fin_swea],
                JE[inicio_swea:fin_swea],
                label=f"{energia} eV",
            )
        ax4.set_ylabel("diff en flux")

        ax5 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
        ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
        plt.plot(t_swia, density)

        ax6 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
        ax6.set_ylabel("Cuentas de \n masa (cm⁻³)")
        ax6.set_xlabel("Tiempo (hdec)")
        ax6.scatter(t, mass, marker="s", c=np.log(counts), cmap="inferno")
        ax6.set_yscale("log")

        for ax in [ax1, ax2, ax4, ax5]:
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.legend()

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.set_xlim(t[0], t[-1])
            ax.grid()

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(
            fig.canvas, (ax1, ax2, ax3, ax4, ax5, ax6), color="black", lw=1
        )

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select MPB: ")
        val = np.asarray(plt.ginput(4))[:, 0]
        print("Selected values: ", val)
        outs = sorted(val)

    print("Happy? Keyboard click for yes, mouse click for no.")
    happy = plt.waitforbuttonpress()

plt.show(block=False)

# with open('../outputs/t1t2t3t4.txt','a') as file:
#     file.write('\n')
#     file.write(f'{year}\t{doy}\t')
#     for k in outs:
#         file.write('{0:1.7}\t'.format(k))
