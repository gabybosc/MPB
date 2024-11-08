import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys

sys.path.append("../")
from funciones import Bpara_Bperp, fechas
from funciones_plot import onpick1
from importar_datos import importar_mag_1s, importar_swea, importar_swia

"""
Le paso un día entero y elijo la órbita y luego elijo los tiempos t1t2t3t4,
que va a guardar en un archivo aparte / en gdocs.
Para eso, grafica |B|, BxByBz, swea, swia y lpw.
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
mag, t, B, posicion = importar_mag_1s(year, month, day, 0.1, 24)

plt.plot(t, np.linalg.norm(B, axis=1))
plt.ylim([0, 70])
plt.show()

ti = float(input("Tiempo inicial hdec\n"))
tf = ti + 2
mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
swia, t_swia, i_density, i_temp, vel_mso = importar_swia(year, month, day, ti, tf)
# lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)

B_norm = np.linalg.norm(B, axis=1)

B_para, B_perp_norm, t_plot = Bpara_Bperp(
    B, t, ti + 0.2, tf - 0.2
)  # estos ti, tf tienen que ser menores que el total de datos

index = np.array((int(year), doy))


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
        ax1.grid()
        ax1.legend()

        ax4 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
        ax4.plot(t, B)
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.set_ylabel("Bx, By, Bz (nT)")
        ax4.legend(["Bx", "By", "Bz"])
        ax4.grid()

        ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
        plt.plot(t, B_norm)
        ax3.grid()
        ax3.set_ylabel("|B| (nT)")
        ax3.set_xlabel("Tiempo (hdec)")

        ax5 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
        if swea != 0:
            ax5.set_ylabel("Energia", picker=True)  # , bbox=dict(facecolor='red'))
            plt.setp(ax5.get_xticklabels(), visible=False)
            im = plt.imshow(
                flux_plot,
                aspect="auto",
                origin="lower",
                extent=(t_swea[0], t_swea[-1], energia[-1], energia[0]),
                cmap="inferno",
                norm=LogNorm(vmin=1e4, vmax=1e9),
            )
            divider = make_axes_locatable(ax5)
            cax = divider.append_axes("top", size="7%", pad="1%")
            cb = plt.colorbar(im, cax=cax, orientation="horizontal")
            cax.xaxis.set_ticks_position("top")

        ax7 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
        plt.setp(ax7.get_xticklabels(), visible=False)
        ax7.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
        plt.plot(t_swia, i_density)
        ax7.grid()

        # ax6 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
        # ax6.set_ylabel("Densidad total \n de e- (cm⁻³)")
        # ax6.set_xlabel("Tiempo (hdec)")
        # plt.semilogy(t_lpw, e_density)
        # ax6.grid()

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(fig.canvas, (ax1, ax3, ax4, ax5, ax7), color="black", lw=1)

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select MPB: ")
        val = np.asarray(plt.ginput(2))[:, 0]
        print("Selected values: ", val)
        outs = np.concatenate((index, val))

    print("Happy? Keyboard click for yes, mouse click for no.")
    happy = plt.waitforbuttonpress()

plt.show(block=False)
print(outs)
