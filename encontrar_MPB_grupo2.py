import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from matplotlib.colors import LogNorm
import matplotlib.dates as md
from funciones import Bpara_Bperp, fechas, hdec_to_UTC, UTC_to_hdec
from funciones_plot import onpick1
from importar_datos import importar_mag_1s, importar_lpw, importar_swea, importar_swia
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt


"""
Le paso un día entero y elijo la órbita y luego marco la posición media de la MPB,
Es en hdec porque si no el valor que devuelve es raro ya que lo saca del gráfico y está en xfmt
Lo escribe en una nueva columna en el archivo que hice con catalogo_maven_parte2
Para eso, grafica |B|, BxByBz, swea, swia y lpw.
"""

np.set_printoptions(precision=4)

catalogo = np.genfromtxt("outputs/hoja_grupo2.txt", dtype="str")

i = 65

for k in range(len(catalogo)):
    i = int(input(f"linea (voy por {i})\n"))

    cat = catalogo[i]
    year, month, day = cat[0].split("-")
    t_bs = UTC_to_hdec(cat[1])

    ti = t_bs - 1.5  # mira +- 1.5h respecto del BS
    if ti < 0:
        ti = 0
    tf = t_bs + 1.5
    if tf > 24:
        tf = 24

    mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
    swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
    swia, t_swia, i_density, i_temp, vel_mso = importar_swia(year, month, day, ti, tf)
    lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)

    B_norm = np.linalg.norm(B, axis=1)

    B_para, B_perp_norm, tpara = Bpara_Bperp(
        B, t, ti + 0.2, tf - 0.2
    )  # estos ti, tf tienen que ser menores que el total de datos

    happy = False
    val = False
    while not happy:
        # while val is False:
        plt.clf()  # clear figure
        fig = plt.figure(
            1, constrained_layout=True
        )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
        fig.subplots_adjust(
            top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
        )
        plt.title("Spacebar when ready to click:")

        ax1 = plt.subplot2grid((3, 2), (0, 0))
        plt.plot(tpara, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
        plt.plot(tpara, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / B")
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r"|$\Delta B$|/ B")
        ax1.set_xlim([t[0], t[-1]])
        ax1.grid()
        ax1.legend()

        ax2 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
        ax2.plot(t, B)
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.set_ylabel("Bx, By, Bz (nT)")
        ax2.legend(["Bx", "By", "Bz"])
        ax2.grid()

        ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
        plt.plot(t, B_norm)
        ax3.grid()
        ax3.set_ylabel("|B| (nT)")
        ax3.set_xlabel("Tiempo (hdec)")

        ax4 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
        # ax6.set_ylabel("Densidad total \n de e- (cm⁻³)")
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.set_xlabel("Tiempo (hdec)")
        ax4.semilogy(t_lpw, e_density)
        ax4.grid()

        ax5 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
        plt.setp(ax5.get_xticklabels(), visible=False)
        ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
        ax5.plot(t_swia, i_density)
        ax5.grid()

        ax6 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
        if swea != 0:
            ax6.set_ylabel("Energia", picker=True)  # , bbox=dict(facecolor='red'))
            plt.setp(ax6.get_xticklabels(), visible=False)
            im = plt.imshow(
                flux_plot,
                aspect="auto",
                origin="lower",
                extent=(t_swea[0], t_swea[-1], energia[-1], energia[0]),
                cmap="inferno",
                norm=LogNorm(vmin=1e4, vmax=1e9),
            )
            divider = make_axes_locatable(ax6)
            cax = divider.append_axes("top", size="7%", pad="1%")
            cb = plt.colorbar(im, cax=cax, orientation="horizontal")
            cax.xaxis.set_ticks_position("top")

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(
            fig.canvas, (ax1, ax2, ax3, ax4, ax5, ax6), color="black", lw=1
        )

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select MPB: ")
        val = plt.ginput(1)[0][0]
        print("Selected values: ", hdec_to_UTC(val))

        print("Happy? Keyboard click for yes, mouse click for no.\n")
        happy = plt.waitforbuttonpress()

    with open("outputs/grupo2.txt", "a") as file:
        file.write(f"{cat[0]}\t{cat[1]}\t{hdec_to_UTC(val)}\t{cat[2]}\t{cat[3]}")
        file.write("\n")
