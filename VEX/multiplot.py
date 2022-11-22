import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
import numpy as np
import sys

sys.path.append("..")
from funciones import find_nearest
from funciones_plot import onpick1


def multi_plot(t, tpara, t_els, B, Bnorm, Bpara, Bperp, energy, JE, cant_puntos):
    energias = [35, 50, 75, 100]
    happy = False
    while not happy:
        val = []
        while len(val) < cant_puntos:
            plt.clf()  # clear figure
            fig = plt.figure(
                1, constrained_layout=True
            )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
            fig.subplots_adjust(
                top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
            )
            # xfmt = md.DateFormatter("%H:%M")
            ax1 = plt.gca()

            plt.title("Spacebar when ready to click:")
            ax1 = plt.subplot2grid((2, 2), (0, 0))
            ax1.plot(t, Bnorm, linewidth=0.5)
            ax1.set_ylabel("|B| (nT)")
            ax1.grid()

            ax2 = plt.subplot2grid((2, 2), (1, 0), sharex=ax1)
            ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=0.5)
            ax2.plot(t, B[:, 1], label="By VSO", linewidth=0.5)
            ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=0.5)
            ax2.set_ylabel("B components (nT)")
            ax2.legend(loc="upper left")
            ax2.grid()

            ax3 = plt.subplot2grid((2, 2), (0, 1), sharex=ax1)
            ax3.plot(tpara, Bpara, linewidth=0.5, label="B ||")
            ax3.plot(tpara, Bperp, linewidth=0.5, label="B perp")
            ax3.set_ylabel("variación de Bpara perp")
            ax3.set_xlabel("Tiempo (hdec)")
            ax3.set_ylim([-0.1, 5])
            ax3.legend(loc="upper left")
            ax3.grid()

            ax4 = plt.subplot2grid((2, 2), (1, 1), sharex=ax1)
            for energia in energias:
                index = np.where(energy == find_nearest(energy, energia))[
                    0
                ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
                plt.semilogy(
                    t_els[index], JE[index], label=f"{energia} eV", linewidth=0.5,
                )
            ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")
            ax4.legend(loc="center right")
            ax4.grid()

            fig.canvas.mpl_connect("pick_event", onpick1)
            MultiCursor(fig.canvas, (ax1, ax2, ax3, ax4), color="black", lw=1)

            zoom_ok = False
            print("\nSpacebar when ready to click:\n")
            while not zoom_ok:
                zoom_ok = plt.waitforbuttonpress(-1)
            print("Click to select MVA region: ")
            val = np.asarray(plt.ginput(cant_puntos))[:, 0]
            print("Selected values: ", val)

        print("Happy? Keyboard click for yes, mouse click for no.")
        happy = plt.waitforbuttonpress()
        return val
