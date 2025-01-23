import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.widgets import MultiCursor
import matplotlib as mpl
from cycler import cycler
import sys
from _importar_datos import (
    importar_MAG,
    importar_ELS_clweb,
    importar_t1t2t3t4,
    importar_fila,
)
from _update_parametros import update_varios
from _old_fit_venus import plot_orbita

sys.path.append("..")
from funciones import datenum, donde, Bpara_Bperp, fechas, tiempos, find_nearest
from funciones_plot import onpick1

np.set_printoptions(precision=4)

"""
comparar con los datos de 1Hz a ver si están bien calibrados
puedo o hacer un avg o un downsampling
"""

year, month, day, doy = fechas()

t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day)
# plt.figure()
# plt.plot(t, B)
# plt.show()

ti, tf = t1 - 0 - 5, t4 + 0.5
t, B, pos, cl, tpos = importar_MAG(year, doy, ti, tf)
Bnorm = np.linalg.norm(B, axis=1)


def tiempos_UTC(yy, mm, dd, t):
    tt = np.array([np.datetime64(datenum(yy, mm, dd, x)) for x in t])
    return tt


def altitude(SZA):
    alt = 0.11 * SZA ** 2 - 0.22 * SZA + 389
    return alt / 6050


def fit():
    sza = np.linspace(0, np.pi, 100)
    alt = 1 + altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]

    return xx, yz


mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

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
        # xfmt = md.DateFormatter("%H:%M")

        plt.title("Spacebar when ready to click:")
        ax1 = plt.subplot2grid((2, 1), (0, 0))
        ax1.plot(t, Bnorm, linewidth=1)
        ax1.set_ylabel("|B| (nT)")
        ax1.set_title(f"VEX MAG {year}-{month}-{day}")
        ax1.grid()

        ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)
        ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=1)
        ax2.plot(t, B[:, 1], label="By VSO", linewidth=1)
        ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=1)
        ax2.set_ylabel("B components (nT)")
        ax2.legend(loc="upper left")
        ax2.grid()
        for xc in [t1, t2, t3, t4]:
            for ax in ax1, ax2:
                ax.axvline(x=xc, color="k", linewidth=1)

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(fig.canvas, [ax1, ax2], color="black", lw=1)

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select Bup and Bdown: ")
        val = np.asarray(plt.ginput(4))[:, 0]
        print("Selected values: ", val)
        outs = sorted(val)

    print("Happy? Keyboard click for yes, mouse click for no.")
    happy = plt.waitforbuttonpress()

plt.show()

Bup = np.mean(B[donde(t, outs[0]): donde(t, outs[1]), :], axis=0)
Bdown = np.mean(B[donde(t, outs[2]): donde(t, outs[3]), :], axis=0)

guardar = input("save? y/n\n")

if guardar == "y":
    with open("../outputs/VEX_times.txt", "a") as file:
        file.write(f"{year}\t{doy}\t")
        for k in range(len(outs)):
            file.write(f"{outs[k]}\t")
        file.write("\n")
    nr, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
        year, month, day
    )
    hoja = hoja_parametros

    cell_Bup = hoja.range(f"T{nr}:V{nr}")
    cell_Bdown = hoja.range(f"W{nr}:Y{nr}")

    update_varios(hoja, cell_Bup, Bup)
    update_varios(hoja, cell_Bdown, Bdown)
