import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mpl
from cycler import cycler
from _importar_datos import importar_MAG_pds
from matplotlib.widgets import MultiCursor
import sys

sys.path.append("..")
from funciones_plot import onpick1
from funciones import datenum, Bpara_Bperp, day_to_doy


np.set_printoptions(precision=4)

"""
comparar con los datos de 1Hz a ver si están bien calibrados
puedo o hacer un avg o un downsampling
"""

year = 2008
lista = np.loadtxt(f"../outputs/VEX{year}_menor65.txt", dtype=str)

i = int(input("indice en lista\n"))

l = lista[i]
year, month, day = l[0].split("-")
year, doy = day_to_doy(year, month, day)
hh = int(l[1].split(":")[0])
ti = hh - 2
tf = hh + 2
if ti < 0:
    ti = 0
if tf > 24:
    tf = 24

t, B, pos = importar_MAG_pds(year, doy, ti, tf)
Bnorm = np.linalg.norm(B, axis=1)
pos_RV = pos / 6050

Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)

# ############ tiempos UTC


def tiempos_UTC(yy, mm, dd, t):
    tt = np.array([np.datetime64(datenum(yy, mm, dd, x)) for x in t])
    return tt


mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

cant_puntos = 4
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
        ax1 = plt.subplot2grid((3, 1), (0, 0))
        ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
        ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)

        ax1.plot(t, Bnorm, linewidth=0.5)
        ax1.set_ylabel("|B| (nT)")

        ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=0.5)
        ax2.plot(t, B[:, 1], label="By VSO", linewidth=0.5)
        ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=0.5)
        ax2.set_ylabel("B components (nT)")
        ax2.legend(loc="upper left")

        ax3.plot(tpara, Bpara, linewidth=0.5, label="B ||")
        ax3.plot(tpara, Bperp, linewidth=0.5, label="B perp")
        ax3.set_ylabel("variación de Bpara perp")
        ax3.set_xlabel("Tiempo (hdec)")
        ax3.set_ylim([-0.1, 1])
        ax3.legend(loc="upper left")

        ax1.grid()
        ax2.grid()
        ax3.grid()
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(fig.canvas, (ax1, ax2, ax3), color="black", lw=1)

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select MVA region: ")
        val = np.asarray(plt.ginput(cant_puntos))[:, 0]
        print("Selected values: ", val)

    print("Happy? Keyboard click for yes, mouse click for no.")
    happy = plt.waitforbuttonpress()

outs = sorted(val)

guardar = input("save? y/n\n")

if guardar == "y":
    with open("../outputs/VEX_times.txt", "a") as file:
        file.write(f"{year}\t{doy}\t")
        for k in range(len(outs)):
            file.write(f"{outs[k]}\t")
        file.write("\n")
