import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from matplotlib.colors import LogNorm
import matplotlib.dates as md
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt

sys.path.append("..")
from funciones import Bpara_Bperp, fechas, hdec_to_UTC, UTC_to_hdec
from funciones_plot import onpick1
from importar_datos import importar_mag_1s, importar_lpw, importar_swea, importar_swia


"""
Le paso un día entero y elijo la órbita y luego marco la posición media de la MPB,
Es en hdec porque si no el valor que devuelve es raro ya que lo saca del gráfico y está en xfmt
Lo escribe en una nueva columna en el archivo que hice con catalogo_maven_parte2
Para eso, grafica |B|, BxByBz, swea, swia y lpw.M   

Después de este correr after_encontrar
"""

np.set_printoptions(precision=4)
grupo = 3
catalogo = np.genfromtxt(f"../outputs/grupo{grupo}/primer_corte.txt", dtype="str")
# ek grupo 4 está bastante bien distrubuido en ángulos

i = int(input("numero de lista\n"))  # hicimos hasta !!

# for i in range(len(catalogo)):
cat = catalogo[i]
year, month, day = cat[0].split("-")

# if year == "2018" and month == "11":  # sólo quiero elegir de este mes

t_bs = UTC_to_hdec(cat[1])
t_mpb = t_bs + 0.0001  # UTC_to_hdec(cat[2])
# print(i, year, month, day, t_bs)

if t_bs < t_mpb:
    ti = t_bs - 2
    tf = t_mpb + 2
else:
    ti = t_mpb - 2
    tf = t_bs + 2
if ti < 0:
    ti = 0
if tf > 24:
    tf = 24
mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
swia, t_swia, i_density, i_temp, vel_mso = importar_swia(year, month, day, ti, tf)
# swia, t_swia, i_density, i_temp, vel_mso = 0, 0, 0, 0, 0

B_norm = np.linalg.norm(B, axis=1)

B_para, B_perp_norm, tpara = Bpara_Bperp(
    B, t, ti + 0.2, tf - 0.2
)  # estos ti, tf tienen que ser menores que el total de datos

happy = False
val = False
while not happy:
    # while val is False:
    plt.clf()  # clear figure
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

    ax1 = plt.subplot2grid((3, 2), (0, 0))
    plt.plot(tpara, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
    plt.plot(
        tpara,
        B_perp_norm,
        "-.",
        linewidth=1,
        label=r"|$\Delta B \perp$| / B",
    )
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r"|$\Delta B$|/ B")
    ax1.set_xlim([t[0], t[-1]])
    ax1.grid()
    ax1.legend()
    ax1.set_title(f"{year}-{month}-{day}")

    ax2 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
    ax2.plot(t, B)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.set_ylabel("Bx, By, Bz (nT)")
    ax2.legend(["Bx", "By", "Bz"])
    ax2.grid()

    ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
    plt.plot(t, B_norm)
    ax3.grid()
    ax3.axvline(x=t_bs, color="c")
    ax3.set_ylabel("|B| (nT)")
    ax3.set_xlabel("Tiempo (hdec)")

    if max(B_norm) > 70:
        ax2.set_ylim([-50, 50])
        ax3.set_ylim([0, 50])

    ax4 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.set_xlabel("Tiempo (hdec)")
    ax4.set_ylabel("proton velocity")
    ax4.plot(t_swia, vel_mso)
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

    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.axvline(x=t_bs, c="m", label="bs")
        ax.axvline(x=t_mpb, c="g", label="mpb")
    fig.canvas.mpl_connect("pick_event", onpick1)
    multi = MultiCursor(fig.canvas, (ax1, ax2, ax3, ax4, ax5, ax6), color="black", lw=1)

    zoom_ok = False
    print("\nSpacebar when ready to click:\n")
    while not zoom_ok:
        zoom_ok = plt.waitforbuttonpress(-1)
    print("Click to select MPB: ")
    val = plt.ginput(2)  # [0][0]
    print("Selected values: ", hdec_to_UTC(val[0][0]), hdec_to_UTC(val[1][0]))

    print("Happy? Keyboard click for yes, mouse click for no.\n")
    happy = plt.waitforbuttonpress()

with open(f"../outputs/grupo{grupo}/segundo_corte.txt", "a") as file:
    file.write(
        f"{cat[0]}\t{hdec_to_UTC(val[0][0])}\t{hdec_to_UTC(val[1][0])}\t{cat[3]}\t{cat[4]}"
    )
    file.write("\n")
