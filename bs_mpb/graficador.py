import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from importar_datos import importar_mag_1s, importar_swea, importar_swia
import sys


sys.path.append("..")
from funciones import Bpara_Bperp, hdec_to_UTC

lista = np.loadtxt("../outputs/new_grupo1.txt")
fig_path = "../../BS_MPB/"

for l in lista:
    year, month, day = l[0].split("-")

    t_bs = hdec_to_UTC(l[1])
    t_mpb = hdec_to_UTC(l[2])

    if t_bs < t_mpb:
        ti = t_bs - 0.2
        tf = t_mpb + 0.2
    else:
        ti = t_mpb - 0.2
        tf = t_bs + 0.2
    if ti < 0:
        ti = 0
    if tf > 24:
        tf = 24

    mag, t, B, pos = importar_mag_1s(year, month, day, ti, tf)
    swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
    swia, t_swia, i_density, i_temp, vel_mso = importar_swia(year, month, day, ti, tf)

    Bnorm = np.linalg.norm(B, axis=1)
    Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)

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

    ax1 = plt.subplot2grid((3, 2), (0, 0))
    ax2 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
    ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
    ax4 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
    ax5 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
    ax6 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
    plt.plot(tpara, Bpara, label=r"|$\Delta B \parallel$| / B")
    plt.plot(tpara, Bperp, "-.", label=r"|$\Delta B \perp$| / B")
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r"|$\Delta B$|/ B")
    ax1.set_xlim([t[0], t[-1]])
    ax1.grid()
    ax1.legend()
    ax1.set_title(f"{year}-{month}-{day}")

    ax2.plot(t, B)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.set_ylabel("Bx, By, Bz (nT)")
    ax2.legend(["Bx", "By", "Bz"])
    ax2.grid()

    plt.plot(t, Bnorm)
    ax3.grid()
    ax3.axvline(x=t_bs, color="c")
    ax3.set_ylabel("|B| (nT)")
    ax3.set_xlabel("Tiempo (hdec)")

    if max(Bnorm) > 70:
        ax2.set_ylim([-50, 50])
        ax3.set_ylim([0, 50])

    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.set_xlabel("Tiempo (hdec)")
    ax4.set_ylabel("proton velocity")
    ax4.plot(t_swia, vel_mso)
    ax4.grid()

    plt.setp(ax5.get_xticklabels(), visible=False)
    ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
    ax5.plot(t_swia, i_density)
    ax5.grid()

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

    figure = plt.gcf()  # get current figure
    figure.set_size_inches(16, 8)
    # when saving, specify the DPI
    plt.savefig(fig_path + f"{year}-{month}-{day}.png", dpi=150)
    # plt.show()
