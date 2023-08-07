import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
import numpy as np
import sys


sys.path.append("..")
from funciones import hdec_to_UTC
from funciones_plot import onpick1


def plot_encontrar(
    frontera,
    fecha,
    tpara,
    B_para,
    B_perp_norm,
    t,
    B,
    t_bs,
    B_norm,
    t_swia,
    vel_mso,
    posicion,
    i_density,
    t_swea,
    JE_pds,
    energias,
):
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
        plt.title(f"Spacebar when ready to click for {frontera}:")

        ax1 = plt.subplot2grid((3, 2), (0, 0))
        plt.plot(tpara, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
        plt.plot(
            tpara,
            B_perp_norm,
            "-.",
            linewidth=1,
            label=r"|$\Delta B \perp$| / B",
        )
        if max(B_para) > 1.2:
            ax1.set_ylim([-0.1, 1])
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r"|$\Delta B$|/ B")
        ax1.set_xlim([t[0], t[-1]])
        ax1.grid()
        ax1.legend()
        ax1.set_title(f"{fecha}  {frontera}")

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

        # ax4 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
        # plt.setp(ax4.get_xticklabels(), visible=False)
        # ax4.set_xlabel("Tiempo (hdec)")
        # ax4.set_ylabel("proton velocity")
        # ax4.plot(t_swia, vel_mso)
        # ax4.grid()
        ax4 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.set_ylabel("posición")
        ax4.plot(t, posicion)
        ax4.grid()

        ax5 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
        plt.setp(ax5.get_xticklabels(), visible=False)
        ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
        ax5.plot(t_swia, i_density)
        ax5.grid()
        if type(i_density) != int:
            if len(i_density) > 0:
                if max(i_density) > 15:
                    ax5.set_ylim([-1, 15])

        ax6 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
        plt.semilogy(t_swea, JE_pds)
        ax6.legend(energias)
        ax6.grid()
        ax6.set_ylabel("Diff. en. flux")

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.axvline(x=t_bs, c="m", label="bs")
        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(
            fig.canvas, (ax1, ax2, ax3, ax4, ax5, ax6), color="black", lw=1
        )

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print(f"Click to select {frontera} limits: ")
        outs = plt.ginput(3)  # [0][0]
        val = sorted(outs)
        val_UTC = [hdec_to_UTC(val[i][0]) for i in range(3)]
        print(
            f"Selected values for {frontera}: ",
            val_UTC,
        )

        print("Happy? Keyboard click for yes, mouse click for no.\n")
        happy = plt.waitforbuttonpress()

    return val_UTC


def BS_MPB(L=0.96, e=0.9, x0=0.78):
    THETA = np.linspace(0, np.pi * 3 / 4, 100)
    PHI = np.linspace(0, 2 * np.pi, 100)

    r1 = L / (1 + e * np.cos(THETA))

    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)
    yz = np.sqrt(Y1**2 + Z1**2)
    return (X1, yz)


def marte(ax, x_bs, yz_bs, x_mpb, yz_mpb):
    ax.plot()
    ax.plot(x_bs, yz_bs, color="#07aec7", linestyle="-.")
    ax.plot(x_mpb, yz_mpb, color="#FF1493", linestyle="-.")
    ax.axis("equal")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 2.5)
    circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
    ax.add_artist(circle)
    ax.set_xlabel(r"$X_{MSO}$ ($R_M$)")
    ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)")


def orbitas(posicion, idx_mpb):
    x = posicion[:, 0]
    y = posicion[:, 1]
    z = posicion[:, 2]
    proyeccion = np.sqrt(y**2 + z**2)
    fig = plt.gcf()
    ax = fig.gca()
    ax.plot(x, proyeccion)
    ax.scatter(x[idx_mpb], proyeccion[idx_mpb])


def plot_encontrar_con_orbita(
    fecha,
    tpara,
    B_para,
    B_perp_norm,
    t,
    B,
    idx_mpb,
    t_mpb,
    B_norm,
    t_swia,
    posicion,
    i_density,
    t_swea,
    JE_pds,
    energias,
):
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
        plt.title(f"Spacebar when ready to click:")

        ax1 = plt.subplot2grid((3, 2), (0, 0))
        plt.plot(tpara, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
        plt.plot(
            tpara,
            B_perp_norm,
            "-.",
            linewidth=1,
            label=r"|$\Delta B \perp$| / B",
        )
        if max(B_para) > 1.2:
            ax1.set_ylim([-0.1, 1])
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r"|$\Delta B$|/ B")
        ax1.set_xlim([t[0], t[-1]])
        ax1.grid()
        ax1.legend()
        ax1.set_title(f"{fecha}")

        ax2 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
        ax2.plot(t, B)
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.set_ylabel("Bx, By, Bz (nT)")
        ax2.legend(["Bx", "By", "Bz"])
        ax2.grid()

        ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
        plt.plot(t, B_norm)
        ax3.grid()
        ax3.axvline(x=t_mpb, color="c")
        ax3.set_ylabel("|B| (nT)")
        ax3.set_xlabel("Tiempo (hdec)")

        if max(B_norm) > 70:
            ax2.set_ylim([-50, 50])
            ax3.set_ylim([0, 50])

        ax4 = plt.subplot2grid((3, 2), (0, 1))
        x_bs, yz_bs = BS_MPB(2.04, 1.03, 0.64)
        x_mpb, yz_mpb = BS_MPB(0.96, 0.9, 0.78)
        marte(ax4, x_bs, yz_bs, x_mpb, yz_mpb)
        orbitas(posicion / 3390, idx_mpb)

        ax5 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
        plt.setp(ax5.get_xticklabels(), visible=False)
        ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
        ax5.plot(t_swia, i_density)
        ax5.grid()
        if type(i_density) != int:
            if len(i_density) > 0:
                if max(i_density) > 15:
                    ax5.set_ylim([-1, 15])

        ax6 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
        plt.semilogy(t_swea, JE_pds)
        ax6.legend(energias)
        ax6.grid()
        ax6.set_ylabel("Diff. en. flux")

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.axvline(x=t_mpb, c="m", label="bs")
        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(
            fig.canvas, (ax1, ax2, ax3, ax4, ax5, ax6), color="black", lw=1
        )

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print(f"Click to select limits: ")
        outs = plt.ginput(3)  # [0][0]
        val = sorted(outs)
        val_UTC = [hdec_to_UTC(val[i][0]) for i in range(3)]
        print(
            f"Selected values: ",
            val_UTC,
        )

        print("Happy? Keyboard click for yes, mouse click for no.\n")
        happy = plt.waitforbuttonpress()

    # return val_UTC
