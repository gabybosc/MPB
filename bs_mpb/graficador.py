import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from cycler import cycler
from mpl_toolkits.axes_grid1 import make_axes_locatable
from figura_encontrar import BS_MPB, marte, orbitas
import sys
from os.path import exists


sys.path.append("..")
from importar_datos import importar_mag_1s, importar_swea, importar_swia
from funciones import Bpara_Bperp, UTC_to_hdec, donde

plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

# grupo = input("grupo\n")
for grupo in [1, 2, 3, 4]:
    lista = np.genfromtxt(
        f"../outputs/grupo{grupo}/jacob_dayside_bs.txt", skip_header=1, dtype=str
    )
    # fig_path = f"../../../../Pictures/BS_MPB/grupo{grupo}_Jacob/"  # casa
    fig_path = f"../../Pictures/BS_MPB/grupo{grupo}_Jacob/"  # iafe

    for l in lista:
        flag = l[-3]
        year, month, day = l[0].split("-")
        t_bs = UTC_to_hdec(l[1])
        t_mpb = (
            UTC_to_hdec(l[2]),
            UTC_to_hdec(l[3]),
            UTC_to_hdec(l[4]),
        )

        if not exists(
            fig_path + f"{year}-{month}-{day}-{int(t_mpb[1])}.png"
        ):  # si no está ya la figura
            if t_mpb[1] < t_bs:
                ti = t_mpb[1] - 0.25
                tf = t_bs + 0.25

            if t_mpb[1] > t_bs:
                ti = t_bs - 0.25
                tf = t_mpb[1] + 0.25

            if ti < 0:
                ti = 0.2
            if tf > 24:
                tf = 24

            mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
            swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
            swia, t_swia, i_density, i_temp, vel_mso = importar_swia(
                year, month, day, ti, tf
            )
            energias = [50 + i * 25 for i in range(6)]
            if type(t_swea) != int:
                JE_pds = np.zeros((len(t_swea), len(energias)))

                for i, e in enumerate(energias):
                    j = donde(energia, e)
                    JE_pds[:, i] = flux_plot[j]
            else:
                JE_pds = 0

            idx_mpb = donde(t, t_mpb[1])
            Bnorm = np.linalg.norm(B, axis=1)
            Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)

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

            ax1.plot(tpara, Bpara, label=r"|$\Delta B \parallel$| / B")
            ax1.plot(tpara, Bperp, "-.", label=r"|$\Delta B \perp$| / B")
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.set_ylabel(r"|$\Delta B$|/ B")
            ax1.set_xlim([t[0], t[-1]])
            if max(Bpara) > 1:
                ax1.set_ylim([-0.1, 1])
            ax1.grid()
            ax1.legend()
            ax1.set_title(f"{year}-{month}-{day}")

            ax2.plot(t, B)
            plt.setp(ax2.get_xticklabels(), visible=False)
            ax2.set_ylabel("Bx, By, Bz (nT)")
            ax2.legend(["Bx", "By", "Bz"])
            ax2.grid()

            ax3.plot(t, Bnorm)
            ax3.grid()
            if max(Bnorm) > 70 and Bnorm[donde(t, t_mpb[1])] < 40:
                ax2.set_ylim([-50, 50])
                ax3.set_ylim([0, 50])
            if Bnorm[donde(t, t_mpb[1])] < 20:
                ax2.set_ylim([-20, 20])
                ax3.set_ylim([0, 30])
            elif max(Bnorm) > 70 and Bnorm[donde(t, t_mpb[1])] > 40:
                ax2.set_ylim([-100, 100])
                ax3.set_ylim([0, 100])
            ax3.set_ylabel("|B| (nT)")
            ax3.set_xlabel("Tiempo (hdec)")

            ax4 = plt.subplot2grid((3, 2), (0, 1))
            plt.setp(ax4.get_xticklabels(), visible=False)
            x_bs, yz_bs = BS_MPB(2.04, 1.03, 0.64)
            x_mpb, yz_mpb = BS_MPB(0.96, 0.9, 0.78)
            marte(ax4, x_bs, yz_bs, x_mpb, yz_mpb)
            orbitas(posicion / 3390, idx_mpb, t)

            plt.setp(ax5.get_xticklabels(), visible=False)
            ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
            ax5.plot(t_swia, i_density)
            if type(i_density) != int:
                if max(i_density) > 30 and i_density[donde(t_swia, t_mpb[1])] < 20:
                    ax5.set_ylim([-0.1, 20])
            ax5.grid()

            ax6.semilogy(t_swea, JE_pds)
            ax6.legend(energias, loc="upper right")
            ax6.grid()
            ax6.set_ylabel("Diff. en. flux")

            for ax in [ax1, ax2, ax3, ax5, ax6]:
                ax.axvline(x=t_mpb[1], c="#FF1493")
                ax.axvline(x=t_bs, c="#07aec7")
                ax.axvspan(xmin=t_mpb[0], xmax=t_mpb[2], facecolor="#79B953", alpha=0.6)

            figure = plt.gcf()  # get current figure
            figure.set_size_inches(16, 8)
            # when saving, specify the DPI
            if flag == "0":
                plt.savefig(
                    fig_path + f"flagged-{year}-{month}-{day}-{int(t_mpb[1])}.png",
                    dpi=150,
                )
            else:
                plt.savefig(
                    fig_path + f"{year}-{month}-{day}-{int(t_mpb[1])}.png", dpi=150
                )
            # plt.show()
