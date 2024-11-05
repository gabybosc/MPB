import matplotlib.pyplot as plt
import numpy as np
from loader import importar_params, importar_posiciones

path = "../../../../datos/bs_mpb/"
"""cone angle (última columna de los catalogos de Jacob actualizados), beta_protones, Mfms, Pdyn, Ls y bueno el Z"""
date, times, pos_bs, pos_mpb, Rsd = importar_posiciones(path)
beta, cone_angle, Mfms, Ls, pdyn = importar_params(path)


def correlaciones(Rsd, color, vmin, vmax, titulo):
    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    plt.scatter(Rsd[1, :], Rsd[0, :], c=color, vmin=vmin, vmax=vmax, cmap="seismic")
    plt.ylabel(r"$R_{SD}$ BS")
    plt.xlabel(r"$R_{SD}$ MPB")
    cbar = plt.colorbar()
    cbar.ax.set_title(titulo)
    plt.grid()


def varianza(n, bins):
    s = 0
    for i in range(len(n)):
        s += n[i] * ((bins[i] + bins[i + 1]) / 2)
    mean = s / np.sum(n)

    t = 0
    for i in range(len(n)):
        t += n[i] * (bins[i] - mean) ** 2
    std = np.sqrt(t / np.sum(n))

    return mean, std


def histograma(hist_mpb1, hist_mpb2, hist_bs1, hist_bs2, label1, label2, dens=True):
    fig = plt.figure()
    fig.subplots_adjust(
        top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
    )
    ax1 = plt.subplot2grid((2, 1), (0, 0))
    ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)

    n1, bins1, patches1 = ax1.hist(hist_mpb1, 20, density=dens, alpha=0.5)
    n2, bins2, patches2 = ax1.hist(hist_mpb2, 20, density=dens, alpha=0.5)
    n3, bins3, patches3 = ax2.hist(hist_bs1, 20, density=dens, alpha=0.5)
    n4, bins4, patches4 = ax2.hist(hist_bs2, 20, density=dens, alpha=0.5)
    ax1.set_xlabel("$R_{SD}^{MPB}$")
    ax2.set_xlabel(r"$R_{SD}^{BS}$")
    if not dens:
        ax1.set_ylabel("Cantidad de eventos")
        ax2.set_ylabel("Cantidad de eventos")
    else:
        ax1.set_ylabel("probability density")
        ax2.set_ylabel("probability density")
    mean1, std1 = varianza(n1, bins1)
    mean2, std2 = varianza(n2, bins2)
    mean3, std3 = varianza(n3, bins3)
    mean4, std4 = varianza(n4, bins4)
    ax1.legend(
        [
            (label1, f"{mean1:.3g}", f"{std1:.3g}"),
            (label2, f"{mean2:.3g}", f"{std2:.3g}"),
        ],
        loc="lower right",
    )
    ax2.legend(
        [
            (label1, f"{mean3:.3g}", f"{std3:.3g}"),
            (label2, f"{mean4:.3g}", f"{std4:.3g}"),
        ],
        loc="lower right",
    )


idx_betamin = [i for i in range(len(beta)) if beta[i] < 1]
idx_betamax = [i for i in range(len(beta)) if beta[i] > 1]

Rsd_betamin = Rsd[:, idx_betamin]
Rsd_betamax = Rsd[:, idx_betamax]

histograma(
    Rsd_betamax[1, :],
    Rsd_betamin[1, :],
    Rsd_betamax[0, :],
    Rsd_betamin[0, :],
    "beta > 1",
    "beta < 1",
    # r"$\beta_p$ > 1",
    # r"$\beta_p$ < 1",
)
correlaciones(Rsd, beta, 0, 6, r"$\beta_p$")
plt.show()

idx_conemin = [i for i in range(len(cone_angle)) if 45 < cone_angle[i] < 135]
idx_conemax = [
    i for i in range(len(cone_angle)) if cone_angle[i] < 45 or cone_angle[i] > 135
]

Rsd_conemin = Rsd[:, idx_conemin]
Rsd_conemax = Rsd[:, idx_conemax]

histograma(
    Rsd_conemax[1, :],
    Rsd_conemin[1, :],
    Rsd_conemax[0, :],
    Rsd_conemin[0, :],
    "cone angle [0-45][135-180]",
    "cone angle [45-135]",
)
correlaciones(Rsd, cone_angle, 0, 180, r"cone angle")
plt.show()

idx_machmin = [i for i in range(len(Mfms)) if Mfms[i] < 5]
idx_machmax = [i for i in range(len(Mfms)) if Mfms[i] > 5]

Rsd_machmin = Rsd[:, idx_machmin]
Rsd_machmax = Rsd[:, idx_machmax]

histograma(
    Rsd_machmax[1, :],
    Rsd_machmin[1, :],
    Rsd_machmax[0, :],
    Rsd_machmin[0, :],
    r"Mfms > 5",
    r"Mfms < 5",
)
correlaciones(Rsd, Mfms, 1, 9, r"Mfms")
plt.show()

idx_pdynmin = [i for i in range(len(pdyn)) if pdyn[i] < 0.5]
idx_pdynmax = [i for i in range(len(pdyn)) if pdyn[i] > 0.5]

Rsd_pdynmin = Rsd[:, idx_pdynmin]
Rsd_pdynmax = Rsd[:, idx_pdynmax]

histograma(
    Rsd_pdynmax[1, :],
    Rsd_pdynmin[1, :],
    Rsd_pdynmax[0, :],
    Rsd_pdynmin[0, :],
    r"$P_{dyn}$ > 0.5",
    r"$P_{dyn}$ < 0.5",
)
correlaciones(Rsd, pdyn, 0, 1, r"$P_{dyn}$")
plt.show()

idx_min = [i for i in range(len(Ls)) if Ls[i] < 180]
idx_max = [i for i in range(len(Ls)) if Ls[i] > 180]

fig, ax = plt.subplots()
ax.set_aspect("equal")
plt.scatter(Rsd[1, idx_min], Rsd[0, idx_min], label="low wave season")
plt.scatter(Rsd[1, idx_max], Rsd[0, idx_max], label="high wave season")
plt.ylabel(r"$R_{SD}^{BS}$")
plt.xlabel(r"$R_{SD}^{MPB}$")
plt.legend()
plt.grid()
plt.show()

Rsd_low = Rsd[:, idx_min]
Rsd_h = Rsd[:, idx_max]

histograma(Rsd_h[1, :], Rsd_low[1, :], Rsd_h[0, :], Rsd_low[0, :], "high", "low")
plt.show()

# bs_sur = [i for i in range(len(pos_bs)) if pos_bs[i, 2] < 0]  # noqa
# bs_norte = [i for i in range(len(pos_bs)) if pos_bs[i, 2] > 0]
# mpb_sur = [i for i in range(len(pos_mpb)) if pos_mpb[i, 2] < 0]
# mpb_norte = [i for i in range(len(pos_mpb)) if pos_mpb[i, 2] > 0]
#
# SS = list(set(bs_sur).intersection(mpb_sur))
# NN = list(set(bs_norte).intersection(mpb_norte))
# SN = list(set(bs_sur).intersection(mpb_norte))
# NS = list(set(bs_norte).intersection(mpb_sur))
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, SS], Rsd[0, SS], label="bs S mpb S")
# plt.scatter(Rsd[1, NN], Rsd[0, NN], label="bs N mpb N")
# plt.scatter(Rsd[1, SN], Rsd[0, SN], label="bs S mpb N")
# plt.scatter(Rsd[1, NS], Rsd[0, NS], label="bs N mpb S")
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.legend()
# plt.grid()
# plt.show()
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, mpb_norte], Rsd[0, mpb_norte], label="mpb N")
# plt.scatter(Rsd[1, mpb_sur], Rsd[0, mpb_sur], label="mpb S")
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.legend()
# plt.grid()
# plt.show()
