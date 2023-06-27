import numpy as np
import matplotlib.pyplot as plt
import func_position as fpos
from os.path import exists
import Vignesfit_functions as fvig
from generar_npys import generar_npys

"""
saca los outliers de la Rsd y Rtd y me los... grafica? pone en una lista? No decidí todavía
Usa los archivos: (si no están, los crea)
pos_bs.npy
pos_mpb.npy
newdates.npy

Genera los archivos:
pos_polar_bs.npy
pos_polar_mpb.npy
Rsd_bs.npy
Rtd_bs.npy
Rsd_mpb.npy
Rtd_mpb.npy
"""


def polarizar(pos, r0):
    """
    Le doy un array en cartesianas y me lo devuelve en polares
    R0 es el foco de la cónica
    """
    Rpolar = np.empty_like(pos)

    for i in range(len(pos)):
        # BS

        x, y, z = pos[i, 0], pos[i, 1], pos[i, 2]

        rho, theta, phi = fpos.cartesian2polar(x, y, z, r0)

        Rpolar[i, 0] = rho
        Rpolar[i, 1] = theta
        Rpolar[i, 2] = phi

    return Rpolar


def Rsd_Rtd(pos, pos_polar, x0, eps):
    """
    Encuentra la standoff distance y terminator distance dado un array de posiciones
    x0, eps son los parámetros de la cónica
    """

    N_events = len(pos[:, 0])
    L = np.empty(N_events)

    for i in range(N_events):
        L[i] = fvig.fit_L(pos_polar[i, 0], pos_polar[i, 1], eps)

    # CALCULATE VIGNES STANDOFF AND TERMINATOR DISTANCE

    Rsd = np.empty(N_events)
    Rtd = np.empty(N_events)

    for i in range(N_events):
        Rsd[i] = fvig.Rsd(L[i], x0, eps)

        Rtd[i] = fvig.Rtd(L[i], pos[i, 1], pos[i, 2], x0, eps)

    return Rsd, Rtd


def chired(ydata, ymod, param_fit):
    """
    Calcula el chi cuadrado reducido de un ajuste lineal respecto a los parámetros
    (ej: un ajuste lineal, que tiene param_fit = 2)
    """
    sigma_sq = np.mean(ydata**2) - (np.mean(ydata)) ** 2

    chi = np.sum(((ydata - ymod) ** 2) / sigma_sq)

    chi_red = chi / (len(ydata) - param_fit)

    return chi_red


g = input("número de grupo\n")

path = f"../outputs/grupo{g}/"

if not exists(path + "pos_bs.npy"):
    generar_npys(path)

pos_bs = np.transpose(np.load(path + "pos_bs.npy"))
pos_mpb = np.transpose(np.load(path + "pos_mpb.npy"))
newdates = np.load(path + "newdates.npy")

r0_BS = np.array([0.64, 0, 0])
Rpolar_BS = polarizar(pos_bs, r0_BS)
np.save(path + "pos_polar_bs.npy", Rpolar_BS)
r0_MPB = np.array([0.78, 0, 0])
Rpolar_MPB = polarizar(pos_mpb, r0_MPB)
np.save(path + "pos_polar_mpb.npy", Rpolar_MPB)
x0_bs = 0.64
eps_bs = 1.03
Rsd_BS, Rtd_BS = Rsd_Rtd(pos_bs, Rpolar_BS, x0_bs, eps_bs)
np.save(path + "Rsd_bs.npy", Rsd_BS)
np.save(path + "Rtd_bs.npy", Rtd_BS)
x0_mpb = 0.86
eps_mpb = 0.92
Rsd_MPB, Rtd_MPB = Rsd_Rtd(pos_mpb, Rpolar_MPB, x0_mpb, eps_mpb)
np.save(path + "Rsd_mpb.npy", Rsd_MPB)
np.save(path + "Rtd_mpb.npy", Rtd_MPB)


# recorto los outliers


idx = [i for i in range(len(Rtd_BS)) if np.abs(Rtd_BS[i]) < 10]

Rtd_BS_mod = np.array(Rtd_BS[idx])
Rsd_BS_mod = np.array(Rsd_BS[idx])
Rtd_MPB_mod = np.array(Rtd_MPB[idx])
Rsd_MPB_mod = np.array(Rsd_MPB[idx])

m_Rtd, b_Rtd = np.polyfit(Rtd_MPB_mod, Rtd_BS_mod, 1, cov=False)
m_Rsd, b_Rsd = np.polyfit(Rsd_MPB_mod, Rsd_BS_mod, 1, cov=False)

names = [i[0] + "-" + i[1] + "-" + i[2] + "h" + str(i[3])[:3] for i in newdates[idx]]
# calculate reduced chi square


chi_Rsd = np.round(chired(Rsd_BS_mod, m_Rsd * Rsd_MPB_mod + b_Rsd, 2), 3)

chi_Rtd = np.round(chired(Rtd_BS_mod, m_Rtd * Rtd_MPB_mod + b_Rtd, 2), 3)


# PLOT


fig = plt.figure(1)
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))
ax1.set_aspect("equal", "box")
ax2.set_aspect("equal", "box")
scatter_sd = ax1.scatter(Rsd_MPB_mod, Rsd_BS_mod)  # measured data
ax1.plot(
    Rsd_MPB_mod,
    m_Rsd * Rsd_MPB_mod + b_Rsd,
    color="red",
    label=r"$\chi^2=${}".format(chi_Rsd),
)  # model

scatter_td = ax2.scatter(Rtd_MPB_mod, Rtd_BS_mod)  # measured data
ax2.plot(
    Rtd_MPB_mod,
    m_Rtd * Rtd_MPB_mod + b_Rtd,
    color="red",
    label=r"$\chi^2=${}".format(chi_Rtd),
)  # model
ax1.set_xlabel(r"standoff distance MPB [$R_M$]")
ax1.set_ylabel(r"standoff distance BS [$R_M$]")
ax1.legend(loc=0)
ax2.set_xlabel(r"terminator distance MPB [$R_M$]")
ax2.set_ylabel(r"terminator distance BS [$R_M$]")
ax2.legend(loc=0)


annot = ax1.annotate(
    "",
    xy=(0, 0),
    xytext=(20, 20),
    textcoords="offset points",
    bbox=dict(boxstyle="round", fc="w"),
    arrowprops=dict(arrowstyle="->"),
)
annot.set_visible(False)

annot = ax2.annotate(
    "",
    xy=(0, 0),
    xytext=(20, 20),
    textcoords="offset points",
    bbox=dict(boxstyle="round", fc="w"),
    arrowprops=dict(arrowstyle="->"),
)
annot.set_visible(False)


def update_annot(ind, scatter):
    pos = scatter.get_offsets()[ind["ind"][0]]

    annot.xy = pos
    text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax1:
        cont_sd, ind_sd = scatter_sd.contains(event)
        if cont_sd:
            update_annot(ind_sd, scatter_sd)
            annot.set_visible(True)
            scatter_sd.set_facecolor("C0")  # Reset facecolor of all points in ax1
            fig.canvas.draw_idle()
            if ind_sd["ind"][0] < len(Rsd_BS):
                ind_td = {"ind": [ind_sd["ind"][0]]}
                update_annot(ind_td, scatter_td)
                scatter_td.set_facecolor(
                    ["red" if i in ind_td["ind"] else "C0" for i in range(len(Rtd_BS))]
                )  # Highlight the corresponding point in ax2
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()
    elif event.inaxes == ax2:
        cont_td, ind_td = scatter_td.contains(event)
        if cont_td:
            update_annot(ind_td, scatter_td)
            annot.set_visible(True)
            scatter_td.set_facecolor("C0")  # Reset facecolor of all points in ax2
            fig.canvas.draw_idle()
            if ind_td["ind"][0] < len(Rtd_MPB):
                ind_sd = {"ind": [ind_td["ind"][0]]}
                update_annot(ind_sd, scatter_sd)
                scatter_sd.set_facecolor(
                    ["red" if i in ind_sd["ind"] else "C0" for i in range(len(Rsd_MPB))]
                )  # Highlight the corresponding point in ax1
        else:
            if vis:
                annot.set_visible(False)
                scatter_sd.set_facecolor("C0")  # Reset facecolor of all points in ax1
                fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
