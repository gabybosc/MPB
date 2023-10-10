import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
from generar_npys import generar_npys_limites, gen_Rsd
from scipy import odr
from scipy.stats import chisquare
import matplotlib as mpl
from cycler import cycler
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

sys.path.append("..")
from funciones import donde, angulo

"""
Grafica R_sd BS vs MPB teniendo en cuenta el error
hace un fit de deming https://en.wikipedia.org/wiki/Deming_regression
que es lo mejor que encontré para que tenga en cuenta los errores
"""

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)


def chired(ydata, ymod, param_fit):
    """
    Calcula el chi cuadrado reducido de un ajuste lineal respecto a los parámetros
    (ej: un ajuste lineal, que tiene param_fit = 2)
    """
    sigma_sq = np.mean(ydata**2) - (np.mean(ydata)) ** 2

    chi = np.sum(((ydata - ymod) ** 2) / sigma_sq)

    chi_red = chi / (len(ydata) - param_fit)

    return np.round(chi_red, 3)


g = input("número de grupo\n")

path = f"../outputs/grupo{g}/"

if not exists(path + "pos_bs_min.npy"):
    generar_npys_limites(path)

lista = np.load(path + "newnew.npy")  # la lista con todo pero en numpy
pos_bs = np.load(path + "pos_bs.npy")
pos_bs_min = np.load(path + "pos_bs_min.npy")
pos_bs_max = np.load(path + "pos_bs_max.npy")
pos_mpb = np.load(path + "pos_mpb.npy")
pos_mpb_min = np.load(path + "pos_mpb_min.npy")
pos_mpb_max = np.load(path + "pos_mpb_max.npy")
newdates = np.load(path + "newdates.npy")

beta = np.array([float(l) for l in lista[:, -1]])
theta = np.array([float(l) for l in lista[:, -2]])

if not exists(path + "pos_polar_bs_min.npy"):
    BS = [pos_bs_min, pos_bs, pos_bs_max]
    MPB = [pos_mpb_min, pos_mpb, pos_mpb_max]
    gen_Rsd(path, BS, MPB)

Rsd_MPB_min = np.load(path + "Rsd_mpb_min.npy")
Rsd_MPB_max = np.load(path + "Rsd_mpb_max.npy")
Rsd_MPB = np.load(path + "Rsd_mpb.npy")
Rsd_BS_min = np.load(path + "Rsd_bs_min.npy")
Rsd_BS_max = np.load(path + "Rsd_bs_max.npy")
Rsd_BS = np.load(path + "Rsd_bs.npy")

Rtd_MPB_min = np.load(path + "Rtd_mpb_min.npy")
Rtd_MPB_max = np.load(path + "Rtd_mpb_max.npy")
Rtd_MPB = np.load(path + "Rtd_mpb.npy")
Rtd_BS_min = np.load(path + "Rtd_bs_min.npy")
Rtd_BS_max = np.load(path + "Rtd_bs_max.npy")
Rtd_BS = np.load(path + "Rtd_bs.npy")


# ojo, los errores son el tamaño, no la posición, entonces tengo que restarlo
# de la posición media
err_bs = np.array((np.abs(Rsd_BS_min - Rsd_BS), np.abs(Rsd_BS_max - Rsd_BS)))
err_mpb = np.array((np.abs(Rsd_MPB_min - Rsd_MPB), np.abs(Rsd_MPB_max - Rsd_MPB)))
# plt.errorbar(Rsd_MPB, Rsd_BS, fmt="o", xerr=err_mpb, yerr=err_bs)
# plt.gca().set_aspect("equal")
# plt.show()

# fiteo la data sin errores para tener coefs que alimentarle al odr
coef = np.polyfit(Rsd_MPB, Rsd_BS, deg=1, w=1 / np.abs(Rsd_BS_min - Rsd_BS))
lin = coef[0] * Rsd_MPB + coef[1]


def f(B, x):
    """Linear function y = m*x + b"""
    return B[0] * x + B[1]


linear = odr.Model(f)

mydata = odr.RealData(
    Rsd_MPB, Rsd_BS, sx=err_mpb[0], sy=err_bs[0]
)  # sx, sy tienen que ser 1D
myodr = odr.ODR(mydata, linear, beta0=[coef[0], coef[1]])
myoutput = myodr.run()
myoutput.pprint()


newfit = myoutput.beta[0] * Rsd_MPB + myoutput.beta[1]
print(f"pendiente del fit c/error {myoutput.beta[0]}")
# chisquare(Rsd_BS, newfit)
# chisquare(Rsd_BS, lin)
chi_Rsd = chired(Rsd_BS, lin, 2)
chi_Rsd_err = chired(Rsd_BS, newfit, 2)
plt.errorbar(Rsd_MPB, Rsd_BS, fmt="o", xerr=err_mpb, yerr=err_bs, zorder=0)
plt.plot(Rsd_MPB, lin, label=f"w/o error, {chi_Rsd}")
plt.plot(Rsd_MPB, newfit, label=f"w/error, {chi_Rsd_err}")
plt.legend()
plt.gca().set_aspect("equal")
plt.xlabel("Standoff Distance MPB (RM)")
plt.ylabel("Standoff Distance BS (RM)")
plt.title(f"Fit con el error mínimo - Grupo {g}")
plt.show()


def plot_error(idx, lbl):
    plt.errorbar(
        Rsd_MPB[idx],
        Rsd_BS[idx],
        fmt="o",
        xerr=err_mpb[:, idx],
        yerr=err_bs[:, idx],
        label=lbl,
    )


"""
Discriminación por beta
"""


# mag = [b for b in beta if b < 1]
# dyn = [b for b in beta if 1 < b < 3]
# dyn2 = [b for b in beta if 3 < b < 5]
# dyn3 = [b for b in beta if b > 5]
# idx_mag = [donde(beta, m) for m in mag]
# idx_dyn = [donde(beta, m) for m in dyn]
# idx_dyn2 = [donde(beta, m) for m in dyn2]
# idx_dyn3 = [donde(beta, m) for m in dyn3]

# plot_error(idx_mag, r"$\beta < 1$")
# plot_error(idx_dyn, r"$1 < \beta < 3$")
# plot_error(idx_dyn2, r"$3 < \beta < 5$")
# plot_error(idx_dyn3, r"$\beta > 5$")
# plt.legend()
# plt.grid()
# plt.gca().set_aspect("equal")
# plt.xlabel("Standoff Distance MPB (RM)")
# plt.ylabel("Standoff Distance BS (RM)")
# plt.title(f"Discriminación por beta - Grupo {g}")
# plt.show()

# con un colormap si no
plt.figure()
plt.errorbar(
    Rsd_MPB,
    Rsd_BS,
    fmt="o",
    xerr=err_mpb,
    yerr=err_bs,
    c="k",
    zorder=0,
    # label=lbl,
)
plt.scatter(
    Rsd_MPB,
    Rsd_BS,
    c=beta,
    vmin=min(beta),
    vmax=5,
)
plt.colorbar()
plt.grid()
plt.gca().set_aspect("equal")
plt.xlabel("Standoff Distance MPB (RM)")
plt.ylabel("Standoff Distance BS (RM)")
plt.title(f"Discriminación por beta (cap en 5) - Grupo {g}")
plt.savefig(path + "Discrimina_beta.png", dpi=300)
# plt.show()

"""
Discriminación por para/perp
"""

para = [t for t in theta if t < 30]
idx_para = [donde(theta, p) for p in para]
perp = [t for t in theta if t > 60]
idx_perp = [donde(theta, p) for p in perp]

plt.figure()
plot_error(idx_para, r"$\theta < 30$")
plot_error(idx_perp, r"$\theta > 70$")

plt.legend()
plt.grid()
plt.gca().set_aspect("equal")
plt.xlabel("Standoff Distance MPB (RM)")
plt.ylabel("Standoff Distance BS (RM)")
plt.title(f"Discriminación por theta - Grupo {g}")
plt.savefig(path + "Discrimina_theta.png", dpi=300)
# plt.show()


"""
Discriminación por SZA
"""
sza_mpb = np.array([angulo(p, [1, 0, 0]) * 180 / np.pi for p in pos_mpb])
sza_bs = np.array([angulo(p, [1, 0, 0]) * 180 / np.pi for p in pos_bs])

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))
ax1.set_aspect("equal", "box")
ax2.set_aspect("equal", "box")
ax1.errorbar(
    Rsd_MPB,
    Rsd_BS,
    fmt="o",
    xerr=err_mpb,
    yerr=err_bs,
    c="k",
    zorder=0,
)
z1 = ax1.scatter(Rsd_MPB, Rsd_BS, c=sza_mpb, vmin=0, vmax=90)
ax2.errorbar(
    Rsd_MPB,
    Rsd_BS,
    fmt="o",
    xerr=err_mpb,
    yerr=err_bs,
    c="k",
    zorder=0,
)
ax2.scatter(Rsd_MPB, Rsd_BS, c=sza_bs, vmin=0, vmax=90)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z1, cax=cax)
plt.grid()
plt.gca().set_aspect("equal")
for ax in [ax1, ax2]:
    ax.set_xlabel("Standoff Distance MPB (RM)")
    ax.set_ylabel("Standoff Distance BS (RM)")
    ax.grid()
ax1.set_title(f"Discriminación por SZA MPB- Grupo {g}")
ax2.set_title(f"Discriminación por SZA BS- Grupo {g}")
plt.savefig(path + "Discrimina_SZA.png", dpi=300)
# plt.show()


"""
Restamos Rsd de MPB - BS: si son muy parecidas o muy diferentes, es que el fit de vignes no está funcionando.
"""
plt.figure()
plt.plot(Rsd_MPB, Rsd_BS - Rsd_MPB, ".")
plt.xlabel("Standoff Distance MPB (RM)")
plt.ylabel("Difference (RM)")
plt.title(f"Difference in Standoff Distance BS - MPB (grupo {g})")
plt.grid()
plt.savefig(path + "Resta_RSD.png", dpi=300)
plt.show()
# recorto los outliers  # esto me falta modificarlo por el momento


# idx = [i for i in range(len(Rtd_BS)) if np.abs(Rtd_BS[i]) < 10]

# Rtd_BS_mod = np.array(Rtd_BS[idx])
# Rsd_BS_mod = np.array(Rsd_BS[idx])
# Rtd_MPB_mod = np.array(Rtd_MPB[idx])
# Rsd_MPB_mod = np.array(Rsd_MPB[idx])

# m_Rtd, b_Rtd = np.polyfit(Rtd_MPB_mod, Rtd_BS_mod, 1, cov=False)
# m_Rsd, b_Rsd = np.polyfit(Rsd_MPB_mod, Rsd_BS_mod, 1, cov=False)

# names = [i[0] + "-" + i[1] + "-" + i[2] + "h" + str(i[3])[:3] for i in newdates[idx]]
# # calculate reduced chi square


# chi_Rsd = np.round(chired(Rsd_BS_mod, m_Rsd * Rsd_MPB_mod + b_Rsd, 2), 3)

# chi_Rtd = np.round(chired(Rtd_BS_mod, m_Rtd * Rtd_MPB_mod + b_Rtd, 2), 3)


# # PLOT


# fig = plt.figure(1)
# fig.subplots_adjust(
#     top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
# )
# ax1 = plt.subplot2grid((1, 2), (0, 0))
# ax2 = plt.subplot2grid((1, 2), (0, 1))
# ax1.set_aspect("equal", "box")
# ax2.set_aspect("equal", "box")
# scatter_sd = ax1.scatter(Rsd_MPB_mod, Rsd_BS_mod)  # measured data
# ax1.plot(
#     Rsd_MPB_mod,
#     m_Rsd * Rsd_MPB_mod + b_Rsd,
#     color="red",
#     label=r"$\chi^2=${}".format(chi_Rsd),
# )  # model

# scatter_td = ax2.scatter(Rtd_MPB_mod, Rtd_BS_mod)  # measured data
# ax2.plot(
#     Rtd_MPB_mod,
#     m_Rtd * Rtd_MPB_mod + b_Rtd,
#     color="red",
#     label=r"$\chi^2=${}".format(chi_Rtd),
# )  # model
# ax1.set_xlabel(r"standoff distance MPB [$R_M$]")
# ax1.set_ylabel(r"standoff distance BS [$R_M$]")
# ax1.legend(loc=0)
# ax2.set_xlabel(r"terminator distance MPB [$R_M$]")
# ax2.set_ylabel(r"terminator distance BS [$R_M$]")
# ax2.legend(loc=0)


# annot = ax1.annotate(
#     "",
#     xy=(0, 0),
#     xytext=(20, 20),
#     textcoords="offset points",
#     bbox=dict(boxstyle="round", fc="w"),
#     arrowprops=dict(arrowstyle="->"),
# )
# annot.set_visible(False)

# annot = ax2.annotate(
#     "",
#     xy=(0, 0),
#     xytext=(20, 20),
#     textcoords="offset points",
#     bbox=dict(boxstyle="round", fc="w"),
#     arrowprops=dict(arrowstyle="->"),
# )
# annot.set_visible(False)


# def update_annot(ind, scatter):
#     pos = scatter.get_offsets()[ind["ind"][0]]

#     annot.xy = pos
#     text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
#     annot.set_text(text)
#     annot.get_bbox_patch().set_alpha(0.4)


# def hover(event):
#     vis = annot.get_visible()
#     if event.inaxes == ax1:
#         cont_sd, ind_sd = scatter_sd.contains(event)
#         if cont_sd:
#             update_annot(ind_sd, scatter_sd)
#             annot.set_visible(True)
#             scatter_sd.set_facecolor("C0")  # Reset facecolor of all points in ax1
#             fig.canvas.draw_idle()
#             if ind_sd["ind"][0] < len(Rsd_BS):
#                 ind_td = {"ind": [ind_sd["ind"][0]]}
#                 update_annot(ind_td, scatter_td)
#                 scatter_td.set_facecolor(
#                     ["red" if i in ind_td["ind"] else "C0" for i in range(len(Rtd_BS))]
#                 )  # Highlight the corresponding point in ax2
#         else:
#             if vis:
#                 annot.set_visible(False)
#                 fig.canvas.draw_idle()
#     elif event.inaxes == ax2:
#         cont_td, ind_td = scatter_td.contains(event)
#         if cont_td:
#             update_annot(ind_td, scatter_td)
#             annot.set_visible(True)
#             scatter_td.set_facecolor("C0")  # Reset facecolor of all points in ax2
#             fig.canvas.draw_idle()
#             if ind_td["ind"][0] < len(Rtd_MPB):
#                 ind_sd = {"ind": [ind_td["ind"][0]]}
#                 update_annot(ind_sd, scatter_sd)
#                 scatter_sd.set_facecolor(
#                     ["red" if i in ind_sd["ind"] else "C0" for i in range(len(Rsd_MPB))]
#                 )  # Highlight the corresponding point in ax1
#         else:
#             if vis:
#                 annot.set_visible(False)
#                 scatter_sd.set_facecolor("C0")  # Reset facecolor of all points in ax1
#                 fig.canvas.draw_idle()


# fig.canvas.mpl_connect("motion_notify_event", hover)

# plt.show()
