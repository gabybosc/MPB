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
Grafica discriminando por beta / SZA / etc
"""

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

path = f"../outputs/allgroups/"
lista = np.load(path + "lista.npy")
pos_bs = np.load(path + "pos_bs.npy")
pos_mpb = np.load(path + "pos_mpb.npy")
newdates = np.load(path + "newdates.npy")
pos_bs_min = np.load(path + "pos_bs_min.npy")
pos_bs_max = np.load(path + "pos_bs_max.npy")
pos_mpb_min = np.load(path + "pos_mpb_min.npy")
pos_mpb_max = np.load(path + "pos_mpb_max.npy")

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


plt.errorbar(Rsd_MPB, Rsd_BS, fmt="o", xerr=err_mpb, yerr=err_bs, zorder=0)
plt.legend()
plt.gca().set_aspect("equal")
plt.xlabel("Standoff Distance MPB (RM)")
plt.ylabel("Standoff Distance BS (RM)")
plt.title(f"Fit con el error mínimo")
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


ancho_MS = Rsd_BS - Rsd_MPB

"""
Discriminación por beta
"""


mag = [b for b in beta if b < 1]
dyn = [b for b in beta if 1 < b < 3]
dyn2 = [b for b in beta if 3 < b < 5]
dyn3 = [b for b in beta if b > 5]
idx_mag = [donde(beta, m) for m in mag]
idx_dyn = [donde(beta, m) for m in dyn]
idx_dyn2 = [donde(beta, m) for m in dyn2]
idx_dyn3 = [donde(beta, m) for m in dyn3]


fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle(r"Cruces divididos por $\beta$")

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

z1 = ax1.scatter(
    range(len(idx_mag)),
    ancho_MS[idx_mag],
    c=beta[idx_mag],
    vmin=0,
    vmax=1,
    label=r"$\beta < 1$",
)
z2 = ax2.scatter(
    range(len(idx_dyn)),
    ancho_MS[idx_dyn],
    c=beta[idx_dyn],
    vmin=1,
    vmax=3,
    label=r"$1 < \beta < 3$",
)
z3 = ax3.scatter(
    range(len(idx_dyn2)),
    ancho_MS[idx_dyn2],
    c=beta[idx_dyn2],
    vmin=3,
    vmax=5,
    label=r"$3 < \beta < 5$",
)
z4 = ax4.scatter(
    range(len(idx_dyn3)),
    ancho_MS[idx_dyn3],
    c=beta[idx_dyn3],
    vmin=5,
    vmax=10,
    label=r"$5 < \beta$",
)


for ax in [ax1, ax2, ax3, ax4]:
    ax.grid()
    # ax.set_aspect("equal", "box")
    # ax.set_xlim([1, 1.75])
    ax.set_ylim([-0.1, 1])
    ax.legend()

ax3.set_xlabel("cruces")
ax3.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
ax1.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z1, cax=cax)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z2, cax=cax)
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z3, cax=cax)
divider = make_axes_locatable(ax4)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z4, cax=cax)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
# plt.title("discriminacion por beta")
plt.show()


fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle(r"Cruces divididos por $\beta$")

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1), sharex=ax1, sharey=ax1)
ax3 = plt.subplot2grid((2, 2), (1, 0), sharex=ax1, sharey=ax1)
ax4 = plt.subplot2grid((2, 2), (1, 1), sharex=ax1, sharey=ax1)

z1 = ax1.scatter(
    Rsd_MPB[idx_mag],
    Rsd_BS[idx_mag],
    c=beta[idx_mag],
    vmin=0,
    vmax=1,
    label=r"$\beta < 1$",
)
z2 = ax2.scatter(
    Rsd_MPB[idx_dyn],
    Rsd_BS[idx_dyn],
    c=beta[idx_dyn],
    vmin=1,
    vmax=3,
    label=r"$1 < \beta < 3$",
)
z3 = ax3.scatter(
    Rsd_MPB[idx_dyn2],
    Rsd_BS[idx_dyn2],
    c=beta[idx_dyn2],
    vmin=3,
    vmax=5,
    label=r"$3 < \beta < 5$",
)
z4 = ax4.scatter(
    Rsd_MPB[idx_dyn3],
    Rsd_BS[idx_dyn3],
    c=beta[idx_dyn3],
    vmin=5,
    vmax=10,
    label=r"$5 < \beta$",
)


for ax in [ax1, ax2, ax3, ax4]:
    ax.grid()
    ax.legend()

ax3.set_xlabel("Rsd_MPB (RM)")
ax4.set_xlabel("Rsd_MPB (RM)")
ax3.set_ylabel("Rsd_BS (RM)")
ax1.set_ylabel("Rsd_BS (RM)")
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z1, cax=cax)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z2, cax=cax)
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z3, cax=cax)
divider = make_axes_locatable(ax4)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z4, cax=cax)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
# plt.title("discriminacion por beta")
plt.show()


"""
Discriminación por para/perp
"""

para = [t for t in theta if t < 30]
idx_para = [donde(theta, p) for p in para]
perp = [t for t in theta if t > 60]
idx_perp = [donde(theta, p) for p in perp]

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle("Cruces quasipara / quasiperp")

ax1 = plt.subplot2grid((2, 1), (0, 0))
ax2 = plt.subplot2grid((2, 1), (1, 0))

z1 = ax1.scatter(
    range(len(idx_para)), ancho_MS[idx_para], vmin=0, vmax=90, label=r"$\theta < 30$"
)

z2 = ax2.scatter(
    range(len(idx_perp)), ancho_MS[idx_perp], vmin=0, vmax=90, label=r"$\theta > 70$"
)

for ax in [ax1, ax2]:
    ax.legend()
    ax.set_ylim([-0.1, 1])
    ax.grid()
ax2.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
ax1.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
ax2.set_xlabel("cruces")

plt.show()


fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle("Cruces quasipara / quasiperp")

ax1 = plt.subplot2grid((2, 1), (0, 0))
ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1, sharey=ax1)

z1 = ax1.scatter(
    Rsd_MPB[idx_para], Rsd_BS[idx_para], vmin=0, vmax=90, label=r"$\theta < 30$"
)

z2 = ax2.scatter(
    Rsd_MPB[idx_perp], Rsd_BS[idx_perp], vmin=0, vmax=90, label=r"$\theta > 70$"
)

for ax in [ax1, ax2]:
    ax.legend()
    ax.grid()
ax2.set_ylabel("Rsd_BS (RM)")
ax1.set_ylabel("Rsd_BS (RM)")
ax2.set_xlabel("Rsd MPB (RM)")

plt.show()


"""
Discriminación por SZA
"""
sza_mpb = np.array([angulo(p, [1, 0, 0]) * 180 / np.pi for p in pos_mpb])
sza_bs = np.array([angulo(p, [1, 0, 0]) * 180 / np.pi for p in pos_bs])

sza1 = [t for t in sza_bs if t < 20]
sza2 = [t for t in sza_bs if 20 < t < 40]
sza3 = [t for t in sza_bs if 40 < t < 60]
sza4 = [t for t in sza_bs if 60 < t]

idx_sza1 = [donde(sza_bs, p) for p in sza1]
idx_sza2 = [donde(sza_bs, p) for p in sza2]
idx_sza3 = [donde(sza_bs, p) for p in sza3]
idx_sza4 = [donde(sza_bs, p) for p in sza4]

sza_mpb1 = [t for t in sza_mpb if t < 20]
sza_mpb2 = [t for t in sza_mpb if 20 < t < 40]
sza_mpb3 = [t for t in sza_mpb if 40 < t < 60]
sza_mpb4 = [t for t in sza_mpb if 60 < t]

idx_sza_mpb1 = [donde(sza_mpb, p) for p in sza_mpb1]
idx_sza_mpb2 = [donde(sza_mpb, p) for p in sza_mpb2]
idx_sza_mpb3 = [donde(sza_mpb, p) for p in sza_mpb3]
idx_sza_mpb4 = [donde(sza_mpb, p) for p in sza_mpb4]

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle("Plots según SZA BS - el color indica el SZA MPB")

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

z1 = ax1.scatter(
    range(len(idx_sza1)),
    ancho_MS[idx_sza1],
    c=sza_mpb[idx_sza1],
    vmin=0,
    vmax=90,
    label="SZA BS < 20",
)
z2 = ax2.scatter(
    range(len(idx_sza2)),
    ancho_MS[idx_sza2],
    c=sza_mpb[idx_sza2],
    vmin=0,
    vmax=90,
    label="20 < SZA BS < 40",
)
z3 = ax3.scatter(
    range(len(idx_sza3)),
    ancho_MS[idx_sza3],
    c=sza_mpb[idx_sza3],
    vmin=0,
    vmax=90,
    label="40 < SZA BS < 60",
)
z4 = ax4.scatter(
    range(len(idx_sza4)),
    ancho_MS[idx_sza4],
    c=sza_mpb[idx_sza4],
    vmin=0,
    vmax=90,
    label="60 < SZA BS",
)


for ax in [ax1, ax2, ax3, ax4]:
    ax.grid()
    ax.set_ylim([-0.1, 1])
    ax.legend()
ax3.set_xlabel("cruces")
ax3.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
ax1.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z1, cax=cax)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.show()


fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle("Plots según SZA BS - el color indica el SZA MPB")

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1), sharex=ax1, sharey=ax1)
ax3 = plt.subplot2grid((2, 2), (1, 0), sharex=ax1, sharey=ax1)
ax4 = plt.subplot2grid((2, 2), (1, 1), sharex=ax1, sharey=ax1)

z1 = ax1.scatter(
    Rsd_MPB[idx_sza1],
    Rsd_BS[idx_sza1],
    c=sza_mpb[idx_sza1],
    vmin=0,
    vmax=90,
    label="SZA BS < 20",
)
z2 = ax2.scatter(
    Rsd_MPB[idx_sza2],
    Rsd_BS[idx_sza2],
    c=sza_mpb[idx_sza2],
    vmin=0,
    vmax=90,
    label="20 < SZA BS < 40",
)
z3 = ax3.scatter(
    Rsd_MPB[idx_sza3],
    Rsd_BS[idx_sza3],
    c=sza_mpb[idx_sza3],
    vmin=0,
    vmax=90,
    label="40 < SZA BS < 60",
)
z4 = ax4.scatter(
    Rsd_MPB[idx_sza4],
    Rsd_BS[idx_sza4],
    c=sza_mpb[idx_sza4],
    vmin=0,
    vmax=90,
    label="60 < SZA BS",
)


for ax in [ax1, ax2, ax3, ax4]:
    ax.grid()
    ax.legend()
ax3.set_xlabel("Rsd MPB (RM)")
ax4.set_xlabel("Rsd MPB (RM)")
ax3.set_ylabel("Rsd_BS (RM)")
ax1.set_ylabel("Rsd_BS (RM)")
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(z1, cax=cax)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.show()


"""
Hemisferio
"""

idx_sur = []
idx_norte = []
for n, z in enumerate(pos_mpb):
    if z[2] < 0:
        idx_sur.append(n)
    else:
        idx_norte.append(n)


fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle("Cruces hemisferio")

ax1 = plt.subplot2grid((2, 1), (0, 0))
ax2 = plt.subplot2grid((2, 1), (1, 0))

z1 = ax1.scatter(range(len(idx_sur)), ancho_MS[idx_sur], vmin=0, vmax=90, label="Sur")

z2 = ax2.scatter(
    range(len(idx_norte)), ancho_MS[idx_norte], vmin=0, vmax=90, label="Norte"
)

for ax in [ax1, ax2]:
    ax.legend()
    ax.set_ylim([-0.1, 1])
    ax.grid()
ax2.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
ax1.set_ylabel("Rsd_BS - Rsd_MPB (RM)")
ax2.set_xlabel("cruces")

plt.show()

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

fig.suptitle("Cruces hemisferio")

ax1 = plt.subplot2grid((2, 1), (0, 0))
ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1, sharey=ax1)

z1 = ax1.scatter(Rsd_MPB[idx_sur], Rsd_BS[idx_sur], vmin=0, vmax=90, label="Sur")

z2 = ax2.scatter(Rsd_MPB[idx_norte], Rsd_BS[idx_norte], vmin=0, vmax=90, label="Norte")

for ax in [ax1, ax2]:
    ax.legend()
    ax.grid()
ax2.set_ylabel("Rsd_BS (RM)")
ax1.set_ylabel("Rsd_BS (RM)")
ax2.set_xlabel("Rsd MPB (RM)")

plt.show()
