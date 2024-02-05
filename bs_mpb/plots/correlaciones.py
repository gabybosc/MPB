import numpy as np
import matplotlib.pyplot as plt
from loader import importar_params, importar_posiciones
from scipy.stats import norm


path = "../../../../datos/bs_mpb/"
"""cone angle (última columna de los catalogos de Jacob actualizados), beta_protones, Mfms, Pdyn, Ls y bueno el Z"""
date, times, pos_bs, pos_mpb, Rsd = importar_posiciones(path)
beta, cone_angle, Mfms, Ls, pdyn = importar_params(path)

# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, :], Rsd[0, :])
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.grid()
# plt.show()
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, :], Rsd[0, :], c=beta, vmin=0, vmax=6)
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.grid()
# cbar = plt.colorbar()
# cbar.ax.set_title(r"$\beta_p$")
# plt.show()
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, :], Rsd[0, :], c=cone_angle, vmin=0, vmax=180)
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.grid()
# cbar = plt.colorbar()
# cbar.ax.set_title("cone angle")
# plt.show()
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, :], Rsd[0, :], c=Mfms, vmin=1, vmax=9)
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# cbar = plt.colorbar()
# plt.grid()
# cbar.ax.set_title("Mfms")
# plt.show()
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, :], Rsd[0, :], c=pdyn, vmin=0, vmax=1)
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.grid()
# cbar = plt.colorbar()
# cbar.ax.set_title(r"$P_{dyn}$")
# plt.show()
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, :], Rsd[0, :], c=Ls, vmin=0, vmax=310)
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.grid()
# cbar = plt.colorbar()
# cbar.ax.set_title(r"$L_s$")
# plt.show()
#
# idx_min = [i for i in range(len(Ls)) if Ls[i] < 180]
# idx_max = [i for i in range(len(Ls)) if Ls[i] > 180]
#
# fig, ax = plt.subplots()
# ax.set_aspect("equal")
# plt.scatter(Rsd[1, idx_min], Rsd[0, idx_min], label="low wave season")
# plt.scatter(Rsd[1, idx_max], Rsd[0, idx_max], label="high wave season")
# plt.ylabel(r"$R_{SD}$ BS")
# plt.xlabel(r"$R_{SD}$ MPB")
# plt.legend()
# plt.grid()
# plt.show()

bs_sur = [i for i in range(len(pos_bs)) if pos_bs[i, 2] < 0]
bs_norte = [i for i in range(len(pos_bs)) if pos_bs[i, 2] > 0]
mpb_sur = [i for i in range(len(pos_mpb)) if pos_mpb[i, 2] < 0]
mpb_norte = [i for i in range(len(pos_mpb)) if pos_mpb[i, 2] > 0]

SS = list(set(bs_sur).intersection(mpb_sur))
NN = list(set(bs_norte).intersection(mpb_norte))
SN = list(set(bs_sur).intersection(mpb_norte))
NS = list(set(bs_norte).intersection(mpb_sur))

fig, ax = plt.subplots()
ax.set_aspect("equal")
plt.scatter(Rsd[1, SS], Rsd[0, SS], label="bs S mpb S")
plt.scatter(Rsd[1, NN], Rsd[0, NN], label="bs N mpb N")
plt.scatter(Rsd[1, SN], Rsd[0, SN], label="bs S mpb N")
plt.scatter(Rsd[1, NS], Rsd[0, NS], label="bs N mpb S")
plt.ylabel(r"$R_{SD}$ BS")
plt.xlabel(r"$R_{SD}$ MPB")
plt.legend()
plt.grid()
plt.show()

fig, ax = plt.subplots()
ax.set_aspect("equal")
plt.scatter(Rsd[1, mpb_norte], Rsd[0, mpb_norte], label="mpb N")
plt.scatter(Rsd[1, mpb_sur], Rsd[0, mpb_sur], label="mpb S")
plt.ylabel(r"$R_{SD}$ BS")
plt.xlabel(r"$R_{SD}$ MPB")
plt.legend()
plt.grid()
plt.show()


Rsd_low = np.vstack(
        (
            np.load("../../outputs/allgroups/Rsd_bs_low.npy").astype(float),
            np.load("../../outputs/allgroups/Rsd_mpb_low.npy").astype(float),
        )
    )
Rsd_h = np.vstack(
        (
            np.load("../../outputs/allgroups/Rsd_bs_high.npy").astype(float),
            np.load("../../outputs/allgroups/Rsd_mpb_high.npy").astype(float),
        )
    )


fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

ax1 = plt.subplot2grid((2, 1), (0, 0))
ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)

ax1.hist(Rsd_h[1, :], 20, density=1, alpha=0.5, label="higi")
ax1.hist(Rsd_low[1, :], 20, density=1, alpha=0.5, label="low")
(muB, sigmaB) = norm.fit(Rsd_h[1, :])
ax1.legend()
ax1.set_xlabel("Standoff Distance (RM) MPB")
ax1.set_ylabel("Cantidad de eventos")
ax2.hist(Rsd_h[0, :], 20, density=1, alpha=0.5, label="high")
ax2.hist(Rsd_low[0, :], 20, density=1, alpha=0.5, label="low")
ax2.set_xlabel("Standoff Distance (RM) BS")
ax2.set_ylabel("Cantidad de eventos")
plt.show()

