import numpy as np
import matplotlib.pyplot as plt
from loader import importar_params, importar_posiciones
from scipy.stats import norm


path = "../../../../datos/bs_mpb/"
"""cone angle (última columna de los catalogos de Jacob actualizados), beta_protones, Mfms, Pdyn, Ls y bueno el Z"""
date, times, pos_bs, pos_mpb, Rsd = importar_posiciones(path)
beta, cone_angle, Mfms, Ls, pdyn = importar_params(path)

plt.figure()
plt.scatter(Rsd[1, :], Rsd[0, :], c=beta, vmin=0, vmax=10)
plt.ylabel(r"$R_{SD}$ BS")
plt.xlabel(r"$R_{SD}$ MPB")
cbar = plt.colorbar()
cbar.ax.set_title(r"$\beta_p$")
plt.show()

plt.figure()
plt.scatter(Rsd[1, :], Rsd[0, :], c=cone_angle, vmin=0, vmax=180)
plt.ylabel(r"$R_{SD}$ BS")
plt.xlabel(r"$R_{SD}$ MPB")
cbar = plt.colorbar()
cbar.ax.set_title("cone angle")
plt.show()

plt.figure()
plt.scatter(Rsd[1, :], Rsd[0, :], c=Mfms, vmin=1, vmax=9)
plt.ylabel(r"$R_{SD}$ BS")
plt.xlabel(r"$R_{SD}$ MPB")
cbar = plt.colorbar()
cbar.ax.set_title("Mfms")
plt.show()

plt.figure()
plt.scatter(Rsd[1, :], Rsd[0, :], c=pdyn, vmin=0, vmax=1)
plt.ylabel(r"$R_{SD}$ BS")
plt.xlabel(r"$R_{SD}$ MPB")
cbar = plt.colorbar()
cbar.ax.set_title(r"$P_{dyn}$")
plt.show()

plt.figure()
plt.scatter(Rsd[1, :], Rsd[0, :], c=Ls, vmin=0, vmax=31e4)
plt.ylabel(r"$R_{SD}$ BS")
plt.xlabel(r"$R_{SD}$ MPB")
cbar = plt.colorbar()
cbar.ax.set_title(r"$L_s$")
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

