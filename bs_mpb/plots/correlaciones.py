import numpy as np
import matplotlib.pyplot as plt
from loader import importar_params, importar_posiciones

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
