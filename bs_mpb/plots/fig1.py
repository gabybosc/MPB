import numpy as np
import matplotlib.pyplot as plt
import sys
import plotly.express as px
from os.path import exists
from loader import importar_params, importar_posiciones

sys.path.append("../..")
path = "../../../../datos/bs_mpb/"
"""cone angle (última columna de los catalogos de Jacob actualizados), beta_protones, Mfms, Pdyn, Ls y bueno el Z"""
date, times, pos_bs, pos_mpb, Rsd = importar_posiciones(path)
beta, cone_angle, Mfms, Ls, pdyn = importar_params(path)

x_bs = pos_bs[0]
yz_bs = np.sqrt(pos_bs[1] ** 2 + pos_bs[2] ** 2)
x_mpb = pos_mpb[0]
yz_mpb = np.sqrt(pos_mpb[1] ** 2 + pos_mpb[2] ** 2)

fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.set_xlim(0, 2.5)
ax.set_ylim(0, 2.5)
circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
ax.add_artist(circle)
ax.set_title("MAVEN MSO coordinates", fontsize=16)
ax.set_xlabel(r"$X_{MSO}$ ($R_M$)", fontsize=14)
ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)", fontsize=14)
scatter_bs = ax.scatter(x_bs, yz_bs, label="BS", c="#FE6779")
scatter_mpb = ax.scatter(x_mpb, yz_mpb, label="MPB", c="#79B953")
plt.legend()
plt.show()
