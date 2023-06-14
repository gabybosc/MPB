import numpy as np
import matplotlib.pyplot as plt
import sys
import plotly.express as px
from os.path import exists

sys.path.append("..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde

"""
Un plot con los cuatro grupos en diferentes subplots
Son los puntos de los cruces en torno al planeta
"""


fig = plt.figure(1, constrained_layout=True)
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)

for g in [1, 2, 3, 4]:
    path = f"../outputs/grupo{g}/"
    pos_bs = np.load(path + "pos_bs.npy")
    pos_mpb = np.load(path + "pos_mpb.npy")
    newdates = np.load(path + "newdates.npy")

    if g == 1:
        ax = plt.subplot2grid((2, 2), (0, 0))
    elif g == 2:
        ax = plt.subplot2grid((2, 2), (1, 0))
    elif g == 3:
        ax = plt.subplot2grid((2, 2), (0, 1))
    elif g == 4:
        ax = plt.subplot2grid((2, 2), (1, 1))

    x_bs = pos_bs[0]
    yz_bs = np.sqrt(pos_bs[1] ** 2 + pos_bs[2] ** 2)
    x_mpb = pos_mpb[0]
    yz_mpb = np.sqrt(pos_mpb[1] ** 2 + pos_mpb[2] ** 2)
    names = [i[0] + "-" + i[1] + "-" + i[2] + "h" + str(i[3])[:3] for i in newdates]

    ax.set_aspect("equal")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 2.5)
    circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
    ax.add_artist(circle)
    ax.set_title(f"Group {g}", fontsize=16)
    if g == 2 or g == 4:
        ax.set_xlabel(r"$X_{MSO}$ ($R_M$)", fontsize=14)
    if g == 1 or g == 2:
        ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)", fontsize=14)
    if g == 1 or g == 3:
        plt.setp(ax.get_xticklabels(), visible=False)

    scatter_bs = ax.scatter(x_bs, yz_bs)
    scatter_mpb = ax.scatter(x_mpb, yz_mpb)
    scatter_bs = ax.scatter(x_bs[1], yz_bs[1], c="C2", marker="s")
    scatter_mpb = ax.scatter(x_mpb[1], yz_mpb[1], c="C2", marker="s")
    scatter_bs = ax.scatter(x_bs[15], yz_bs[15], c="C3", marker="d")
    scatter_mpb = ax.scatter(x_mpb[15], yz_mpb[15], c="C3", marker="d")
    print(g, newdates[1], newdates[15])


plt.show()
