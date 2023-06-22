import numpy as np
import matplotlib.pyplot as plt
import sys
import plotly.express as px
from os.path import exists

sys.path.append("..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde


path = f"../outputs/grupo1/"
bs1 = np.load(path + "pos_bs.npy")
mpb1 = np.load(path + "pos_mpb.npy")
fechas1 = np.load(path + "newdates.npy")

path = f"../outputs/grupo2/"
bs2 = np.load(path + "pos_bs.npy")
mpb2 = np.load(path + "pos_mpb.npy")
fechas2 = np.load(path + "newdates.npy")

path = f"../outputs/grupo3/"
bs3 = np.load(path + "pos_bs.npy")
mpb3 = np.load(path + "pos_mpb.npy")
fechas3 = np.load(path + "newdates.npy")

path = f"../outputs/grupo4/"
bs4 = np.load(path + "pos_bs.npy")
mpb4 = np.load(path + "pos_mpb.npy")
fechas4 = np.load(path + "newdates.npy")


pos_bs = np.hstack((bs2, bs3, bs4))
pos_mpb = np.hstack((mpb2, mpb3, mpb4))
newdates = np.vstack((fechas2, fechas3, fechas4))

# from plot_orbitas import marte, BS_MPB
# X_bs, YZ_bs = BS_MPB(2.04, 1.03, 0.64)
# X_mpb, YZ_mpb = BS_MPB(0.96, 0.9, 0.78)
# # marte(ax, X_bs, YZ_bs, X_mpb, YZ_mpb)

x_bs = pos_bs[0]
yz_bs = np.sqrt(pos_bs[1] ** 2 + pos_bs[2] ** 2)
x_mpb = pos_mpb[0]
yz_mpb = np.sqrt(pos_mpb[1] ** 2 + pos_mpb[2] ** 2)
names = [i[0] + "-" + i[1] + "-" + i[2] + "h" + str(i[3])[:3] for i in newdates]

fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.set_xlim(0, 3)
ax.set_ylim(0, 2.5)
circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
ax.add_artist(circle)
ax.set_title("MAVEN MSO coordinates", fontsize=16)
ax.set_xlabel(r"$X_{MSO}$ ($R_M$)", fontsize=14)
ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)", fontsize=14)

scatter_bs = ax.scatter(x_bs, yz_bs)
scatter_mpb = ax.scatter(x_mpb, yz_mpb)

annot = ax.annotate(
    "",
    xy=(0, 0),
    xytext=(20, 20),
    textcoords="offset points",
    bbox=dict(boxstyle="round", fc="w"),
    arrowprops=dict(arrowstyle="->"),
)
annot.set_visible(False)


def update_annot(ind, scatter):
    if scatter == scatter_bs:
        pos = scatter_bs.get_offsets()[ind["ind"][0]]
    elif scatter == scatter_mpb:
        pos = scatter_mpb.get_offsets()[ind["ind"][0]]
    else:
        return

    annot.xy = pos
    text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont_bs, ind_bs = scatter_bs.contains(event)
        cont_mpb, ind_mpb = scatter_mpb.contains(event)
        if cont_bs:
            update_annot(ind_bs, scatter_bs)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        elif cont_mpb:
            update_annot(ind_mpb, scatter_mpb)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
