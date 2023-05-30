import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde


lista = np.genfromtxt("../outputs/grupo1.txt", dtype=str)
fig_path = "../../Pictures/BS_MPB/"

pos_bs = []
pos_mpb = []


for l in lista:
    year, month, day = l[0].split("-")

    t_bs = UTC_to_hdec(l[1])
    t_mpb = UTC_to_hdec(l[2])

    if t_bs < t_mpb:
        ti = t_bs - 0.2
        tf = t_mpb + 0.2
    else:
        ti = t_mpb - 0.2
        tf = t_bs + 0.2
    if ti < 0:
        ti = 0
    if tf > 24:
        tf = 24

    mag, t, B, pos = importar_mag_1s(year, month, day, ti, tf)
    idx_bs = donde(t, t_bs)
    idx_mpb = donde(t, t_mpb)

    pos_bs.append(pos[idx_bs])
    pos_mpb.append(pos[idx_mpb])


from plot_orbitas import marte, orbitas, BS_MPB

x_bs, yz_bs = BS_MPB(2.04, 1.03, 0.64)
x_mpb, yz_mpb = BS_MPB(0.96, 0.9, 0.78)
marte(x_bs, yz_bs, x_mpb, yz_mpb)
x_bs = np.transpose(pos_bs)[0] / 3390
yz_bs = np.sqrt(
    (np.transpose(pos_bs)[1] / 3390) ** 2 + (np.transpose(pos_bs)[2] / 3390) ** 2
)
x_mpb = np.transpose(pos_mpb)[0] / 3390
yz_mpb = np.sqrt(
    (np.transpose(pos_mpb)[1] / 3390) ** 2 + (np.transpose(pos_mpb)[2] / 3390) ** 2
)
plt.scatter(x_bs, yz_bs)
plt.scatter(x_mpb, yz_mpb)
plt.show()
# orbitas(posicion_cut, year, month, day)

# t, posicion_cut, year, month, day = datos()
# t, posicion_cut, year, month, day = datos_fijos(2015, 10, 10, 12, 13)
# orbitas(posicion_cut, year, month, day)

# t, posicion_cut, year, month, day = datos_fijos(2015, 10, 12, 18.75, 19.75)
# orbitas(posicion_cut, year, month, day)
