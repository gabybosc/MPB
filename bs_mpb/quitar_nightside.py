import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde

lista = np.genfromtxt("../outputs/new_grupo4.txt", dtype=str)

pos_bs = []
pos_mpb = []
newdates = []
final = []

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

    if pos[idx_mpb][0] > 0:
        pos_bs.append(pos[idx_bs] / 3390)
        pos_mpb.append(pos[idx_mpb] / 3390)
        newdates.append((year, month, day))
        final.append(l)


pos_bs = np.transpose(pos_bs)
pos_mpb = np.transpose(pos_mpb)

from plot_orbitas import marte, BS_MPB

x_bs, yz_bs = BS_MPB(2.04, 1.03, 0.64)
x_mpb, yz_mpb = BS_MPB(0.96, 0.9, 0.78)
marte(x_bs, yz_bs, x_mpb, yz_mpb)
x_bs = pos_bs[0]
yz_bs = np.sqrt(pos_bs[1] ** 2 + pos_bs[2] ** 2)
x_mpb = pos_mpb[0]
yz_mpb = np.sqrt(pos_mpb[1] ** 2 + pos_mpb[2] ** 2)
plt.scatter(x_bs, yz_bs)
plt.scatter(x_mpb, yz_mpb)
plt.show()

with open("../outputs/newnew_grupo4.txt", "a") as file:
    for f in final:
        file.write(f"{f[0]}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}")
        file.write("\n")
