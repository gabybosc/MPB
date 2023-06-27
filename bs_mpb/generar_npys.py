import numpy as np
import matplotlib.pyplot as plt
import sys
import plotly.express as px
from os.path import exists

sys.path.append("..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde

"""
Genera los archivos
pos_bs.npy
pos_mpb.npy
newdates.npy
newnew.npy
"""


def generar_npys(path):
    lista = np.genfromtxt(path + f"bs_mpb_final.txt", dtype=str)

    pos_bs = []  # la posición media del bs
    pos_mpb = []  # la posición media de la mpb
    newdates = []  # día y hora aprox del cruce
    final = []  # la lista original pero en npy

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

        pos_bs.append(pos[idx_bs] / 3390)
        pos_mpb.append(pos[idx_mpb] / 3390)
        newdates.append((year, month, day, t_bs))
        final.append(l)

    pos_bs = np.transpose(pos_bs)
    pos_mpb = np.transpose(pos_mpb)
    np.save(path + "newnew.npy", final)
    np.save(path + "newdates.npy", newdates)
    np.save(path + "pos_mpb.npy", pos_mpb)
    np.save(path + "pos_bs.npy", pos_bs)


# grupo = input("número de grupo\n")
# path = f"../outputs/grupo{grupo}/"
# generar_npys(grupo)
