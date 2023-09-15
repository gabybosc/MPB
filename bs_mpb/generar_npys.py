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
    lista = np.genfromtxt(path + f"FINAL_ESTA_SI.txt", dtype=str)

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

    pos_bs = np.array(pos_bs)
    pos_mpb = np.array(pos_mpb)
    np.save(path + "newnew.npy", final)
    np.save(path + "newdates.npy", newdates)
    np.save(path + "pos_mpb.npy", pos_mpb)
    np.save(path + "pos_bs.npy", pos_bs)


def generar_npys_limites(path):
    lista = np.genfromtxt(path + f"FINAL_ESTA_SI.txt", skip_header=1, dtype=str)

    pos_bs_min = []  # la posición min del bs
    pos_bs = []
    pos_bs_max = []
    pos_mpb_min = []  # la posición media del bs
    pos_mpb = []  # la posición media de la mpb
    pos_mpb_max = []
    newdates = []  # día y hora aprox del cruce
    final = []  # la lista original pero en npy


    # ate_bs	date_mpb	BS_min	BS	BS_max	flag_BS	MPB_min	MPB	MPB_max	flag_MPB	theta	beta
    # 0       1           2       3       4   5       6       7       8       9
    for l in lista:
        year, month, day = l[0].split("-")

        t_bs_min = UTC_to_hdec(l[2])
        t_bs = UTC_to_hdec(l[3])
        t_bs_max = UTC_to_hdec(l[4])
        t_mpb_min = UTC_to_hdec(l[6])
        t_mpb = UTC_to_hdec(l[7])
        t_mpb_max = UTC_to_hdec(l[8])

        if t_bs < t_mpb:
            ti = t_bs - 0.5
            tf = t_mpb + 0.5
        else:
            ti = t_mpb - 0.5
            tf = t_bs + 0.5
        if ti < 0:
            ti = 0
        if tf > 24:
            tf = 24

        mag, t, B, pos = importar_mag_1s(year, month, day, ti, tf)
        idx_bs = donde(t, t_bs)
        idx_bs_min = donde(t, t_bs_min)
        idx_bs_max = donde(t, t_bs_max)
        idx_mpb_min = donde(t, t_mpb_min)
        idx_mpb_max = donde(t, t_mpb_max)
        idx_mpb = donde(t, t_mpb)

        p = pos[idx_bs] / 3390
        pos_bs.append(p)
        # va a guardar como pos_min la que sea menor posición
        # no la que haya sido cruzada primero
        if np.linalg.norm(pos[idx_bs_min] / 3390) < np.linalg.norm(p):
            pos_bs_min.append(pos[idx_bs_min] / 3390)
            pos_bs_max.append(pos[idx_bs_max] / 3390)
        else:
            pos_bs_min.append(pos[idx_bs_max] / 3390)
            pos_bs_max.append(pos[idx_bs_min] / 3390)

        pp = pos[idx_mpb] / 3390
        pos_mpb.append(pp)
        if np.linalg.norm(pos[idx_mpb_min] / 3390) < np.linalg.norm(pp):
            pos_mpb_min.append(pos[idx_mpb_min] / 3390)
            pos_mpb_max.append(pos[idx_mpb_max] / 3390)
        else:
            pos_mpb_min.append(pos[idx_mpb_max] / 3390)
            pos_mpb_max.append(pos[idx_mpb_min] / 3390)
        newdates.append((year, month, day, t_bs))
        final.append(l)

    pos_bs = np.array(pos_bs)
    pos_bs_min = np.array(pos_bs_min)
    pos_bs_max = np.array(pos_bs_max)
    pos_mpb_min = np.array(pos_mpb_min)
    pos_mpb_max = np.array(pos_mpb_max)
    pos_mpb = np.array(pos_mpb)
    # return (pos_bs, pos_bs_max, pos_bs_min, pos_mpb, pos_bs_min, pos_mpb_max)
    np.save(path + "newnew.npy", final)
    np.save(path + "newdates.npy", newdates)
    np.save(path + "pos_bs.npy", pos_bs)
    np.save(path + "pos_bs_min.npy", pos_bs_min)
    np.save(path + "pos_bs_max.npy", pos_bs_max)
    np.save(path + "pos_mpb.npy", pos_mpb)
    np.save(path + "pos_mpb_min.npy", pos_mpb_min)
    np.save(path + "pos_mpb_max.npy", pos_mpb_max)


grupo = input("número de grupo\n")
path = f"../outputs/grupo{grupo}/"
generar_npys_limites(path)
