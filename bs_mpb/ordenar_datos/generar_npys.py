import numpy as np
import matplotlib.pyplot as plt
import sys
import plotly.express as px
from os.path import exists
import Vignesfit_functions as fvig
import func_position as fpos

sys.path.append("../..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde

"""
Genera los archivos
pos_bs.npy
pos_mpb.npy
newdates.npy
newnew.npy
"""


def polarizar(pos, r0):
    """
    Le doy un array en cartesianas y me lo devuelve en polares
    R0 es el foco de la cónica
    """
    Rpolar = np.empty_like(pos)

    for i in range(len(pos)):
        # BS

        x, y, z = pos[i, 0], pos[i, 1], pos[i, 2]

        rho, theta, phi = fpos.cartesian2polar(x, y, z, r0)

        Rpolar[i, 0] = rho
        Rpolar[i, 1] = theta
        Rpolar[i, 2] = phi

    return Rpolar


def Rsd_Rtd(pos, pos_polar, x0, eps):
    """
    Encuentra la standoff distance y terminator distance dado un array de posiciones
    x0, eps son los parámetros de la cónica
    """

    N_events = len(pos[:, 0])
    L = np.empty(N_events)

    for i in range(N_events):
        L[i] = fvig.fit_L(pos_polar[i, 0], pos_polar[i, 1], eps)

    # CALCULATE VIGNES STANDOFF AND TERMINATOR DISTANCE

    Rsd = np.empty(N_events)
    Rtd = np.empty(N_events)

    for i in range(N_events):
        Rsd[i] = fvig.Rsd(L[i], x0, eps)

        Rtd[i] = fvig.Rtd(L[i], x0, eps)

    return Rsd, Rtd


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
    lista = np.genfromtxt(path + f"catalogo_final.txt", skip_header=1, dtype=str)

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
        print(l)
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


def gen_Rsd(path, BS, MPB):
    pos_bs_min, pos_bs, pos_bs_max = BS
    pos_mpb_min, pos_mpb, pos_mpb_max = MPB
    r0_BS = np.array([0.64, 0, 0])
    Rpolar_BS_min = polarizar(pos_bs_min, r0_BS)
    Rpolar_BS_max = polarizar(pos_bs_max, r0_BS)
    Rpolar_BS = polarizar(pos_bs, r0_BS)
    np.save(path + "pos_polar_bs_min.npy", Rpolar_BS_min)
    np.save(path + "pos_polar_bs_max.npy", Rpolar_BS_max)
    np.save(path + "pos_polar_bs.npy", Rpolar_BS)

    r0_MPB = np.array([0.78, 0, 0])
    Rpolar_MPB_min = polarizar(pos_mpb_min, r0_MPB)
    Rpolar_MPB_max = polarizar(pos_mpb_max, r0_MPB)
    Rpolar_MPB = polarizar(pos_mpb, r0_MPB)
    np.save(path + "pos_polar_mpb_min.npy", Rpolar_MPB_min)
    np.save(path + "pos_polar_mpb_max.npy", Rpolar_MPB_max)
    np.save(path + "pos_polar_mpb.npy", Rpolar_MPB)

    x0_bs = 0.64
    eps_bs = 1.03
    Rsd_BS_min, Rtd_BS_min = Rsd_Rtd(pos_bs_min, Rpolar_BS_min, x0_bs, eps_bs)
    Rsd_BS_max, Rtd_BS_max = Rsd_Rtd(pos_bs_max, Rpolar_BS_max, x0_bs, eps_bs)
    Rsd_BS, Rtd_BS = Rsd_Rtd(pos_bs, Rpolar_BS, x0_bs, eps_bs)
    np.save(path + "Rsd_bs_min.npy", Rsd_BS_min)
    np.save(path + "Rtd_bs_min.npy", Rtd_BS_min)
    np.save(path + "Rsd_bs_max.npy", Rsd_BS_max)
    np.save(path + "Rtd_bs_max.npy", Rtd_BS_max)
    np.save(path + "Rsd_bs.npy", Rsd_BS)
    np.save(path + "Rtd_bs.npy", Rtd_BS)

    x0_mpb = 0.78
    eps_mpb = 0.90
    Rsd_MPB_min, Rtd_MPB_min = Rsd_Rtd(pos_mpb_min, Rpolar_MPB_min, x0_mpb, eps_mpb)
    Rsd_MPB_max, Rtd_MPB_max = Rsd_Rtd(pos_mpb_max, Rpolar_MPB_max, x0_mpb, eps_mpb)
    Rsd_MPB, Rtd_MPB = Rsd_Rtd(pos_mpb, Rpolar_MPB, x0_mpb, eps_mpb)
    np.save(path + "Rsd_mpb_min.npy", Rsd_MPB_min)
    np.save(path + "Rtd_mpb_min.npy", Rtd_MPB_min)
    np.save(path + "Rsd_mpb_max.npy", Rsd_MPB_max)
    np.save(path + "Rtd_mpb_max.npy", Rtd_MPB_max)
    np.save(path + "Rsd_mpb.npy", Rsd_MPB)
    np.save(path + "Rtd_mpb.npy", Rtd_MPB)
