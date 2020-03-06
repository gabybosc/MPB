import numpy as np
import os
import sys

sys.path.append("..")
from funciones import find_nearest


def importar_mag(year, month, day, ti, tf):
    # path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
    path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la desktop.
    # Estaría bueno ponerle un if para que detecte en cuál estoy.
    if os.path.isfile(path + "mag_filtrado.txt"):
        mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    else:
        mag = np.loadtxt(path + "MAG.asc")
        B = mag[:, 6:9]

    hh = mag[:, 3]
    mm = mag[:, 4]
    ss = mag[:, 5]

    t = hh + mm / 60 + ss / 3600  # hdec

    posicion = np.zeros((len(t), 3))
    for i in range(9, 12):
        posicion[:, i - 9] = mag[:, i]

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    posicion_cut = posicion[inicio:fin]
    return mag, t_cut, B_cut, posicion_cut


# ##############################################################################################SWEA
def importar_swea(year, month, day, ti, tf):
    # path = f'../../../datos/clweb/{year}-{month}-{day}/'
    path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    swea = np.loadtxt(path + "SWEA.asc")

    energy = swea[:, 7]
    JE_total = swea[:, -1]

    t = np.unique(swea[:, 3] + swea[:, 4] / 60 + swea[:, 5] / 3600)  # hdec

    energias = [50 + i * 50 for i in range(3)]

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    energy_cut = energy[inicio:fin]
    JE_cut = JE_total[inicio:fin]

    return swea, t, energias, energy_cut, JE_cut


# ##############################################################################################SWIA


def importar_swia(year, month, day, ti, tf):
    # path = f'../../../datos/clweb/{year}-{month}-{day}/'
    path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    swia = np.loadtxt(path + "SWIA.asc")

    density = swia[:, -1]

    t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    density_cut = density[inicio:fin]

    return swia, t_cut, density_cut


# ############################################################################################## LPW


def importar_swia_vel(year, month, day, ti, tf):
    # path = f'../../../datos/clweb/{year}-{month}-{day}/'
    path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    swia = np.loadtxt(path + "SWIA_vel.asc")

    density = swia[:, 6]
    vel_mso_xyz = swia[:, 7:9]

    t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    density_cut = density[inicio:fin]
    vel_cut = vel_mso_xyz[inicio:fin]

    return swia, t_cut, density_cut, vel_cut


# ############################################################################################## LPW


def importar_lpw(year, month, day, ti, tf):
    # path = f'../../../datos/clweb/{year}-{month}-{day}/'
    path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    lpw = np.loadtxt(path + "LPW.asc")

    e_density = lpw[:, -1]

    t = lpw[:, 3] + lpw[:, 4] / 60 + lpw[:, 5] / 3600

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    e_density_cut = e_density[inicio:fin]

    return lpw, t_cut, e_density_cut
