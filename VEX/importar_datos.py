import numpy as np
import sys
from glob import glob

sys.path.append("..")
from funciones import donde

np.set_printoptions(precision=4)


def importar_VEX_mag_AMDA(year, month, day, ti, tf):
    path_B = glob(f"../../../datos/VEX/{day}{month}{year}/*.txt")
    path_pos = glob(f"../../../datos/VEX/{day}{month}{year}_pos/*.txt")
    B = np.genfromtxt(path_B[0], usecols=[1, 2, 3])
    pos = np.genfromtxt(path_pos[0], usecols=[1, 2, 3])
    tt = np.genfromtxt(path_B[0], usecols=0, dtype="str")
    fecha = np.array([x.split("T") for x in tt])
    hora = np.array([x.split(":") for x in fecha[:, 1]])
    hh = np.array([int(x) for x in hora[:, 0]])
    mm = np.array([int(x) for x in hora[:, 1]])
    ss = np.array([float(x) for x in hora[:, 2]])
    t = hh + mm / 60 + ss / 3600  # hdec
    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    return t_cut, B_cut, pos


def importar_MAG_pds(year, doy, ti, tf):
    path = f"../../../datos/VEX/{year}/VEX_MAG_{year}{doy}.tab"

    B = np.genfromtxt(path, skip_header=1, usecols=[5, 6, 7])
    pos = np.genfromtxt(path, skip_header=1, usecols=[8, 9, 10])
    tt = np.genfromtxt(path, skip_header=1, usecols=0, dtype="str")

    fecha = np.array([x.split("T") for x in tt])
    hora = np.array([x.split(":") for x in fecha[:, 1]])
    hh = np.array([int(x) for x in hora[:, 0]])
    mm = np.array([int(x) for x in hora[:, 1]])
    ss = np.array([float(x) for x in hora[:, 2]])
    t = hh + mm / 60 + ss / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    return t_cut, B_cut, pos
