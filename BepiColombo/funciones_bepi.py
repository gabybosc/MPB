import numpy as np
import sys

sys.path.append("..")
from funciones import datenum, donde, Bpara_Bperp


def tiempos_UTC(yy, mm, dd, t):
    tt = np.array([np.datetime64(datenum(yy, mm, dd, x)) for x in t])
    return tt


def altitude(SZA):
    alt = 0.11 * SZA**2 - 0.22 * SZA + 389
    return alt / 6050


def fit():
    sza = np.linspace(0, np.pi, 100)
    alt = 1 + altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]

    return xx, yz
