import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_MAG_pds

import sys

sys.path.append("..")
from funciones import datenum, donde, fechas, tiempos

np.set_printoptions(precision=4)


def altitude(SZA):
    alt = 0.11 * SZA ** 2 - 0.22 * SZA + 389
    return alt / 6050


def fit_Xu():
    """ Devuelve el fit de Xu 2021    """
    sza = np.linspace(0, np.pi, 100)
    alt = 1 + altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def fit_R(R, sza):
    """ Ajusta el fit de Xu 2021 para la posicion de la nave"""
    a = (np.linalg.norm(R) - 1) * 6050 - 0.22 * sza - 389
    sza_array = np.linspace(0, np.pi / 2, 100)
    alt = 1 + (a + 0.22 * (sza_array * 180 / np.pi) + 389) / 6050

    y_alt = np.array([alt[i] * np.sin(sza_array[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza_array[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def plot_orbita(pos_RV, orbita, xx, yz):
    fig, ax = plt.subplots()
    ax.plot(pos_RV[:, 0], orbita)
    ax.scatter(
        pos_RV[0, 0], orbita[0], s=50, zorder=2, marker="o", color="r", label="start"
    )
    ax.scatter(
        pos_RV[-1, 0], orbita[-1], s=50, zorder=2, marker="x", color="k", label="end"
    )
    ax.plot(xx, yz, color="#5647b4", linestyle="-.")
    ax.axis("equal")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 2.5)
    circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
    ax.add_artist(circle)
    ax.set_title("VENUS VSO coordinates", fontsize=16)
    ax.set_xlabel(r"$X_{VSO}$ ($R_V$)", fontsize=14)
    ax.set_ylabel(r"$(Y²_{VSO} + Z²_{VSO} )^{1/2}$ ($R_V$)", fontsize=14)
    plt.legend()
    # plt.show()


# year, month, day, doy = fechas()
# ti, tf = tiempos()
#
# t, B, pos = importar_MAG_pds(year, doy, ti, tf)
# # t, B, pos = importar_MAG_pds(2011, 213, 0, 10)
# pos_RV = pos / 6050
# xx, yz = fit()
# orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)
# plot_orbita(pos_RV, orbita, xx, yz)