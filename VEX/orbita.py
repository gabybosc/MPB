import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as dt
from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
from _importar_datos import importar_MAG
import sys

sys.path.append("..")
from funciones import (
    donde,
    UTC_to_hdec,
    hdec_to_UTC,
    SZA,
)

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

"""
Ahora quiero que el fit sea considerando el punto por el que cruza
"""


def altitude(sza, a=0.11, b=-0.22, c=389):
    """El SZA en grados!!"""
    alt = a * sza ** 2 + b * sza + c
    # alt es la altitud, en realidad yo quiero que la función me devuelva la coord r medida desde el (0,0) en RV
    return 1 + alt / 6050


def normal(sza):
    """sza en rad"""
    r = 0.0597 * sza ** 2 - 0.002 * sza + 1.12
    dr = 0.1194 * sza - 0.002
    dx = dr * np.cos(sza) - r * np.sin(sza)
    dy = dr * np.sin(sza) + r * np.cos(sza)
    normal = np.array([-dy, dx])
    n = normal / np.linalg.norm(normal)

    if n[0] < 0:
        n = -n
    return n


def fit_Xu():
    sza = np.linspace(0, np.pi * 0.5, 100)
    alt = altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def fit_2d(a=0.11, b=-0.22, c=389):
    sza = np.linspace(0, np.pi / 2, 100)
    alt = altitude(sza * 180 / np.pi, a, b, c)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def plot_2D(pos_RV, idx):
    i1, i2, i3, i4 = idx[0], idx[1], idx[2], idx[3]
    # i5 = donde(tpos, 3)

    xx, yz = fit_2d()
    orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)

    fig, ax = plt.subplots()
    ax.plot(pos_RV[:, 0], orbita)
    ax.plot(xx, yz, color="#5647b4", linestyle="-.")
    ax.scatter(
        pos_RV[i1, 0], np.sqrt(pos_RV[i1, 1] ** 2 + pos_RV[i1, 2] ** 2), label="02:00"
    )
    ax.scatter(
        pos_RV[i2, 0], np.sqrt(pos_RV[i2, 1] ** 2 + pos_RV[i2, 2] ** 2), label="02:20"
    )
    ax.scatter(
        pos_RV[i3, 0], np.sqrt(pos_RV[i3, 1] ** 2 + pos_RV[i3, 2] ** 2), label="02:30"
    )
    ax.scatter(
        pos_RV[i4, 0], np.sqrt(pos_RV[i4, 1] ** 2 + pos_RV[i4, 2] ** 2), label="02:40"
    )
    # ax.scatter(pos_RV[i5, 0], np.sqrt(pos_RV[i5, 1] ** 2 + pos_RV[i5, 2] ** 2))

    ax.axis("equal")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
    ax.add_artist(circle)
    ax.set_title("VENUS VSO coordinates", fontsize=16)
    ax.set_xlabel(r"$X_{VSO}$ ($R_V$)", fontsize=14)
    ax.set_ylabel(r"$(Y²_{VSO} + Z²_{VSO} )^{1/2}$ ($R_V$)", fontsize=14)
    plt.legend()


# year, doy = 2011, 120  # fechas()
# year, doy = 2008, 301
year, doy = 2014, 116
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")

# t1, t2, t3, t4 = [2.755924008, 2.774626456, 2.785536217, 2.804238665]
# t1, t2, t3, t4 = (
#     8.541556527954604,
#     8.544405851015947,
#     8.551476393427427,
#     8.556111111111111,
# )

t1, t2, t3, t4 = [
    2.314023652,
    2.31986403,
    2.331544785,
    2.337131233,
]

t, B, posicion, cl, tpos = importar_MAG(year, doy, t1 - 0.5, t4 + 0.5)
Bnorm = np.linalg.norm(B, axis=1)

plot_2D(
    posicion / 6050,
    [donde(tpos, 2), donde(tpos, 2.33), donde(tpos, 2.5), donde(tpos, 2.66)],
)
plt.title(f"VEX {year}-{month}-{day}")
plt.show()
