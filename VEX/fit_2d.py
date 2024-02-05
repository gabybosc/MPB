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
    alt = a * sza**2 + b * sza + c
    # alt es la altitud, en realidad yo quiero que la función me devuelva la coord r medida desde el (0,0) en RV
    return 1 + alt / 6050


def normal(sza):
    """sza en rad"""
    r = 0.0597 * sza**2 - 0.002 * sza + 1.12
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


def c_parametro(posicion, pos_MPB):
    r = np.linalg.norm(posicion[pos_MPB, :]) - 6050  # para convertirla en altitud
    theta = SZA(posicion, pos_MPB)
    c = r - 0.11 * theta**2 + 0.22 * theta

    return c


def b_parametro(posicion, pos_MPB):
    r = np.linalg.norm(posicion[pos_MPB, :]) - 6050
    theta = SZA(posicion, pos_MPB)
    b = (r - 0.11 * theta**2 - 389) / theta

    return b


def a_parametro(posicion, pos_MPB):
    """theta tiene que estar en grados"""
    r = np.linalg.norm(posicion[pos_MPB, :]) - 6050
    theta = SZA(posicion, pos_MPB)
    a = (r + 0.22 * theta - 389) / theta**2
    return a


def fit_2d(a=0.11, b=-0.22, c=389):
    sza = np.linspace(0, np.pi / 2, 100)
    alt = altitude(sza * 180 / np.pi, a, b, c)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def plot_2D(pos_RV, R, n, c):
    xx, yz = fit_2d(c=c)
    orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)
    nmva = np.array([0.391, 0.92])

    fig, ax = plt.subplots()
    ax.plot(pos_RV[:, 0], orbita)
    ax.scatter(
        pos_RV[0, 0], orbita[0], s=50, zorder=2, marker="o", color="r", label="start"
    )
    ax.scatter(
        pos_RV[-1, 0], orbita[-1], s=50, zorder=2, marker="x", color="k", label="end"
    )
    ax.plot(xx, yz, color="#5647b4", linestyle="-.")
    ax.quiver(R[0], R[1], n[0], n[1], scale=10, label="fit")
    ax.quiver(R[0], R[1], nmva[0], nmva[1], color="C2", scale=10, label="MVA")
    ax.scatter(R[0], R[1])
    ax.axis("equal")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 2.5)
    circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
    ax.add_artist(circle)
    ax.set_title("VENUS VSO coordinates", fontsize=16)
    ax.set_xlabel(r"$X_{VSO}$ ($R_V$)", fontsize=14)
    ax.set_ylabel(r"$(Y²_{VSO} + Z²_{VSO} )^{1/2}$ ($R_V$)", fontsize=14)
    plt.legend()


year, doy = 2011, 120  # fechas()
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")

ti_MVA, tf_MVA = 2.778339444, 2.780409722
t1, t2, t3, t4 = [2.755924008, 2.774626456, 2.785536217, 2.804238665]
t, B, posicion, cl, tpos = importar_MAG(year, doy, t1 - 0.5, t4 + 0.5)
Bnorm = np.linalg.norm(B, axis=1)

inicio_MVA = donde(tpos, ti_MVA)
fin_MVA = donde(tpos, tf_MVA)
pos_MPB = int(0.5 * (fin_MVA + inicio_MVA))

R = posicion[pos_MPB, :]
R_2d = [R[0], np.sqrt(R[1] ** 2 + R[2] ** 2)]

alt = np.linalg.norm(posicion[pos_MPB, :]) - 1
sza_rad = SZA(posicion, pos_MPB) / 180 * np.pi
sza = SZA(posicion, pos_MPB)
n = normal(sza_rad)
# n = normal_cartesianas(1 + alt / 6050, sza_rad)
c = c_parametro(posicion, pos_MPB)

n_mva = [0.391, 0.92]  # [0.391,	-0.129, 0.911]

angulo_mva = np.arccos(np.clip(np.dot(n_mva, n), -1.0, 1.0))

print(
    f"El ángulo entre las normales 2D de MVA y del fit es {angulo_mva * 180 / np.pi:.3g}º"
)

plot_2D(posicion, R_2d, n, c)
plt.title(f"{year}-{month}-{day}")
plt.show()

"""
descomentar lo siguiente si quiero chequear qué pasa cambiando diferentes parámetros

a = a_parametro(posicion, pos_MPB)
b = b_parametro(posicion, pos_MPB)
R_mpr = posicion[mpr, :] / 6050
R_bs = posicion[bs, :] / 6050
R_mpr2d = [R_mpr[0], np.sqrt(R_mpr[1] ** 2 + R_mpr[2] ** 2)]
R_bs2d = [R_bs[0], np.sqrt(R_bs[1] ** 2 + R_bs[2] ** 2)]
xx, yz = fit_Xu()
XX_a, YZ_a = fit_2d(a=a)
XX_b, YZ_b = fit_2d(b=b)
XX_c, YZ_c = fit_2d(c=c)
orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)
fig, ax = plt.subplots()
ax.plot(pos_RV[:, 0], orbita)
ax.scatter(R_mpr2d[0], R_mpr2d[1], label="MPR (13:50)")
ax.scatter(R_2d[0], R_2d[1], label="MPB (13:55)")
ax.scatter(R_bs2d[0], R_bs2d[1], label="BS (14:00)")
ax.plot(xx, yz, color="#5647b4", linestyle="-.", label="Xu")
ax.plot(XX_c, YZ_c, color="C3", linestyle="-.", label="c")
ax.plot(XX_b, YZ_b, color="C4", linestyle="-.", label="b")
ax.plot(XX_a, YZ_a, color="C5", linestyle="-.", label="a")
ax.axis("equal")
ax.set_xlim(0, 2.5)
ax.set_ylim(0, 2.5)
circle = plt.Circle((0, 0), 1, color="#eecb8b", clip_on=True)
ax.add_artist(circle)
ax.set_title("VENUS VSO coordinates", fontsize=16)
ax.set_xlabel(r"$X_{VSO}$ ($R_V$)", fontsize=14)
ax.set_ylabel(r"$(Y²_{VSO} + Z²_{VSO} )^{1/2}$ ($R_V$)", fontsize=14)
plt.legend()
plt.show()
"""
