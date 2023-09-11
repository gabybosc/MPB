import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
from leer_datos import importar_bepi
import sys

sys.path.append("..")
from funciones import (
    donde,
    fechas,
    tiempos,
    Mij,
    error,
    corrientes,
    ancho_mpb,
    Bpara_Bperp,
    find_nearest,
    next_available_row,
    SZA,
)

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

"""
Tengo que hacer el fit este que es 2D en 3D para que tenga sentido calcular la normal
"""


def altitude(SZA):
    """El SZA tiene que estar en grados!!"""

    alt = 0.11 * SZA**2 - 0.22 * SZA + 389
    # alt es simplemente la altitud *desde la corteza* en km, por eso le sumo 1 RV
    return 1 + alt / 6050


def fit_2d():
    sza = np.linspace(0, np.pi / 2, 100)
    alt = altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def fit_3d():
    theta = np.linspace(0, np.pi * 2 / 4, 100)  # es el SZA
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    alt = altitude(THETA * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(THETA[i]) for i in range(len(alt))])
    x = np.array([alt[i] * np.cos(THETA[i]) for i in range(len(alt))])

    y = y_alt * np.cos(PHI)
    z = y_alt * np.sin(PHI)

    return x, y, z


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


def normal(p):
    """la normal de un paraboloide de revolución en un punto p en cartesianas"""
    if p[0] == 0:
        a = 1
    else:
        a = (p[1] ** 2 + p[2] ** 2) / p[0] ** 2
        # "a" es el parámetro pero cambia punto a punto
    print(p)
    norm = np.array([2 * p[0], 2 * p[1] / a**2, 2 * p[2] / a**2])
    norm = norm / np.linalg.norm(norm)
    return norm


def punto(theta, phi):
    sza = theta * 180 / np.pi
    alt = altitude(sza)

    x = alt * np.cos(theta)
    y = alt * np.sin(theta) * np.cos(phi)
    z = alt * np.sin(theta) * np.sin(phi)

    R = [x, y, z]
    return R


def plot_3d(x, y, z, R, norm):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{MSO} (R_m)$")
    ax.set_ylabel(r"$Y_{MSO} (R_m)$")
    ax.set_zlabel(r"$Z_{MSO} (R_m)$")

    u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
    ax.plot_wireframe(
        np.cos(u) * np.sin(v),
        np.sin(u) * np.sin(v),
        np.cos(v),
        color="#c1440e",
        linewidth=0.5,
    )
    ax.plot_surface(
        x,
        y,
        z,
        rstride=4,
        cstride=4,
        alpha=0.5,
        edgecolor="none",
        cmap=plt.get_cmap("Blues_r"),
    )
    ax.quiver(
        R[0],
        R[1],
        R[2],
        norm[0],
        norm[1],
        norm[2],
        color="k",
        # length=0.5,
        label="Normal del MVA",
    )
    plt.show()


x, y, z = fit_3d()
R = punto(0, 0)  # theta y phi se cuentan desde el [1, 0, 0]
norm = normal(R)
plot_3d(x, y, z, R, norm)


# def desparametrizar(a, b, c):
#     L = b**2 / a
#     e = np.sqrt(1 - b**2/a**2)
#     x0 = e * a - c

#     return L, e, x0


year, month, day = 2021, "08", 10  # fechas()
ti_MVA, tf_MVA = 13.911111, 13.922222
t1, t2, t3, t4 = [13.8974, 13.9078, 13.9283, 13.9469]
t, B, posicion = importar_bepi(t1 - 0.5, t4 + 0.5)
Bnorm = np.linalg.norm(B, axis=1)

Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], t1 - 0.2, t4 + 0.2)

# # val = multi_plot(t, tpara, t_els, B, Bnorm, Bpara, Bperp, energy, JE_total, 2)

inicio_MVA = donde(t, ti_MVA)
fin_MVA = donde(t, tf_MVA)


# xx, yz = fit()
# pos_RV = pos / 6050
# orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)
# sza = SZA(pos, inicio_MVA)
# R = pos_RV[inicio_MVA]


# def fit_R(R, sza):
#     a = (np.linalg.norm(R) - 1) * 6050 - 0.22 * sza - 389
#     sza_array = np.linspace(0, np.pi / 2, 100)
#     alt = 1 + (a + 0.22 * (sza_array * 180 / np.pi) + 389) / 6050

#     y_alt = np.array([alt[i] * np.sin(sza_array[i]) for i in range(len(alt))])
#     x_alt = np.array([alt[i] * np.cos(sza_array[i]) for i in range(len(alt))])

#     yz = y_alt[x_alt >= 0]
#     xx = x_alt[x_alt >= 0]
#     return xx, yz


# xx, yz = fit()
# xxR, yzR = fit_R(R, sza)
# pos_RV = pos / 6050
# orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)
# plot_orbita(pos_RV, orbita, xx, yz)
# plt.plot(xxR, yzR, color="m", linestyle="--")
# plt.show()


# # ahroa busco la normal

# # plt.scatter(x_alt, y_alt, c="m")
# # # plt.scatter(pos_RV[inicio_MVA, 0], np.sqrt(pos_RV[inicio_MVA, 1]**2 + pos_RV[inicio_MVA, 2]**2))
# # plt.scatter(R[0], np.sqrt(R[1] ** 2 + R[2] ** 2))

# p = 0.22 * (sza - 1)
# tg = np.array([1, p])

# plt.quiver(
#     pos_RV[inicio_MVA, 0],
#     np.sqrt(pos_RV[inicio_MVA, 1] ** 2 + pos_RV[inicio_MVA, 2] ** 2),
#     tg[0],
#     tg[1],
# )
# m = -1 / p
# y = m * ()
# k = altitude(sza) - 1 / p * sza
# n = -1 / p * sza + k
