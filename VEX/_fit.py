import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
import datetime as dt

import sys

sys.path.append("..")
from funciones import SZA
from funciones_plot import equal_axes

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

"""
Fit 2D hardcodeado para el 28 de octubre 
"""


def altitude(sza, a=0.11, b=-0.22, c=389):
    """El SZA en grados!!"""
    alt = a * sza ** 2 + b * sza + c
    # alt es la altitud, en realidad yo quiero que la función me devuelva la coord r medida desde el (0,0) en RV
    return 1 + alt / 6050


def fit_Xu():
    sza = np.linspace(0, np.pi * 0.5, 100)
    alt = altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def calcular_normal(sza):
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


def c_parametro(posicion, pos_MPB):
    r = np.linalg.norm(posicion[pos_MPB, :]) - 6050  # para convertirla en altitud
    theta = SZA(posicion, pos_MPB)
    c = r - 0.11 * theta ** 2 + 0.22 * theta

    return c


def b_parametro(posicion, pos_MPB):
    r = np.linalg.norm(posicion[pos_MPB, :]) - 6050
    theta = SZA(posicion, pos_MPB)
    b = (r - 0.11 * theta ** 2 - 389) / theta

    return b


def a_parametro(posicion, pos_MPB):
    """theta tiene que estar en grados"""
    r = np.linalg.norm(posicion[pos_MPB, :]) - 6050
    theta = SZA(posicion, pos_MPB)
    a = (r + 0.22 * theta - 389) / theta ** 2
    return a


def fit_2d(a=0.11, b=-0.22, c=389):
    sza = np.linspace(0, np.pi / 2, 100)
    alt = altitude(sza * 180 / np.pi, a, b, c)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]
    return xx, yz


def fit_3d(c=389):
    theta = np.linspace(0, np.pi * 2 / 4, 100)  # es el SZA
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    alt = altitude(THETA * 180 / np.pi, c=c)

    y_alt = np.array([alt[i] * np.sin(THETA[i]) for i in range(len(alt))])
    x = np.array([alt[i] * np.cos(THETA[i]) for i in range(len(alt))])

    y = y_alt * np.cos(PHI)
    z = y_alt * np.sin(PHI)

    return x, y, z


def plot_2D(pos_RV, R, n, c):
    xx, yz = fit_2d(c=c)
    orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)
    nmva = [0.391, 0.856]  # [0.517,0.103,0.850]

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


def plot_3d(x, y, z, R, norm, nmva):
    # nmva = [0.517, 0.103, 0.850]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{VSO} (R_V)$")
    ax.set_ylabel(r"$Y_{VSO} (R_V)$")
    ax.set_zlabel(r"$Z_{VSO} (R_V)$")
    ax.plot_wireframe(x, y, z, color="gray", alpha=0.5, linewidth=0.5)
    ax.quiver(
        R[0],
        R[1],
        R[2],
        norm[0],
        norm[1],
        norm[2],
        color="#ffa600",
        length=0.5,
        label="Normal del fit",
    )
    ax.quiver(
        R[0],
        R[1],
        R[2],
        nmva[0],
        nmva[1],
        nmva[2],
        color="#de425b",
        length=0.5,
        label="Normal del mva",
    )
    u, v = np.mgrid[0: 2 * np.pi: 20j, 0: np.pi: 10j]
    ax.plot_wireframe(
        np.cos(u) * np.sin(v),
        np.sin(u) * np.sin(v),
        np.cos(v),
        color="#eecb8b",
        linewidth=0.5,
    )

    # Ajustar la escala de los ejes para que sean iguales
    max_range = (
            np.array([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]).max() / 2.0
    )
    mid_x = (x.max() + x.min()) * 0.5
    mid_y = (y.max() + y.min()) * 0.5
    mid_z = (z.max() + z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    equal_axes(ax, x, y, z)
    plt.legend()
    plt.show()


def rotacion(phi, norm2d):
    n3 = np.array([norm2d[0], norm2d[1], 0])
    mat = np.array(
        [[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]]
    )
    normal_3d = np.dot(mat, n3)

    return normal_3d


def hallar_phi(R):
    x, y, z = R[0], R[1], R[2]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arccos(x / r)
    phi = np.sign(z) * np.arccos(y / np.sqrt(y ** 2 + z ** 2))

    return r, theta, phi


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
