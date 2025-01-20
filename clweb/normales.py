import numpy as np
import matplotlib.pyplot as plt
from plot_orbitas import datos_fijos
import sys
from mpl_toolkits.mplot3d import axes3d  # lo usa aunque esté marcado como que no!

from mpl_toolkits.mplot3d import Axes3D

sys.path.append("..")
from funciones import find_nearest, UTC_to_hdec

"""
Plotea la fig 5 del poster de la AGU2019
"""


def MPB():
    # Datos para la superficie
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    x0 = 0.78
    e = 0.9
    L = 0.96

    r1 = L / (1 + e * np.cos(THETA))
    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)

    # Crear la figura y el eje 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Crear la superficie wireframe
    ax.plot_wireframe(X1, Y1, Z1, color="gray", alpha=0.5, linewidth=0.5)

    # Ajustar la escala de los ejes para que sean iguales
    max_range = (
        np.array([X1.max() - X1.min(), Y1.max() - Y1.min(), Z1.max() - Z1.min()]).max()
        / 2.0
    )
    mid_x = (X1.max() + X1.min()) * 0.5
    mid_y = (Y1.max() + Y1.min()) * 0.5
    mid_z = (Z1.max() + Z1.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Opcional: usar set_box_aspect para una proporción 1:1:1
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel(r"$X_{MSO} (R_M)$", fontsize=14)
    ax.set_ylabel(r"$Y_{MSO} (R_M)$", fontsize=14)
    ax.set_zlabel(r"$Z_{MSO} (R_M)$", fontsize=14)


def parametros_elipse(R, x0=0.78, e=0.9, L=0.96):  # R es la posicion del cruce en RM
    # conica que toma los parámetros de Vignes
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    # Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0, 0, 0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))

    asc = L0 / (1 - e**2)  # semieje mayor
    bsc = np.sqrt(asc * L0)
    csc = e * asc - x0  # donde está centrada. Hay que ver el signo

    return (asc, bsc, csc)


def normal_vignes(R, asc, bsc, csc, colour="#ffa600"):
    plt.gcf()
    ax = plt.gca()
    norm_vignes = np.array(
        [(R[0] + csc) * 2 / asc**2, R[1] * 2 / (bsc) ** 2, R[2] * 2 / (bsc) ** 2]
    )  # la normal de vignes
    norm_vignes = norm_vignes / np.linalg.norm(norm_vignes)  # normalizado
    ax.quiver(
        R[0],
        R[1],
        R[2],
        norm_vignes[0],
        norm_vignes[1],
        norm_vignes[2],
        color=colour,
        length=0.5,
        label="Ajuste",
    )


def normal_MVA(R, x3, colour="#de425b", lab="MVA"):
    plt.gcf()
    ax = plt.gca()
    ax.quiver(
        R[0], R[1], R[2], x3[0], x3[1], x3[2], color=colour, length=0.5, label=lab
    )
    return ()


MPB()

x3 = [
    np.array([0.956, 0.048, -0.290]),
    np.array([0.956, -0.286, 0.070]),
    np.array([0.815, -0.575, 0.076]),
    np.array([0.871, -0.476, -0.117]),
    np.array([0.9198, -0.3021, 0.2505]),
    np.array([0.981, -0.032, 0.193]),
]
dates = [
    [2015, 10, 12],
    [2015, 10, 10],
    [2016, "04", "05"],
    [2016, "03", 31],
    [2016, "03", 16],
    [2017, 11, 24],
]
times = [[18.75, 19.75], [12, 13], [5, 6], [12.5, 13.5], [17.75, 18.75], [11.75, 12.75]]
tmean = ["19:19:13", "12:40:47", "05:16:30", "13:04:38", "18:13:50", "12:15:27"]

for i in range(len(x3)):
    t, posicion_cut, year, month, day = datos_fijos(
        dates[i][0], dates[i][1], dates[i][2], times[i][0], times[i][1]
    )
    t_medio = UTC_to_hdec(tmean[i])
    index = np.where(t == find_nearest(t, t_medio))[0][0]
    R = posicion_cut[index, :]  # la posicion de la nave en RM
    asc, bsc, csc = parametros_elipse(R)

    # normal_MVA(R, x3, '2015-oct-12')
    normal_MVA(R, x3[i])
    normal_vignes(R, asc, bsc, csc)
    if i == 0:
        plt.legend()

plt.show()
