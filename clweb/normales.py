import numpy as np
import matplotlib.pyplot as plt
from plot_orbitas import datos_fijos
import sys
from mpl_toolkits.mplot3d import axes3d  # lo usa aunque esté marcado como que no!

sys.path.append("..")
from funciones import find_nearest, UTC_to_hdec
from funciones_plot import set_axes_equal


"""
Plotea la fig 5 del poster de la AGU2019
"""


def MPB(x0=0.78, e=0.9, L=0.96):
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    r1 = L / (1 + e * np.cos(THETA))

    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)

    # ahora plotea
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{MSO} (R_m)$", fontsize=14)
    ax.set_ylabel(r"$Y_{MSO} (R_m)$", fontsize=14)
    # ax.set_zlabel(r"$Z_{MSO} (R_m)$", fontsize=14)
    ax.plot_surface(
        X1,
        Y1,
        Z1,
        rstride=4,
        cstride=4,
        alpha=0.5,
        edgecolor="gray",
        cmap=plt.get_cmap("Blues_r"),
    )

    u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
    ax.plot_wireframe(
        np.cos(u) * np.sin(v),
        np.sin(u) * np.sin(v),
        np.cos(v),
        color="#c1440e",
        linewidth=0.5,
    )
    plt.show()
    # set_axes_equal(ax)

    # plt.savefig(f'../outputs/figs_MPB/ajuste_{fecha}.png')
    #


def parametros_elipse(R, x0=0.78, e=0.9, L=0.96):  # R es la posicion del cruce en RM
    # conica que toma los parámetros de Vignes
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    # Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0, 0, 0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))

    asc = L0 / (1 - e ** 2)  # semieje mayor
    bsc = np.sqrt(asc * L0)
    csc = e * asc - x0  # donde está centrada. Hay que ver el signo

    return (asc, bsc, csc)


def normal_vignes(R, asc, bsc, csc):
    plt.gcf()
    ax = plt.gca()
    norm_vignes = np.array(
        [(R[0] + csc) * 2 / asc ** 2, R[1] * 2 / (bsc) ** 2, R[2] * 2 / (bsc) ** 2]
    )  # la normal de vignes
    norm_vignes = norm_vignes / np.linalg.norm(norm_vignes)  # normalizado
    ax.quiver(
        R[0],
        R[1],
        R[2],
        norm_vignes[0],
        norm_vignes[1],
        norm_vignes[2],
        color="r",
        length=0.5,
    )


def normal_MVA(R, x3, colour="k", lab="MVA normal"):
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

plt.show()
