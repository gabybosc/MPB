import numpy as np
import matplotlib.pyplot as plt
from plot_orbitas import datos_fijos
from funciones import find_nearest, UTC_to_hdec


def MPB(x0=0.78, e=0.9, L=0.96):
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    r1 = L / (1 + e * np.cos(THETA))

    # ahora plotea
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{MSO} (R_m)$")
    ax.set_ylabel(r"$Y_{MSO} (R_m)$")
    ax.set_zlabel(r"$Z_{MSO} (R_m)$")
    ax.set_aspect("equal")
    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)
    plot = ax.plot_surface(
        X1,
        Y1,
        Z1,
        rstride=4,
        cstride=4,
        alpha=0.5,
        edgecolor="none",
        cmap=plt.get_cmap("Blues_r"),
    )

    asc = L0 / (1 - e ** 2)  # semieje mayor
    bsc = np.sqrt(asc * L0)
    csc = e * asc - x0  # donde está centrada. Hay que ver el signo

    u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
    ax.plot_wireframe(
        np.cos(u) * np.sin(v),
        np.sin(u) * np.sin(v),
        np.cos(v),
        color="r",
        linewidth=0.5,
    )
    ax.legend()
    set_axes_equal(ax)

    # plt.savefig(f'../outputs/figs_MPB/ajuste_{fecha}.png')
    #
    return (asc, bsc, csc)


def normal_vignes(R, asc, bsc, csc):
    fig = plt.gcf()
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
        color="b",
        length=0.5,
        label="Normal del ajuste",
    )  # asi se plotea un vector

    return ()


def normal_MVA(R, x3):
    fig = plt.gcf()
    ax = plt.gca()
    ax.quiver(
        R[0],
        R[1],
        R[2],
        x3[0],
        x3[1],
        x3[2],
        color="k",
        length=0.5,
        label="Normal del MVA",
    )
    return ()


asc, bsc, csc = MPB()

t, posicion_cut, year, month, day = datos_fijos(2015, 10, 12, 18.75, 19.75)
t_medio = UTC_to_hdec("19:19:13")
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index, :]  # la posicion de la nave en RM
