import numpy as np
import matplotlib.pyplot as plt
from funciones_plot import set_axes_equal
from funciones import Mij, error
from scipy.stats import norm

from mpl_toolkits.mplot3d import Axes3D  # de acá importo la proyección 3D

"""
Acá van las funciones que se relacionan con el MVA, el fit o el bootstrap
"""


def flechas(ax, inicio, fin, c="k", lth=0.5, lbl=None):
    ax.quiver(
        inicio[0],
        inicio[1],
        inicio[2],
        fin[0],
        fin[1],
        fin[2],
        color=c,
        length=lth,
        label=lbl,
    )


def planeta(ax):
    u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
    ax.plot_wireframe(
        np.cos(u) * np.sin(v),
        np.sin(u) * np.sin(v),
        np.cos(v),
        color="r",
        linewidth=0.5,
    )


def ajuste_conico(R, x0=0.78, e=0.9):
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    # ###### Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0, 0, 0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
    r1 = L0 / (1 + e * np.cos(THETA))

    asc = L0 / (1 - e**2)  # semieje mayor
    bsc = np.sqrt(asc * L0)
    csc = e * asc - x0  # donde está centrada. Hay que ver el signo

    norm_vignes = np.array(
        [(R[0] + csc) * 2 / asc**2, R[1] * 2 / bsc**2, R[2] * 2 / bsc**2]
    )
    norm_vignes = norm_vignes / np.linalg.norm(norm_vignes)

    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)

    return (X1, Y1, Z1, L0, norm_vignes)


def ajuste_conico_plot(posicion, index, orbita, x3, x0=0.78, e=0.9):
    """conica que toma los parámetros de Vignes y devuelve la normal
    para plotear usar normales.py"""

    R = posicion[index, :] / 3390  # la posicion de la nave en RM
    X1, Y1, Z1, L0, norm_vignes = ajuste_conico(R)
    # ahora plotea
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{MSO} (R_m)$")
    ax.set_ylabel(r"$Y_{MSO} (R_m)$")
    ax.set_zlabel(r"$Z_{MSO} (R_m)$")
    ax.plot(orbita[:, 0], orbita[:, 1], orbita[:, 2], color="green", label="Órbita")
    ax.scatter(R[0], R[1], R[2], label="MAVEN", color="k", s=40)
    ax.plot_surface(
        X1,
        Y1,
        Z1,
        rstride=4,
        cstride=4,
        alpha=0.5,
        edgecolor="none",
        cmap=plt.get_cmap("Blues_r"),
    )

    flechas(ax, R, norm_vignes, c="b", lbl="Normal del Ajuste")
    flechas(ax, R, x3, lbl="Normal del MVA")

    planeta(ax)
    ax.legend()
    set_axes_equal(ax)

    return norm_vignes, X1, Y1, Z1, R, L0


def bootstrap(N, B_cut):
    """Hace el bootstrap para M puntos de un campo B, N veces."""
    M_cut = len(B_cut)
    out = np.zeros(N)
    out_phi = np.zeros((N, 2))
    normal_ran = np.zeros((N, 3))
    for a in range(N):
        index = np.random.choice(
            B_cut.shape[0], M_cut, replace=True
        )  # elige M índices de B, puede repetir (replace = True)
        B_random = B_cut[index, :]  # me da un B hecho a partir de estos índices random

        Mij_random = Mij(B_random)

        [lamb_ran, x_ran] = np.linalg.eigh(Mij_random)  # uso eigh porque es simetrica
        idx_ran = lamb_ran.argsort()[::-1]
        lamb_ran = lamb_ran[idx_ran]
        x_ran = x_ran[:, idx_ran]
        # ojo que a veces me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
        x3_ran = x_ran[:, 2]
        if x3_ran[0] < 0:
            x3_ran = -x3_ran
        normal_ran[a, :] = x3_ran

        B3_ran = np.dot(B_random, x3_ran)
        phi, delta_B3 = error(lamb_ran, B_random, x_ran)
        out[a] = np.mean(B3_ran)
        out_phi[a, 0] = phi[2, 0]
        out_phi[a, 1] = phi[2, 1]

    normal = np.mean(normal_ran, axis=0)
    normal = normal / np.linalg.norm(normal)

    return normal, phi, delta_B3, out, out_phi


def plot_velocidades(X1, Y1, Z1, R, norm_vignes, x3, v_media, v_para, v_para_MVA):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1, 1, 1, projection="3d")
    ax1.set_xlabel("x mso")
    ax1.set_ylabel("y mso")
    ax1.set_zlabel("z mso")
    ax1.set_aspect("equal")
    ax1.plot_surface(
        X1, Y1, Z1, rstride=4, cstride=4, cmap=plt.get_cmap("Blues_r"), alpha=0.5
    )
    ax1.scatter(R[0], R[1], R[2])

    flechas(ax1, R, norm_vignes, c="b", lbl="Normal del Ajuste")
    flechas(ax1, R, v_media, c="g", lbl="Velocidad")
    flechas(ax1, R, x3, lbl="Normal del MVA")
    flechas(ax1, R, v_para, c="r", lbl="vel paralela")
    flechas(ax1, R, v_para_MVA, c="m", lbl="vel paralela MVA")

    ax1.legend()

    planeta(ax1)
    set_axes_equal(ax1)


def plot_FLorentz(X1, Y1, Z1, R, J_v, B_upstream, B_downstream, fuerza_mva, x3):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1, 1, 1, projection="3d")
    ax2.set_xlabel(r"$X_{MSO} (R_M)$")
    ax2.set_xlim(left=2, right=0)
    ax2.set_ylabel(r"$Y_{MSO} (R_M)$")
    ax2.set_zlabel(r"$Z_{MSO} (R_M)$")
    ax2.set_aspect("equal")
    ax2.scatter(R[0], R[1], R[2])
    ax2.plot_surface(
        X1, Y1, Z1, rstride=4, cstride=4, alpha=0.5, cmap=plt.get_cmap("Blues_r")
    )

    planeta(ax2)

    flechas(ax2, R, J_v, c="b", lth=1e-3, lbl="J vol")
    flechas(ax2, R, B_upstream, c="g", lth=1e-2, lbl="B upstream")
    flechas(ax2, R, x3, lbl="Normal del MVA")
    flechas(ax2, R, B_downstream, c="r", lth=1e-2, lbl="B downstream")
    flechas(ax2, R, fuerza_mva, c="m", lth=2e13, lbl="Fuerza MVA")

    set_axes_equal(ax2)
    ax2.legend(loc="upper right", bbox_to_anchor=(1.1, 1.05))


def plot_bootstrap(out, out_phi):
    """Plotea  el histograma del bootstrap y fitea la gaussiana.
    Devuelve los valores de la gaussiana"""
    plt.figure()
    plt.subplot(311)
    n, bins, patches = plt.hist(out, 50, density=1, alpha=0.5)
    (muB, sigmaB) = norm.fit(out)
    y = norm.pdf(bins, muB, sigmaB)
    plt.plot(bins, y)
    plt.xlabel(r"$\langle B_3 \rangle$ (nT)")

    plt.subplot(312)
    n, bins, patches = plt.hist(out_phi[:, 0] * 57.2958, 50, density=1, alpha=0.5)
    (mu31, sigma31) = norm.fit(out_phi[:, 0] * 57.2958)
    y = norm.pdf(bins, mu31, sigma31)
    plt.plot(bins, y)
    plt.xlabel(r"$\Delta \phi_{31}$ (º)")

    plt.subplot(313)
    n, bins, patches = plt.hist(out_phi[:, 1] * 57.2958, 50, density=1, alpha=0.5)
    (mu32, sigma32) = norm.fit(out_phi[:, 1] * 57.2958)
    y = norm.pdf(bins, mu32, sigma32)
    plt.plot(bins, y)
    plt.xlabel(r"$\Delta \phi_{32}$ (º)")
    plt.tight_layout()

    return muB, sigmaB, mu31, sigma31, mu32, sigma32
