import numpy as np
import matplotlib.pyplot as plt
from funciones_plot import set_axes_equal
from funciones import Mij, error
from scipy.stats import norm

# from mpl_toolkits.mplot3d import Axes3D #de acá importo la proyección 3D

"""
Acá van las funciones que se relacionan con el MVA, el fit o el bootstrap
"""


def ajuste_conico(posicion, index, orbita, x3, x0=0.78, e=0.9, L=0.96):
    """conica que toma los parámetros de Vignes y devuelve la normal
    para plotear usar normales.py"""
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    R = posicion[index, :] / 3390  # la posicion de la nave en RM
    # ###### Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0, 0, 0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
    r1 = L0 / (1 + e * np.cos(THETA))

    # ahora plotea
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{MSO} (R_m)$")
    ax.set_ylabel(r"$Y_{MSO} (R_m)$")
    ax.set_zlabel(r"$Z_{MSO} (R_m)$")
    ax.set_aspect("equal")
    ax.plot(orbita[:, 0], orbita[:, 1], orbita[:, 2], color="green", label="Órbita")
    ax.scatter(R[0], R[1], R[2], label="MAVEN", color="k", s=40)
    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)
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

    asc = L0 / (1 - e ** 2)  # semieje mayor
    bsc = np.sqrt(asc * L0)
    csc = e * asc - x0  # donde está centrada. Hay que ver el signo

    norm_vignes = np.array(
        [(R[0] + csc) * 2 / asc ** 2, R[1] * 2 / bsc ** 2, R[2] * 2 / bsc ** 2]
    )
    norm_vignes = norm_vignes / np.linalg.norm(norm_vignes)

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
    )
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

    return norm_vignes, X1, Y1, Z1, R, L0


def bootstrap(N, B_cut, M_cut):
    """Hace el bootstrap para M puntos de un campo B, N veces."""
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
        phi, delta_B3 = error(lamb_ran, B_random, M_cut, x_ran)
        out[a] = np.mean(B3_ran)
        out_phi[a, 0] = phi[2, 0]
        out_phi[a, 1] = phi[2, 1]

    normal = np.mean(normal_ran, axis=0)
    normal = normal / np.linalg.norm(normal)

    return normal, phi, delta_B3, out, out_phi


def normal_fit(posicion, index, x0=0.78, e=0.9):
    """conica que toma los parámetros de Vignes y devuelve la normal
    para plotear usar ajuste_conico """
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    R = posicion[index, :] / 3390  # la posicion de la nave en RM
    # ###### Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0, 0, 0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))

    asc = L0 / (1 - e ** 2)  # semieje mayor
    bsc = np.sqrt(asc * L0)
    csc = e * asc - x0  # donde está centrada. Hay que ver el signo

    norm_vignes = np.array(
        [(R[0] + csc) * 2 / asc ** 2, R[1] * 2 / bsc ** 2, R[2] * 2 / bsc ** 2]
    )  # la normal de vignes
    norm_vignes = norm_vignes / np.linalg.norm(norm_vignes)  # normalizado

    return norm_vignes, L0


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

    ax1.quiver(
        R[0],
        R[1],
        R[2],
        norm_vignes[0],
        norm_vignes[1],
        norm_vignes[2],
        color="g",
        length=0.5,
        label="fit normal",
    )
    ax1.quiver(
        R[0],
        R[1],
        R[2],
        v_media[0],
        v_media[1],
        v_media[2],
        color="b",
        length=0.5,
        label="velocity",
    )
    ax1.quiver(
        R[0], R[1], R[2], x3[0], x3[1], x3[2], color="k", length=0.5, label="MVA normal"
    )
    ax1.quiver(
        R[0],
        R[1],
        R[2],
        v_para[0],
        v_para[1],
        v_para[2],
        color="r",
        length=0.5,
        label="v parallel",
    )
    ax1.quiver(
        R[0],
        R[1],
        R[2],
        v_para_MVA[0],
        v_para_MVA[1],
        v_para_MVA[2],
        color="m",
        length=0.5,
        label="v parallel MVA",
    )
    ax1.legend()

    u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
    ax1.plot_wireframe(
        np.cos(u) * np.sin(v),
        np.sin(u) * np.sin(v),
        np.cos(v),
        color="r",
        linewidth=0.5,
    )

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

    u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
    ax2.plot_wireframe(
        np.cos(u) * np.sin(v),
        np.sin(u) * np.sin(v),
        np.cos(v),
        color="r",
        linewidth=0.5,
    )

    ax2.quiver(
        R[0],
        R[1],
        R[2],
        J_v[0],
        J_v[1],
        J_v[2],
        color="r",
        length=1e-3,
        label="Corriente en volumen",
    )
    ax2.quiver(
        R[0],
        R[1],
        R[2],
        B_upstream[0],
        B_upstream[1],
        B_upstream[2],
        color="b",
        length=1e-2,
        label="B upstream",
    )
    ax2.quiver(
        R[0],
        R[1],
        R[2],
        B_downstream[0],
        B_downstream[1],
        B_downstream[2],
        color="g",
        length=1e-2,
        label="B downstream",
    )
    ax2.quiver(
        R[0],
        R[1],
        R[2],
        fuerza_mva[0],
        fuerza_mva[1],
        fuerza_mva[2],
        color="m",
        length=2e13,
        label="Fuerza MVA",
    )
    ax2.quiver(
        R[0],
        R[1],
        R[2],
        x3[0],
        x3[1],
        x3[2],
        color="k",
        length=0.5,
        label="Normal del MVA",
        linewidths=0.5,
    )

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
