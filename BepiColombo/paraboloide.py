import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
from Titan.leer_datos import importar_bepi
import sys

sys.path.append("..")
from funciones import (
    donde,
    Bpara_Bperp,
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


def altitude(sza, a=0.11, b=-0.22, c=389):
    """El SZA en grados!!"""
    alt = a * sza ** 2 + b * sza + c
    # alt es la altitud, en realidad yo quiero que la función me devuelva la coord r medida desde el (0,0) en RV
    return 1 + alt / 6050


def c_parametro(posicion, pos_MPB):
    r = np.linalg.norm(posicion[pos_MPB, :]) - 6050  # para convertirla en altitud
    theta = SZA(posicion, pos_MPB)
    c = r - 0.11 * theta ** 2 + 0.22 * theta

    return c


def paraboloide(a=0.11, b=-0.22, c=389):
    theta = np.linspace(0, np.pi * 0.5, 100)  # es el SZA
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    alt = altitude(THETA * 180 / np.pi, a, b, c)

    y_alt = np.array([alt[i] * np.sin(THETA[i]) for i in range(len(alt))])
    x = np.array([alt[i] * np.cos(THETA[i]) for i in range(len(alt))])

    y = y_alt * np.cos(PHI)
    z = y_alt * np.sin(PHI)

    return x, y, z


def plot_3d(pos_RV, R, n, c):
    x, y, z = paraboloide(c=c)
    nmva = np.array([0.7467, -0.6497, 0.1428])

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{MSO} (R_m)$")
    ax.set_ylabel(r"$Y_{MSO} (R_m)$")
    ax.set_zlabel(r"$Z_{MSO} (R_m)$")
    ax.plot(pos_RV[:, 0], pos_RV[:, 1], pos_RV[:, 2])
    ax.quiver(R[0], R[1], R[2], n[0], n[1], n[2])
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

    plt.show()


# x, y, z = paraboloide()
# R = punto(0, 0)  # theta y phi se cuentan desde el [1, 0, 0]
# norm = normal(R)
# plot_3d(x, y, z, R, norm)


# def desparametrizar(a, b, c):
#     L = b**2 / a
#     e = np.sqrt(1 - b**2/a**2)
#     x0 = e * a - c

#     return L, e, x0


year, month, day = 2021, "08", 10
ti_MVA, tf_MVA = 13.911111, 13.922222
t1, t2, t3, t4 = [13.8974, 13.9078, 13.9283, 13.9469]
t, B, posicion = importar_bepi(t1 - 0.5, t4 + 0.5)
Bnorm = np.linalg.norm(B, axis=1)

Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], t1 - 0.2, t4 + 0.2)

# # val = multi_plot(t, tpara, t_els, B, Bnorm, Bpara, Bperp, energy, JE_total, 2)

inicio_MVA = donde(t, ti_MVA)
fin_MVA = donde(t, tf_MVA)
pos_MPB = int((inicio_MVA + fin_MVA) * 0.5)
R = posicion[inicio_MVA, :] / 6050
# n = normal(R)
nmva = np.array([0.7467, -0.6497, 0.1428])
c = c_parametro(posicion, inicio_MVA)
plot_3d(posicion / 6050, R, nmva, c)

# xx, yz = fit()
# pos_RV = pos / 6050
# orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)
# sza = SZA(pos, inicio_MVA)
# R = pos_RV[inicio_MVA]


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
