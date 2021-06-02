import numpy as np
from os import listdir
import glob as glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import scipy.signal as signal
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys

sys.path.append("..")

from funciones import find_nearest
from funciones_plot import set_axes_equal

plt.ion()

"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos en los cuales se cumplen:
SZA < 45º (o incluso <30º)
Altitud entre 300 y 1300 km
Z_MSO > 0 (si después se pueden volver a bajar los datos y hacer Z_pc > 0, mejor)

Se fija dónde es que coincide la posicion de MAVEN con el fit de vignes y mira estas condiciones en ese punto.

tenemos datos desde 10/2014 hasta 02/2018
"""


# mag = np.loadtxt('../../../MAVEN/mag_1s/2016/03/mvn_mag_l2_2016085ss1s_20160325_v01_r01.sts', skiprows=148) #en la compu del iafe
# mag = np.loadtxt('../../datos/MAG_1s/mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148) #en mi compu
mag = np.loadtxt(
    "../../../datos/MAG_1s/2016/mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts",
    skiprows=160,
)  # en mi compu

calendario_2016 = np.zeros(
    (5, 5)
)  # la primera columna es el día del año, la segunda es el número de orbita, la tercera dice si cumple el SZA, la

B = np.zeros((len(mag[:, 0]), 3))
for j in range(7, 10):
    B[:, j - 7] = mag[:, j]

    B_norm = np.linalg.norm(B, axis=1)

    posicion = np.zeros((len(mag[:, 0]), 3))
    for j in range(11, 14):
        posicion[:, j - 11] = mag[:, j]

        orbita = posicion / 3390  # radios marcianos

una_vuelta = int(len(orbita) / 5)

# Ajuste de Vignes:
x0 = 0.78
e = 0.9
L = 0.96

theta = np.linspace(0, 3 * np.pi / 4, 100)
phi = np.linspace(0, np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(theta))

X = x0 + r * np.cos(THETA)
Y = r * np.sin(THETA) * np.cos(PHI)
Z = r * np.sin(THETA) * np.sin(PHI)

R = np.transpose(np.array([X.flatten(), Y.flatten(), Z.flatten()]))
resta = np.zeros((len(R), 3))
"""
Son dos loops: el loop en i barre toda la superficie y la resta para cada punto de la órbita. El loop en j agarra esa resta y ve dónde es que es mínima (busca el máximo acercamiento entre la órbita y la superficie). Luego, guarda el mínimo para cada punto de la órbita. Finalmente, busca el mínimo de mínimos.
Hace esto cada 100 puntos porque si no tarda mucho.

Necesito que haga la resta sólo en Z positivo, porque si no, crea cruces ficticios en el sur.
"""
orbitas = [
    orbita[:una_vuelta],
    orbita[una_vuelta : una_vuelta * 2],
    orbita[una_vuelta * 2 : una_vuelta * 3],
    orbita[una_vuelta * 3 : una_vuelta * 4],
    orbita[una_vuelta * 4 :],
]
paso = 50
for n, l in enumerate(orbitas):
    calendario_2016[0, n] = mag[1, 1]
    calendario_2016[1, n] = n + 1
    pos = l * 3390
    X_MSO = pos[:, 0]
    Z_MSO = pos[:, 2]
    idx_min = np.zeros(int(una_vuelta / paso))
    max_acercamiento = np.zeros(int(una_vuelta / paso))
    minimo = 0
    for j in range(0, int(una_vuelta) - paso, paso):
        if Z_MSO[j] > 0 and X_MSO[j] > 0:
            for i in range(len(R)):
                resta[i, :] = l[j, :] - R[i, :]
            A = np.linalg.norm(resta, axis=1)
            idx_min[int(j / paso)] = np.argmin(A)
            max_acercamiento[int(j / paso)] = A[int(idx_min[int(j / paso)])]
    if sum(max_acercamiento) == 0:
        calendario_2016[2, n] = 0
        calendario_2016[3, n] = 0
        calendario_2016[4, n] = 0

    else:
        minimo = np.where(
            max_acercamiento == np.min(max_acercamiento[np.nonzero(max_acercamiento)])
        )[0][
            0
        ]  # busca el minimo que no sea cero
        # print(minimo, max_acercamiento[minimo])
        # print(Z_MSO[minimo * 100], X_MSO[minimo * 100])

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        ax.set_xlabel(r"$X_{MSO} (R_m)$")
        ax.set_ylabel(r"$Y_{MSO} (R_m)$")
        ax.set_zlabel(r"$Z_{MSO} (R_m)$")
        # ax.set_aspect("equal")
        ax.plot(l[:, 0], l[:, 1], l[:, 2], color="C2", label="Órbita")
        ax.scatter(l[minimo * paso, 0], l[minimo * paso, 1], l[minimo * paso, 2])
        plot = ax.plot_surface(
            X,
            Y,
            Z,
            rstride=4,
            cstride=4,
            alpha=0.5,
            edgecolor="none",
            cmap=plt.get_cmap("Blues_r"),
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
        set_axes_equal(ax)  # para que tenga forma de esfera la esfera
        plt.show(block=False)

        """
        Clasificación.
        """

        idx = minimo * paso
        SZA = (
            np.arccos(
                np.clip(
                    np.dot(pos[int(idx)] / np.linalg.norm(pos[int(idx)]), [1, 0, 0]),
                    -1.0,
                    1.0,
                )
            )
            * 180
            / np.pi
        )
        altitud = np.linalg.norm(pos[int(idx), :]) - 3390
        Z_MSO_cruce = pos[int(idx), 2]
        if SZA < 30:
            calendario_2016[2, n] = 1
        if altitud < 1300 and altitud > 300:
            calendario_2016[3, n] = 1
        if Z_MSO_cruce > 0:
            calendario_2016[4, n] = 1

np.savetxt(
    "prueba.txt",
    np.transpose(calendario_2016),
    fmt="%10d",
    header="Día        Órbita        SZA        altitud       Z_MSO",
    newline="\r\n",
)
