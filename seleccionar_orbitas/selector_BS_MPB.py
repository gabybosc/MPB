import numpy as np
import glob as glob
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

from funciones import SZA, donde


"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos
en los cuales se cumplen que tanto el BS como la MPB tienen SZA < 90º en el
hemisferio norte.
Se fija dónde es que coincide la posicion de MAVEN con el fit de vignes y
mira estas condiciones en ese punto.

tenemos datos desde 10/2014 hasta 02/2018
"""

year = 2017
# path = glob.glob(f"../../../datos/MAG_1s/{year}/*.sts")
path = glob.glob("../../../datos/MAG_1s/prueba/*.sts")
cantidad_datos = len(path)
calendario_2017 = np.zeros(
    (3, cantidad_datos)
)  # la primera columna es el día del año, la segunda es la hora,
# la tercera vale 1 si el SZA de la MPB es < 90


# Ajuste de Vignes de la MPB:
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

R = np.transpose(np.array([X.flatten(), Y.flatten(), Z.flatten()]))  # en RM

# Importamos los datos y vemos uno por uno qué pasa.
for i, j in enumerate(
    path
):  # loop en todos los archivos .sts para cada año. i me da el índice en la lista, j me da el archivo
    print(f"{i * 100 / len(path):.3g}%\n")
    mag = np.loadtxt(j, skiprows=160)

    calendario_2017[0, i] = mag[1, 1]  # doy

    posicion = mag[:, 11:14]  # en km
    B = mag[:, 7:10]
    hh = mag[:, 2]
    orbita = posicion / 3390  # es la posicion en RM
    una_vuelta = int(len(orbita) / 5)

    """
    Son dos loops: el loop en i barre toda la superficie y la resta para cada
    punto de la órbita. El loop en j agarra esa resta y ve dónde es que es mínima
    (busca el máximo acercamiento entre la órbita y la superficie).
    Luego, guarda el mínimo para cada punto de la órbita.
    Finalmente, busca el mínimo de mínimos.
    Hace esto cada 100 puntos y sólo donde Z y X MSO son positivas,
    total es donde está mi cruce. (esto además me disminuye los falsos positivos)

    """

    # vamos a trabajar con la primera vuelta

    resta = [orbita[:una_vuelta:100] - R[i, :] for i in range(len(R))]

    a = np.linspace(0, len(orbita[:una_vuelta]), len(R))
    plt.plot(orbita[:una_vuelta, 0])
    plt.plot(a, R[:, 0], ".")
    plt.show()
    # a cada uno de los 173 datos de posición les resto todos los R
    # y veo cuál es el mínimo.
    # orbitas = [
    #     orbita[:una_vuelta],
    #     orbita[una_vuelta : una_vuelta * 2],
    #     orbita[una_vuelta * 2 : una_vuelta * 3],
    #     orbita[una_vuelta * 3 : una_vuelta * 4],
    #     orbita[una_vuelta * 4 :],
    # ]
    # resta = np.zeros((len(R), 3))

    for indice, l in enumerate(orbitas):
        pos_km = l * 3390
        X_MSO = pos_km[:, 0]
        Z_MSO = pos_km[:, 2]
        idx_min = np.zeros(int(una_vuelta / 100))
        max_acercamiento = np.zeros(int(una_vuelta / 100))
        minimo = 0
        for k in range(int(una_vuelta) - 100):
            if k % 100 == 0 and Z_MSO[k] > 0 and X_MSO[k] > 0:
                for m in range(len(R)):
                    resta[m, :] = l[k, :] - R[m, :]
                A = np.linalg.norm(resta, axis=1)
                idx_min[int(k / 100)] = np.argmin(A)
                max_acercamiento[int(k / 100)] = A[int(idx_min[int(k / 100)])]
        minimo = donde(
            max_acercamiento, np.min(max_acercamiento[np.nonzero(max_acercamiento)])
        )  # busca el minimo que no sea cero
        calendario_2017[1, i] = hh[minimo]

        """
        Clasificación por SZA
        """
        idx = minimo * 100

        SZA_MPB = SZA(posicion, idx)

        if SZA_MPB < 90:
            calendario_2017[2, i] = 1


# ahora, sobre los que cumplen esto para la MPB, busca que lo cumplan para el BS
# Ajuste de Vignes del BS:
x0 = 0.64
e = 1.03
L = 2.04

theta = np.linspace(0, 3 * np.pi / 4, 100)
phi = np.linspace(0, np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(theta))

X = x0 + r * np.cos(THETA)
Y = r * np.sin(THETA) * np.cos(PHI)
Z = r * np.sin(THETA) * np.sin(PHI)

R = np.transpose(np.array([X.flatten(), Y.flatten(), Z.flatten()]))
