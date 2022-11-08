import numpy as np
import glob as glob
from datetime import datetime
import datetime as dt
import matplotlib.pyplot as plt
from importar_datos import importar_MAG_pds

import sys

sys.path.append("..")

# from funciones import find_nearest
plt.ion()

"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos
en los cuales se cumplen:
SZA < 80º
B en el SW no es mayormente Bx

Se fija dónde es que coincide la posicion de MAVEN con el fit de vignes y mira
estas condiciones en ese punto.


"""

year = int(input("Year\n"))
# path = glob.glob('../../../MAVEN/mag_1s/{}/*/*.sts'.format(year)) #en desktop
path = glob.glob(f"../../../datos/VEX/{year}/*.tab")  # en laptop
cantidad_datos = len(path)
calendario = np.zeros(
    (5 * cantidad_datos, 6)
)  # la primera columna es el día del año, la segunda es el número de orbita,
# la tercera dice si cumple el SZA, la cuarta la altitud y la quinta el ZMSO


# Ajuste de Xu2021:
def altitude(SZA):
    alt = 0.11 * SZA ** 2 - 0.22 * SZA + 389
    return alt / 6050


def Xu2021():
    sza = np.linspace(0, np.pi, 100)
    alt = 1 + altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yy = y_alt[x_alt >= 0]  # me quedo solo con el dayside
    xx = x_alt[x_alt >= 0]
    return xx, yy


doy = 213
ti = 0
tf = 11
t, B, pos = importar_MAG_pds(year, doy, ti, tf)

for i, j in enumerate(path):
    mag = np.loadtxt(j, skiprows=160)
    # arr[i,:] = np.array([mag[1,1],len(mag)])
    posicion = np.zeros((len(mag[:, 0]), 3))
    for k in range(11, 14):
        posicion[:, k - 11] = mag[:, k]
    orbita = posicion / 3390
    una_vuelta = int(len(orbita) / 5)

    B = np.zeros((len(mag[:, 0]), 3))
    for k in range(7, 10):
        B[:, k - 7] = mag[:, k]

    """
    Son dos loops: el loop en i barre toda la superficie y la resta para cada
    punto de la órbita. El loop en j agarra esa resta y ve dónde es que es mínima
    (busca el máximo acercamiento entre la órbita y la superficie).
    Luego, guarda el mínimo para cada punto de la órbita. Finalmente,
    busca el mínimo de mínimos. Hace esto cada 500 puntos y sólo donde Z y X MSO
    son positivas, total es donde está mi cruce. (esto además me disminuye los
    falsos positivos)

    """
    orbitas = [
        orbita[:una_vuelta],
        orbita[una_vuelta : una_vuelta * 2],
        orbita[una_vuelta * 2 : una_vuelta * 3],
        orbita[una_vuelta * 3 : una_vuelta * 4],
        orbita[una_vuelta * 4 :],
    ]
    resta = np.zeros((len(R), 3))
    paso = 50  # se traduce en 150km si se mueve a 3km/s

    for indice, l in enumerate(orbitas):
        calendario[i * 5 + indice, 0] = mag[1, 1]
        calendario[i * 5 + indice, 1] = indice + 1
        pos = l * 3390
        X_MSO = pos[:, 0]
        Z_MSO = pos[:, 2]
        idx_min = np.zeros(int(una_vuelta / paso))
        max_acercamiento = np.zeros(int(una_vuelta / paso))
        minimo = 0
        for k in range(0, int(una_vuelta) - 100, paso):
            if Z_MSO[k] > 0 and X_MSO[k] > 0:
                for m in range(len(R)):
                    resta[m, :] = l[k, :] - R[m, :]
                A = np.linalg.norm(resta, axis=1)
                idx_min[int(k / paso)] = np.argmin(A)
                max_acercamiento[int(k / paso)] = A[int(idx_min[int(k / paso)])]
        if (
            sum(max_acercamiento) == 0
        ):  # si es cero, va a fallar todo el script, así que digo que esa órbita es mala y listo
            calendario[i * 5 + indice, 2] = 0
            calendario[i * 5 + indice, 3] = 0
            calendario[i * 5 + indice, 4] = 0
        else:
            minimo = np.where(
                max_acercamiento
                == np.min(max_acercamiento[np.nonzero(max_acercamiento)])
            )[0][
                0
            ]  # busca el minimo que no sea cero

        idx = minimo * paso
        """
        Clasificación por altitud (debería cumplirse siempre, pero aparentemente no pasa)
        """
        altitud = np.linalg.norm(pos[int(idx), :]) - 3390

        if altitud < 1300 and altitud > 300:
            calendario[i * 5 + indice, 2] = 1

        """
        Clasificación por SZA despues poner entre 30 y 60 y luego 60 y 90
        """

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
        # print(j, indice, SZA)
        if SZA < 30:
            calendario[i * 5 + indice, 3] = 1
        elif SZA > 60:
            calendario[i * 5 + indice, 4] = 1
        else:
            calendario[i * 5 + indice, 5] = 1

        """
        Clasificación por Z_MSO (ya se cumple porque lo pedí en el if)
        """
        Z_MSO = pos[int(idx), 2]
        if Z_MSO > 0:
            calendario[i * 5 + indice, 4] = 1


np.savetxt(
    f"../outputs/tiempos_{year}.txt",
    sorted(calendario, key=lambda x: x[0]),
    fmt="%10d",
    header="Día        Órbita        Altitud        SZA < 30     SZA entre 30 y 60    SZA > 60",
    newline="\r\n",
)

clasific_30 = calendario[:, 2] + calendario[:, 3]  # si suma 2, tiene SZA menor a 30.
clasific_60 = calendario[:, 2] + calendario[:, 4]  # si suma 2, tiene SZA entre 30 y 60.
clasific_90 = calendario[:, 2] + calendario[:, 5]  # si suma 2, tiene SZA mayor a 90.

fechas_SZA30 = np.zeros((len(clasific_30), 2))
fechas_SZA60 = np.zeros((len(clasific_60), 2))
fechas_SZA90 = np.zeros((len(clasific_90), 2))

for i in range(len(calendario)):
    if clasific_30[i] == 2:
        fechas_SZA30[i, 0] = calendario[i, 0]
        fechas_SZA30[i, 1] = calendario[i, 1]
    else:  # está adentro de un else para que sólo chequee el if cuando no es menor a 30
        if clasific_60[i] == 2:
            fechas_SZA60[i, 0] = calendario[i, 0]
            fechas_SZA60[i, 1] = calendario[i, 1]
        else:
            if clasific_90[i] == 2:
                fechas_SZA90[i, 0] = calendario[i, 0]
                fechas_SZA90[i, 1] = calendario[i, 1]

# para que me de un array ordenado y sin ceros:
fechas_SZA30_trim = fechas_SZA30[~np.all(fechas_SZA30 == 0, axis=1)]
fechas_SZA60_trim = fechas_SZA60[~np.all(fechas_SZA60 == 0, axis=1)]
fechas_SZA90_trim = fechas_SZA90[~np.all(fechas_SZA90 == 0, axis=1)]

np.savetxt(
    "SZA30_{}.txt".format(year),
    sorted(fechas_SZA30_trim, key=lambda x: x[0]),
    fmt="%10d",
    header="Las fechas de {} en las cuales tenemos SZA < 30º".format(year),
    newline="\r\n",
)

np.savetxt(
    "SZA60_{}.txt".format(year),
    sorted(fechas_SZA60_trim, key=lambda x: x[0]),
    fmt="%10d",
    header="Las fechas de {} en las cuales tenemos SZA > 30º y < 60º".format(year),
    newline="\r\n",
)

np.savetxt(
    "SZA90_{}.txt".format(year),
    sorted(fechas_SZA90_trim, key=lambda x: x[0]),
    fmt="%10d",
    header="Las fechas de {} en las cuales tenemos SZA > 60º".format(year),
    newline="\r\n",
)

doy = np.loadtxt("SZA30_{}.txt".format(year), skiprows=1)[:, 0]
month = np.zeros(len(doy))
for d in range(len(doy)):
    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(
        doy[d] - 1
    )  # para convertir el doy en date
    month[d] = date_orbit.strftime("%m")

plt.hist(month, 12, range=(1, 13))
plt.xlim(left=1, right=13)
plt.xlabel("Mes")
plt.ylabel("Cantidad de órbitas")
plt.title("Cantidad mensual de órbitas con SZA < 30º en {}".format(year))
