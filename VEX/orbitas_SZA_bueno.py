import numpy as np
import glob as glob
from datetime import datetime
import datetime as dt
import matplotlib.pyplot as plt
from pathlib import Path


import sys

sys.path.append("..")

from funciones import donde, SZA

"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos
en los cuales se cumplen:
SZA < 80º
B en el SW no es mayormente Bx

Se fija dónde es que coincide la posicion de VEX con el fit de Xu2021 y mira
estas condiciones en ese punto.


"""

year = int(input("Year\n"))
# path = glob.glob('../../../MAVEN/mag_1s/{}/*/*.sts'.format(year)) #en desktop
path = glob.glob(f"../../../datos/VEX/{year}/*.tab")  # en laptop
cantidad_datos = len(path)
calendario = np.zeros((cantidad_datos, 4))
# la primera columna es el día del año, Las otras tres clasifican segun sza

# Ajuste de Xu2021:


def altitude(sza):
    alt = 0.11 * sza**2 - 0.22 * sza + 389
    return alt / 6050


def Xu2021_3D():
    p = 50
    sza = np.linspace(0, np.pi, p)
    alt = 1 + altitude(sza * 180 / np.pi)

    theta = np.linspace(0, 3 * np.pi / 4, p)
    phi = np.linspace(0, np.pi, p)
    THETA, PHI = np.meshgrid(theta, phi)
    X = alt * np.cos(THETA)
    Y = alt * np.sin(THETA) * np.cos(PHI)
    Z = alt * np.sin(THETA) * np.sin(PHI)

    R = np.transpose(np.array([X.flatten(), Y.flatten(), Z.flatten()]))
    return X, Y, Z, R


def Xu2021_2D():
    p = 50
    sza = np.linspace(0, np.pi, p)
    alt = 1 + altitude(sza * 180 / np.pi)

    yz_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = yz_alt[x_alt >= 0]  # me quedo solo con el dayside
    xx = x_alt[x_alt >= 0]
    return xx, yz


def importar_MAG(path):
    # chequea que el archivo no está vacío
    if Path(path).stat().st_size > 1000:
        B = np.genfromtxt(path, skip_header=1, usecols=[5, 6, 7])
        pos = np.genfromtxt(path, skip_header=1, usecols=[8, 9, 10])
        tt = np.genfromtxt(path, skip_header=1, usecols=0, dtype="str")

        fecha = np.array([x.split("T") for x in tt])
        hora = np.array([x.split(":") for x in fecha[:, 1]])
        hh = np.array([int(x) for x in hora[:, 0]])
        mm = np.array([int(x) for x in hora[:, 1]])
        ss = np.array([float(x) for x in hora[:, 2]])
        t = hh + mm / 60 + ss / 3600  # hdec
        doy = path.split("/")[-1].split("_")[-1].split(".")[0][4:]

        return t, B, pos, doy


for i, j in enumerate(path):
    print(j)
    t, B, pos, doy = importar_MAG(j)

    orbita = pos / 6050
    una_vuelta = len(orbita)

    """
    Son dos loops: el loop en i barre toda la superficie y la resta para cada
    punto de la órbita. El loop en j agarra esa resta y ve dónde es que es mínima
    (busca el máximo acercamiento entre la órbita y la superficie).
    Luego, guarda el mínimo para cada punto de la órbita. Finalmente,
    busca el mínimo de mínimos. Hace esto cada 500 puntos y sólo donde Z y X MSO
    son positivas, total es donde está mi cruce. (esto además me disminuye los
    falsos positivos)

    """

    X, Y, Z, R = Xu2021_3D()
    resta = np.zeros((len(R), 3))
    paso = 50  # se traduce en 150km si se mueve a 3km/s

    calendario[i, 0] = doy
    X_MSO = pos[:, 0]
    Z_MSO = pos[:, 2]
    idx_min = np.zeros(int(una_vuelta / paso))
    max_acercamiento = np.zeros(int(una_vuelta / paso))
    minimo = 0
    for k in range(0, int(una_vuelta) - 100, paso):
        if Z_MSO[k] > 0 and X_MSO[k] > 0:
            for m in range(len(R)):
                resta[m, :] = orbita[k, :] - R[m, :]
            A = np.linalg.norm(resta, axis=1)
            idx_min[int(k / paso)] = np.argmin(A)
            max_acercamiento[int(k / paso)] = A[int(idx_min[int(k / paso)])]
    if (
        sum(max_acercamiento) == 0
    ):  # si es cero, va a fallar todo el script, así que digo que esa órbita es mala y listo
        calendario[i, 1] = 0
        calendario[i, 2] = 0
        calendario[i, 3] = 0
    else:
        minimo = donde(
            max_acercamiento, np.min(max_acercamiento[np.nonzero(max_acercamiento)])
        )
        # busca el minimo que no sea cero

    idx = minimo * paso

    """
    Clasificación por SZA
    """

    sza = SZA(pos, idx)
    # print(j, indice, SZA)
    if sza < 70:
        calendario[i, 1] = 1
    elif sza > 80:
        calendario[i, 3] = 1
    else:
        calendario[i, 2] = 1


# # descomentar la siguiente sección si quiero ver el plot que muestre que funciona bien el script
# from funciones_plot import equal_axes
#
# from mpl_toolkits.mplot3d import Axes3D  # de acá importo la proyección 3D
#
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection="3d")
# ax.set_xlabel(r"$X_{MSO} (R_m)$", fontsize=14)
# ax.set_ylabel(r"$Y_{MSO} (R_m)$", fontsize=14)
# ax.set_zlabel(r"$Z_{MSO} (R_m)$", fontsize=14)
# ax.plot_surface(
#     X,
#     Y,
#     Z,
#     rstride=4,
#     cstride=4,
#     alpha=0.5,
#     # edgecolor="gray",
#     cmap=plt.get_cmap("Blues_r"),
# )
# ax.plot(pos[:, 0]/6050, pos[:, 1]/6050, pos[:, 2]/6050)
# ax.scatter(pos[idx, 0]/6050, pos[idx, 1]/6050, pos[idx, 2]/6050, label="VEX", color="k", s=40)


np.savetxt(
    f"../outputs/VEX_{year}.txt",
    sorted(calendario, key=lambda x: x[0]),
    fmt="%10d",
    header="Día\tSZA < 70\tSZA entre 70 y 80\tSZA > 80",
    newline="\r\n",
)

# clasific_30 = calendario[:, 2] + calendario[:, 3]  # si suma 2, tiene SZA menor a 30.
# clasific_60 = calendario[:, 2] + calendario[:, 4]  # si suma 2, tiene SZA entre 30 y 60.
# clasific_90 = calendario[:, 2] + calendario[:, 5]  # si suma 2, tiene SZA mayor a 90.
#
# fechas_SZA30 = np.zeros((len(clasific_30), 2))
# fechas_SZA60 = np.zeros((len(clasific_60), 2))
# fechas_SZA90 = np.zeros((len(clasific_90), 2))
#
# for i in range(len(calendario)):
#     if clasific_30[i] == 2:
#         fechas_SZA30[i, 0] = calendario[i, 0]
#         fechas_SZA30[i, 1] = calendario[i, 1]
#     else:  # está adentro de un else para que sólo chequee el if cuando no es menor a 30
#         if clasific_60[i] == 2:
#             fechas_SZA60[i, 0] = calendario[i, 0]
#             fechas_SZA60[i, 1] = calendario[i, 1]
#         else:
#             if clasific_90[i] == 2:
#                 fechas_SZA90[i, 0] = calendario[i, 0]
#                 fechas_SZA90[i, 1] = calendario[i, 1]
#
# # para que me de un array ordenado y sin ceros:
# fechas_SZA30_trim = fechas_SZA30[~np.all(fechas_SZA30 == 0, axis=1)]
# fechas_SZA60_trim = fechas_SZA60[~np.all(fechas_SZA60 == 0, axis=1)]
# fechas_SZA90_trim = fechas_SZA90[~np.all(fechas_SZA90 == 0, axis=1)]
#
# np.savetxt(
#     "SZA30_{}.txt".format(year),
#     sorted(fechas_SZA30_trim, key=lambda x: x[0]),
#     fmt="%10d",
#     header="Las fechas de {} en las cuales tenemos SZA < 30º".format(year),
#     newline="\r\n",
# )
#
# np.savetxt(
#     "SZA60_{}.txt".format(year),
#     sorted(fechas_SZA60_trim, key=lambda x: x[0]),
#     fmt="%10d",
#     header="Las fechas de {} en las cuales tenemos SZA > 30º y < 60º".format(year),
#     newline="\r\n",
# )
#
# np.savetxt(
#     "SZA90_{}.txt".format(year),
#     sorted(fechas_SZA90_trim, key=lambda x: x[0]),
#     fmt="%10d",
#     header="Las fechas de {} en las cuales tenemos SZA > 60º".format(year),
#     newline="\r\n",
# )
#
# doy = np.loadtxt("SZA30_{}.txt".format(year), skiprows=1)[:, 0]
# month = np.zeros(len(doy))
# for d in range(len(doy)):
#     date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(
#         doy[d] - 1
#     )  # para convertir el doy en date
#     month[d] = date_orbit.strftime("%m")
#
# plt.hist(month, 12, range=(1, 13))
# plt.xlim(left=1, right=13)
# plt.xlabel("Mes")
# plt.ylabel("Cantidad de órbitas")
# plt.title("Cantidad mensual de órbitas con SZA < 30º en {}".format(year))
