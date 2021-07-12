import numpy as np
import matplotlib.pyplot as plt
import sys

from mpl_toolkits.mplot3d import Axes3D  # de acá importo la proyección 3D

sys.path.append("..")

from funciones_plot import set_axes_equal
from funciones import SZA, donde, angulo


"""
Voy a agarrar un archivo solo y sobre ese calcular el SZA de MPB y BS
"""

year = 2016
# path = "../../../datos/MAG_1s/2016/mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts"
path = "../../../datos/MAG_1s/prueba/mvn_mag_l2_2017034ss1s_20170203_v01_r01.sts"

# Ajuste de Vignes de la MPB:
x0 = 0.78
e = 0.9
L = 0.96

theta = np.linspace(0, 3 * np.pi / 4, 100)
phi = np.linspace(0, np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(theta))

X = x0 + r * np.cos(THETA)  # en RM
Y = r * np.sin(THETA) * np.cos(PHI)
Z = r * np.sin(THETA) * np.sin(PHI)

R = np.transpose([X.flatten(), Y.flatten(), Z.flatten()])  # en RM


# importo los datos

mag = np.loadtxt(path, skiprows=160)

calendario_2017 = np.zeros(3)
calendario_2017[0] = mag[1, 1]  # doy

posicion = mag[:, 11:14]  # en km
B = mag[:, 7:10]
hh = mag[:, 2]

# vamos a trabajar sobre la primera órbita nomás
# creo que en general alcanzaría con ver la primera órbita del día, no voy a
# perder mucho y me ahorro tiempo.
# voy a descartar todos los puntos x,z < 0 para no tener que buscar ahí.
orbita = posicion / 3390  # es la posicion en RM
una_vuelta = int(len(orbita) / 5)

orb = orbita[:una_vuelta]

x_pos = np.array([orb[i, :] for i in range(len(orb)) if orb[i, 0] > 0])
posicion_RM = np.array([x_pos[i, :] for i in range(len(x_pos)) if x_pos[i, 2] > 0])

# hasta acá simplemente recorté la órbita

"""
Hace una lista con los índices de la órbita donde se parezca a R dentro de una
tolerancia rtol=0.01 (esto se puede modificar).
Entonces x, y, z me devuelven listas de bool de si se parece algún punto sobre
esa coordenada. Luego, xx, yy, zz me dicen el índice de esta nueva lista (x, y, z)
donde se parecen. Ojo que no es el índice de posicion_RM!!
Finalmente, posicion_RM[lista, :] me va a dar los valores de la órbita. Ahí uso
la función SZA y listo.
"""

lista = []
for i in range(len(posicion_RM)):
    x = np.isclose(posicion_RM[i, 0], R[:, 0], rtol=0.01)
    y = np.isclose(posicion_RM[i, 1], R[:, 1], rtol=0.01)
    z = np.isclose(posicion_RM[i, 2], R[:, 2], rtol=0.01)

    xx = np.where(x)[0]
    yy = np.where(y)[0]
    zz = np.where(z)[0]

    # si hay alguna coincidencia
    if any(np.isin(xx, yy)) and any(np.isin(zz, yy)) and any(np.isin(xx, zz)):
        lista.append(i)

if lista == []:
    print("ninguno cumple, en el for debería pasar al siguiente caso")
# me da una lista de las posiciones donde R se parece a la orbita: R[lista, :]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection="3d")
ax.set_xlabel(r"$X_{MSO} (R_m)$")
ax.set_ylabel(r"$Y_{MSO} (R_m)$")
ax.set_zlabel(r"$Z_{MSO} (R_m)$")
ax.plot_surface(
    X,
    Y,
    Z,
    rstride=4,
    cstride=4,
    alpha=0.5,
    edgecolor="none",
    cmap=plt.get_cmap("Blues_r"),
)
ax.plot(
    posicion_RM[:, 0],
    posicion_RM[:, 1],
    posicion_RM[:, 2],
    color="green",
    label="Órbita",
)
ax.scatter(posicion_RM[lista, 0], posicion_RM[lista, 1], posicion_RM[lista, 2])
plt.show()


"""
Clasificación por SZA
"""

SZA_MPB = SZA(posicion_RM, lista[0])  # tengo que elegir un solo elemento de lista
# me da mal si pongo la lista entera

# Ajuste de Vignes del BS:

x0 = 0.64
e = 1.03
L = 2.04

theta = np.linspace(0, 3 * np.pi / 4, 100)
phi = np.linspace(0, np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(theta))

X = x0 + r * np.cos(THETA)  # en RM
Y = r * np.sin(THETA) * np.cos(PHI)
Z = r * np.sin(THETA) * np.sin(PHI)

R = np.transpose([X.flatten(), Y.flatten(), Z.flatten()])  # en RM

lista_BS = []
if SZA_MPB < 90:
    for i in range(len(posicion_RM)):
        x = np.isclose(posicion_RM[i, 0], R[:, 0], rtol=0.01)
        y = np.isclose(posicion_RM[i, 1], R[:, 1], rtol=0.01)
        z = np.isclose(posicion_RM[i, 2], R[:, 2], rtol=0.01)

        xx = np.where(x)[0]
        yy = np.where(y)[0]
        zz = np.where(z)[0]

        # si hay alguna coincidencia
        if any(np.isin(xx, yy)) and any(np.isin(zz, yy)) and any(np.isin(xx, zz)):
            lista_BS.append(i)

    SZA_BS = SZA(posicion_RM, lista_BS[0])

    print(SZA_MPB, SZA_BS)
