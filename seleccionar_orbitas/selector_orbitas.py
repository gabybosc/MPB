import numpy as np
import glob as glob
import sys

sys.path.append("..")
from funciones import SZA


"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos
en los cuales se cumplen:
SZA <30º
Z_MSO > 0
(si después se pueden volver a bajar los datos y hacer Z_pc > 0, mejor)

Se fija dónde es que coincide la posicion de MAVEN con el fit de vignes y
mira estas condiciones en ese punto.

Esta es una versión actualizada en 2021 cuando me tocó filtrar también BS.
"""

year = 2016
path = glob.glob(f"../../../datos/MAG_1s/{year}/*.sts")
# path = glob.glob("../../../datos/MAG_1s/prueba/*.sts")

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


calendario = np.zeros((len(path), 2))
for i, j in enumerate(path):

    mag = np.loadtxt(j, skiprows=160)
    calendario[i, 0] = int(mag[1, 1])  # doy

    posicion = mag[:, 11:14]  # en km
    B = mag[:, 7:10]
    hh = mag[:, 2]

    """Vamos a trabajar sobre la primera órbita nomás ya que creo que en general
    no voy a perder mucho y me ahorro tiempo.
    Voy a descartar todos los puntos x,z < 0 para no tener que buscar ahí total no
    me interesa ni el lado de noche ni el sur."""

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
    for k in range(len(posicion_RM)):
        x = np.isclose(posicion_RM[k, 0], R[:, 0], rtol=0.01)
        y = np.isclose(posicion_RM[k, 1], R[:, 1], rtol=0.01)
        z = np.isclose(posicion_RM[k, 2], R[:, 2], rtol=0.01)

        xx = np.where(x)[0]
        yy = np.where(y)[0]
        zz = np.where(z)[0]

        # si hay alguna coincidencia
        if any(np.isin(xx, yy)) and any(np.isin(zz, yy)) and any(np.isin(xx, zz)):
            lista.append(k)

        # me da una lista de las posiciones donde R se parece a la orbita: R[lista, :]
    if lista != []:  # si está vacía, va a pasar al caso siguiente.

        SZA_MPB = SZA(
            posicion_RM, lista[0]
        )  # tengo que elegir un solo elemento de lista
        # me da mal si pongo la lista entera
        if SZA_MPB < 30:
            calendario[i, 1] = 1
        else:
            calendario[i, 1] = 0
    else:
        calendario[i, 1] = 0


a = np.array(
    [calendario[i, 0] for i in range(len(calendario)) if calendario[i, 1] == 1]
)
a.sort()

np.savetxt(f"orbitas_SZA30_{year}.txt", a, delimiter="\t", fmt="%i")
