import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

import sys

sys.path.append("..")

from funciones import donde

path = "../../../datos/simulacion_leonardo/"
campo = np.loadtxt(path + "campo.txt")
posicion = np.loadtxt(path + "pos_mhd.txt")

"""
Para hacer un remesh voy a: poner x,y,z + B y ordenar por x creciente, y creciente,
z creciente (en ese orden).
Tengo que sacar los ceros porque son espacios sin datos.
"""

ceros = [i for i in range(len(posicion)) if posicion[i, 1] == 0]

pos_cut = np.delete(posicion, ceros, axis=0)
B_cut = np.delete(campo, ceros, axis=0)

val = np.concatenate((pos_cut, B_cut), axis=1)

X = np.array([val[i, :] for i in range(len(val)) if 0 < val[i, 0] < 3])
Y = np.array([X[i, :] for i in range(len(X)) if np.abs(X[i, 1]) <= 3])
Z = np.array([Y[i, :] for i in range(len(Y)) if np.abs(Y[i, 2]) <= 3])
# Z es el cut final y el que voy a usar el resto del tiempo

# vamos a reordenar.

datos = sorted(Z, key=itemgetter(0))
# si quiero que sortee por más de una columna, hago itemgetter(0,1)

# x_val = [x[0] for x in z]
# y_val = [x[1] for x in z]
#
# plt.scatter(zcut[:, 0], zcut[:, 1])
# plt.plot(x_val, y_val, c="C1")
# plt.show()

# hago un histograma para ver la separación entre datos.
# plt.hist(Z[:, 0], bins="auto")  # este muestra que efectivamente hay más cerca de cero
# plt.show()

"""
Voy a tomar un deltaX pequeño, según lo que vi en el histograma. Con ese valor
obtengo un paso y un tamaño para mi nuevo mesh. Luego, hago un loop donde voy
llenando un nuevo eje x más pequeño (llamado A) con los valores promediados
entre xmin y xmin+deltaX (es decir, entre idx_min y idx_min+paso).
"""

X = [x[0] for x in datos]  # es el eje x, debería estar ordenado.
xmin = min(X)
deltaX = 0.00004
minimo = donde(X, xmin)  # deberia ser X[0]
maximo = donde(X, xmin + deltaX)
pasos = len(X[minimo:maximo])

N = int(len(X) / pasos)  # cantidad de X
A = np.zeros(N)
for i in range(N):
    maximo = minimo + pasos
    A[i] = np.mean(X[minimo:maximo])
    minimo = maximo
