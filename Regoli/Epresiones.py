import numpy as np
from scipy.spatial import KDTree
from mpl_toolkits.mplot3d import Axes3D  # de acá importo la proyección 3D
import matplotlib.pyplot as plt
# import sys

# sys.path.append("..")
# from funciones import donde_kD

# fill the cube with random points
path = "../../../datos/simulacion_leonardo/"
data = np.load(path + "cubo_2RM.npy")

posicion = data[:, :3]
m = 3
n = len(posicion)

kdtree = KDTree(posicion)

# prueba con un punto cualquiera sobre la grilla
# un punto cualquiera, el más cercano al [1.5,1.5,1.5]
distancia, indice = kdtree.query(1.5 * np.ones((1,m)))
point = posicion[indice[0]]

# Coords of k-Nearest Neighbors
k = 100  # elijo 100 para después elegir un subgrupo que me guste sobre estos
dists, idxs = kdtree.query(point, k)
vecinos = data[idxs,:]
pos_vecinos = posicion[idxs,:]

# grafico los vecinos en 3D a ver si me da un cubito
# quiero ver los puntos sobre el eje x cuán lejos en y,z están
x = np.where(posicion[:,0] == point[0])[0]
posx = posicion[x]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection="3d")
ax.set_xlabel(r"$X_{MSO} (R_m)$")
ax.set_ylabel(r"$Y_{MSO} (R_m)$")
ax.set_zlabel(r"$Z_{MSO} (R_m)$")
ax.scatter(point[0], point[1], point[2], s=40)
ax.scatter(pos_vecinos[:,0],pos_vecinos[:,1],pos_vecinos[:,2])
# ax.scatter(posx[:,0],posx[:,1],posx[:,2])
plt.show()


deltax = pos_vecinos

presion_vecino = vecinos[:, 10]
