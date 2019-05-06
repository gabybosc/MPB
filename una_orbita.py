import numpy as np
from os import listdir
import glob as glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import scipy.signal as signal
from funciones import find_nearest, set_axes_equal
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
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
mag = np.loadtxt('../../datos/MAG_1s/mvn_mag_l2_2016092ss1s_20160401_v01_r01.sts', skiprows=148) #en mi compu

B = np.zeros((len(mag[:,0]), 3))
for j in range(7,10):
    B[:,j-7] = mag[:, j]

    B_norm = np.linalg.norm(B, axis=1)

    posicion = np.zeros((len(mag[:,0]), 3))
    for j in range(11,14):
        posicion[:,j-11] = mag[:, j]

        orbitas = posicion / 3390 #radios marcianos

una_vuelta = int(len(orbitas)/5)

# Ajuste de Vignes:
x0 = 0.78
e = 0.9
L = 0.96

theta = np.linspace(0, 3*np.pi/4, 100)
phi = np.linspace(0, np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(theta))

X = x0 + r * np.cos(THETA)
Y = r * np.sin(THETA) * np.cos(PHI)
Z = r * np.sin(THETA) * np.sin(PHI)

R = np.transpose(np.array([X.flatten(),Y.flatten(),Z.flatten()]))

"""
Son dos loops: el loop en i barre toda la superficie y la resta para cada punto de la órbita. El loop en j agarra esa resta y ve dónde es que es mínima (busca el máximo acercamiento entre la órbita y la superficie). Luego, guarda el mínimo para cada punto de la órbita. Finalmente, busca el mínimo de mínimos.
Hace esto cada 100 puntos porque si no tarda mucho.

Necesito que haga la resta sólo para Z positivo, si no, puede crear cruces ficticios en el sur (si la posicion de maven en el sur se parece a la posición de la MPB en el norte también lo cuenta, ya que tomo la norma)
"""
resta = np.zeros((len(R),3))
idx_min = np.zeros(int(una_vuelta/100))
max_acercamiento = np.zeros(int(una_vuelta/100))
X_MSO = posicion[int(una_vuelta*2):int(una_vuelta*3), 0]
Z_MSO = posicion[int(una_vuelta*2):int(una_vuelta*3),2]
orbita = posicion[int(una_vuelta*2):int(una_vuelta*3),:] /3390
for j in range(una_vuelta-100):
    if j%100 == 0 and Z_MSO[j] >0 and X_MSO[j]>0: #hace la cuenta cada 100 pasos.
        for i in range(len(R)):
            resta[i, :] = np.abs(orbita[j,:] - R[i,:])
        A = np.linalg.norm(resta, axis=1)
        idx_min[int(j/100)] = np.argmin(A)
        max_acercamiento[int(j/100)] = A[int(idx_min[int(j/100)])]
minimo = np.where( max_acercamiento==np.min(max_acercamiento[np.nonzero(max_acercamiento)]))[0][0]
print(max_acercamiento[minimo])

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.set_xlabel(r'$X_{MSO} (R_m)$')
ax.set_ylabel(r'$Y_{MSO} (R_m)$')
ax.set_zlabel(r'$Z_{MSO} (R_m)$')
ax.set_aspect('equal')
# ax.plot(orbita[:,0], orbita[:,1], orbita[:,2], color='C2', label='Órbita')
ax.plot(orbita[:,0], orbita[:,1], orbita[:,2], color='C2', label='Órbita')
ax.scatter(orbita[minimo*100, 0],orbita[minimo*100, 1],orbita[minimo*100, 2])
plot = ax.plot_surface(
    X, Y, Z, rstride=4, cstride=4, alpha=0.5, edgecolor='none', cmap=plt.get_cmap('Blues_r'))
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)
ax.legend()
set_axes_equal(ax) #para que tenga forma de esfera la esfera


"""
Clasificación por SZA, es el que menos varía. Si el SZA medio es < 45, probablemente todos los SZA lo sean.
"""

idx = minimo * 100
SZA = np.arccos(np.clip(np.dot(posicion[int(idx)]/np.linalg.norm(posicion[int(idx)]), [1,0,0]), -1.0, 1.0))* 180/np.pi

altitud = np.linalg.norm(posicion[int(idx), :]) - 3390

Z_MSO = posicion[int(idx), 2]

if Z_MSO > 0:
    print('Tiene SZA = {0:1.3g}, altitud = {1:1.3g} y Z_MSO > 0'.format(SZA, altitud))
