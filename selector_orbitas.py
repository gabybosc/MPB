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

# path = glob.glob('../../../MAVEN/mag_1s/2018/*/*.sts')
path =  glob.glob('../../datos/MAG_1s/*.sts')
cantidad_datos = len(path)
calendario_2018 = np.zeros((5, cantidad_datos)) #la primera columna es el día del año, la segunda es el número de orbita, la tercera dice si cumple el SZA, la cuarta la altitud y la quinta el ZMSO


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

#Importamos los datos y vemos uno por uno qué pasa.
for i,j in enumerate(path): #loop en todos los archivos .sts para cada año. i me da el índice en la lista, j me da el archivo
    mag = np.loadtxt(j, skiprows=160)

    calendario_2018[0,i] = mag[1,1]

    posicion = np.zeros((len(mag[:,0]), 3))
    for k in range(11,14):
        posicion[:,k-11] = mag[:, k]
    orbita = posicion/3390
    una_vuelta = int(len(orbita)/5)


    B = np.zeros((len(mag[:,0]), 3))
    for k in range(7,10):
        B[:,k-7] = mag[:, k]

    """
    Son dos loops: el loop en i barre toda la superficie y la resta para cada punto de la órbita. El loop en j agarra esa resta y ve dónde es que es mínima (busca el máximo acercamiento entre la órbita y la superficie). Luego, guarda el mínimo para cada punto de la órbita. Finalmente, busca el mínimo de mínimos.
    Hace esto cada 100 puntos y sólo donde Z y X MSO son positivas, total es donde está mi cruce. (esto además me disminuye los falsos positivos)

    """
    orbitas = [orbita[:una_vuelta], orbita[una_vuelta:una_vuelta*2], orbita[una_vuelta*2:una_vuelta*3], orbita[una_vuelta*3:una_vuelta*4], orbita[una_vuelta*4:]]
    resta = np.zeros((len(R),3))
    for indice, l in enumerate(orbitas):
        pos = l * 3390
        X_MSO = pos[:, 0]
        Z_MSO = pos[:, 2]
        idx_min = np.zeros(int(una_vuelta/100))
        max_acercamiento = np.zeros(int(una_vuelta/100))
        minimo = 0
        for k in range(int(una_vuelta)-100):
            if k%100 == 0 and Z_MSO[k] > 0 and X_MSO[k] > 0:
                for m in range(len(R)):
                    resta[m, :] = l[k,:] - R[m,:]
                A = np.linalg.norm(resta, axis=1)
                idx_min[int(k/100)] = np.argmin(A)
                max_acercamiento[int(k/100)] = A[int(idx_min[int(k/100)])]
        minimo = np.where( max_acercamiento==np.min(max_acercamiento[np.nonzero(max_acercamiento)]))[0][0] #busca el minimo que no sea cero
        calendario_2018[1,i] = indice+1

        """
        Clasificación por SZA
        """
        idx = minimo * 100

        SZA = np.arccos(np.clip(np.dot(pos[int(idx)]/np.linalg.norm(posicion[int(idx)]), [1,0,0]), -1.0, 1.0))* 180/np.pi

        if SZA < 30:
            calendario_2018[2,i] = 1

        """
        Clasificación por altitud
        """
        altitud = np.linalg.norm(pos[int(idx),:]) - 3390

        if altitud < 1300 and altitud > 300:
            calendario_2018[3,i] = 1

        """
        Clasificación por Z_MSO
        """
        Z_MSO = pos[int(idx),2]
        if Z_MSO > 0:
            calendario_2018[4,i] = 1

np.savetxt('tiempos_2018.txt', np.transpose(calendario_2018), fmt='%10d' ,header= "       dia        SZA        altitud       Z_MSO", newline="\r\n")

    # SZA = np.arccos(np.clip(np.dot(pos[int(idx)]/np.linalg.norm(pos[int(idx)]), [1,0,0]), -1.0, 1.0))* 180/np.pi
    # altitud = np.linalg.norm(pos[int(idx), :]) - 3390
    # Z_MSO_cruce = pos[int(idx),2]
    #
    # print('Tiene SZA = {0:1.3g}, altitud = {1:1.3g} y Z_MSO = {2:1.3g}'.format(SZA, altitud, Z_MSO_cruce))
    # # if Z_MSO > 0:
    #     print('Tiene SZA = {0:1.3g}, altitud = {1:1.3g} y Z_MSO > 0'.format(SZA, altitud))
