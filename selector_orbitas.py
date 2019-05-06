import numpy as np
from os import listdir
import glob as glob
from datetime import datetime, timedelta
import scipy.signal as signal
from funciones import find_nearest


"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos en los cuales se cumplen:
SZA < 45º (o incluso <30º)
Altitud entre 300 y 1300 km
Z_MSO > 0 (si después se pueden volver a bajar los datos y hacer Z_pc > 0, mejor)

Voy a hacer el fit de vignes para la MPB y ver si cuando cruza el fit cada día se cumplen estas tres cosas.
También le voy a pedir que me diga la hora a la cual pasa.

tenemos datos desde 10/2018 hasta 02/2018
"""

path = glob.glob('../../../MAVEN/mag_1s/2018/*/*.sts')
cantidad_datos = len(path)
calendario_2018 = np.zeros((4, cantidad_datos)) #la primera columna es el día del año, la segunda dice si cumple el SZA, la tercera la altitud y la cuarta el ZMSO


# Ajuste de Vignes:
x0 = 0.78
e = 0.9
L = 0.96

theta = np.linspace(0, np.pi *3/4, 100)
phi = np.linspace(0, 2 * np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(THETA))

#Importamos los datos y vemos uno por uno qué pasa.
for i,j in enumerate(path): #loop en todos los archivos .sts para cada año. i me da el índice en la lista, j me da el archivo
    mag = np.loadtxt(j, skiprows=160)

    calendario_2018[0,i] = mag[1,1]

    posicion = np.zeros((len(mag[:,0]), 3))
    for k in range(11,14):
        posicion[:,k-11] = mag[:, k]

    B = np.zeros((len(mag[:,0]), 3))
    for k in range(7,10):
        B[:,k-7] = mag[:, k]

    cruce_MPB = np.where(posicion == find_nearest(posicion, r))
    """
    Filtra los datos del campo y busca los picos mayores a 40 nT
    """
    B_norm = np.linalg.norm(B, axis=1)

    b,a = signal.butter(3,0.01,btype='lowpass')
    filtered = signal.filtfilt(b, a, B_norm)
    peaks = signal.find_peaks(filtered, 40) #peaks[0] son los índices, peaks[1] son los valores
    mpb = peaks[0]-500 #la posicion media de la mpb es aprox 500 puntos antes de cada pico
    """
    Clasificación por SZA
    """
    SZA = np.zeros(len(mpb))

    for j in range(len(mpb)):
        SZA[j] = np.arccos(np.clip(np.dot(posicion[mpb[j]]/np.linalg.norm(posicion[mpb[j]]), [1,0,0]), -1.0, 1.0))* 180/np.pi

    if any(SZA < 45):
        calendario_2018[1,i] = 1

    """
    Clasificación por altitud
    """
    altitud = np.linalg.norm(posicion[mpb], axis=1) - 3390

    if any(altitud < 1300) and any(altitud > 300):
        calendario_2018[2,i] = 1

    """
    Clasificación por Z_MSO
    """
    Z_MSO = posicion[mpb,2]
    if any(Z_MSO > 0):
        calendario_2018[3,i] = 1

np.savetxt('tiempos_2018.txt', np.transpose(calendario_2018), fmt='%10d' ,header= "       dia        SZA        altitud       Z_MSO", newline="\r\n")

"""
Ahora, una vez que tengo todo, vamos a ver qué días cumplen la condición.
Si suma 3, es que cumple las tres cosas.
"""

clasific = np.sum(calendario_2018[1:,:], axis =0) #si suma 3, es un día que cumple todo.

fechas_buenas = np.zeros(len(clasific))

for i in range(len(clasific)):
    if clasific[i] == 3:
        fechas_buenas[i] = calendario_2018[0,i]

#para que me de un array ordenado y sin ceros:
fechas_buenas.sort()
fechas_buenas = np.trim_zeros(fechas_buenas)

np.savetxt('fechas_buenas_2018.txt', fechas_buenas, fmt='%10d' ,header= "Las fechas de 2018 en las cuales se cumple todo", newline="\r\n")
