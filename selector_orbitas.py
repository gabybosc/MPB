import numpy as np
from os import listdir
import glob as glob
from datetime import datetime, timedelta
import scipy.signal as signal


"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos en los cuales se cumplen:
SZA < 45º (o incluso <30º)
Altitud entre 300 y 1300 km
Z_MSO > 0 (si después se pueden volver a bajar los datos y hacer Z_pc > 0, mejor)

Nos interesa que se cumplan estas cosas en la MPB, entonces le voy a pedir que encuentre un cruce en cada día y analice estas tres condiciones
Para esto, le voy a decir qeu encuentre el mayor pico en todo el módulo del campo, habiendolo filtrado antes, y analice las tres características en ese punto.
Si se cumplen en una órbita, lo más seguro es que se cumplan en todas las del mismo día.

El filtro es pasa bajo con ventana Butterworth y frecuencia de corte de 0.01 Hz de orden 3

tenemos datos desde 10/2018 hasta 02/2018
"""

path = glob.glob('../../../MAVEN/mag_1s/2018/*/*.sts')
cantidad_datos = len(path)
calendario_2018 = np.zeros((4, cantidad_datos)) #la primera columna es el día del año, la segunda dice si cumple el SZA, la tercera la altitud y la cuarta el ZMSO


for i,j in enumerate(path): #loop en todos los archivos .sts para cada año. i me da el índice en la lista, j me da el archivo
    mag = np.loadtxt(j, skiprows=160)

    calendario_2018[0,i] = mag[1,1]

    posicion = np.zeros((len(mag[:,0]), 3))
    for k in range(11,14):
        posicion[:,k-11] = mag[:, k]

    B = np.zeros((len(mag[:,0]), 3))
    for k in range(7,10):
        B[:,k-7] = mag[:, k]

    """
    Filtra los datos del campo y busca los picos mayores a 40 nT
    """
    B_norm = np.linalg.norm(B, axis=1)

    b,a = signal.butter(3,0.01,btype='lowpass')
    filtered = signal.filtfilt(b, a, B_norm)
    peaks = signal.find_peaks(filtered, 40) #peaks[0] son los índices, peaks[1] son los valores

    """
    Clasificación por SZA
    """
    SZA = np.arccos(np.clip(np.dot(posicion[peaks[0]]/np.linalg.norm(posicion[peaks[0]]), [1,0,0]), -1.0, 1.0))* 180/np.pi

    if any(SZA < 45):
        calendario_2018[1,i] = 1

    """
    Clasificación por altitud
    """
    altitud = np.linalg.norm(posicion[peaks[0]], axis=1) - 3390

    if any(altitud < 1300) and any(altitud > 300):
        calendario_2018[2,i] = 1

    """
    Clasificación por Z_MSO
    """
    Z_MSO = posicion[peaks[0],2]
    if any(Z_MSO > 0):
        calendario_2018[3,i] = 1

np.savetxt('tiempos_2018.txt', np.transpose(calendario_2018), fmt='%10d' ,header= "       dia          SZA          altitud        Z_MSO", newline="\r\n")
