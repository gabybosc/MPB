import numpy as np
from os import listdir
import glob as glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import scipy.signal as signal

"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos en los cuales se cumplen:
SZA < 45º (o incluso <30º)
Altitud entre 300 y 1300 km
Z_MSO > 0 (si después se pueden volver a bajar los datos y hacer Z_pc > 0, mejor)

Primero filtra los datos de campo utilizando un filtro pasa bajo con ventana Butterworth y frecuencia de corte
de 0.01 Hz de orden 3 a fin de atenuar las variaciones del campo de frecuencias mayores a 0.125 Hz.

tenemos datos desde 10/2014 hasta 02/2018
"""

# path = glob.glob('../../../MAVEN/mag_1s/2016/*/*.sts')
# cantidad_datos = len(path)
# calendario_2014 = np.zeros(cantidad_datos)

mag = np.loadtxt('../../../MAVEN/mag_1s/2016/03/mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148)

mag[:,1]

B = np.zeros((len(mag[:,0]), 3))
for j in range(7,10):
    B[:,j-7] = mag[:, j]

B_norm = np.linalg.norm(B, axis=1)

b,a = signal.butter(3,0.01,btype='lowpass')
filtered = signal.filtfilt(b, a, B_norm)
peaks = signal.find_peaks(filtered, 40) #todos los picos mayores a 40 nT

posicion = np.zeros((len(mag[:,0]), 3))
for j in range(11,14):
    posicion[:,j-11] = mag[:, j]

"""
Clasificación por SZA, es el que menos varía. Si el SZA medio es < 45, probablemente todos los SZA lo sean.
"""
B_norm = np.linalg.norm(B, axis=1)

b,a = signal.butter(3,0.01,btype='lowpass')
filtered = signal.filtfilt(b, a, B_norm)
peaks = signal.find_peaks(filtered, 40)

SZA = np.zeros(len(peaks[0]))

for j in range(len(peaks[0])):
    SZA[j] = np.arccos(np.clip(np.dot(posicion[peaks[0][j]]/np.linalg.norm(posicion[peaks[0][j]]), [1,0,0]), -1.0, 1.0))* 180/np.pi

altitud = np.linalg.norm(posicion[peaks[0]], axis=1) - 3390
