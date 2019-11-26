import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import datetime as dt
from shutil import copyfile
from funciones import datenum, fechas
import matplotlib.dates as md
import matplotlib.cm as cm

"""
Hace un filtro butterworth para quitar el ruido de la señal que es de aproximadamente 180 ms.
Primero usa buttord para encontrar el orden. Es un filtro digital con lo cual pide que las frecuencias estén normalizadas respecto de la frec de Nyquist. En este caso Nyquist es 32Hz/2 = 16Hz.
Como las frecuencias están normalizadas, da lo mismo usar f o w.
N es el orden del filtro, en general voy a querer que N sea cercano a 10.
"""

year, month, day, doy = fechas()

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
mag = np.loadtxt(path + 'MAG.asc')

hh = mag[:,3]
mm = mag[:,4]
ss = mag[:,5]

t = hh + mm/60 + ss/3600 #hdec

M = np.size(t) #el numero de datos

#el campo
B = np.zeros((M, 3))
B = mag[:, 6:9]

Bnorm = mag[:,-1]


Tseg = 180E-3 #180ms
fs = 1/Tseg /16 #f normalizada, da lo mismo si es omega o frec
fp = 3 /16
N, Wn = signal.buttord(fp, fs, 3, 50)
b,a = signal.butter(N, Wn,'low')
Bx_filtrado = signal.filtfilt(b, a, B[:,0])
By_filtrado = signal.filtfilt(b, a, B[:,1])
Bz_filtrado = signal.filtfilt(b, a, B[:,2])

B_filtrado = np.linalg.norm([Bx_filtrado,By_filtrado,Bz_filtrado], axis=0)

plt.plot(Bnorm, label='sin filtro')
plt.plot(B_filtrado,linewidth = 0.5, label = f'fs = {fs:.3g}, fp = {fp:.3g}')
plt.legend()
plt.show(block=False)

happy = input('If happy press Y\n')

if happy == 'y' or 'Y':
    with open(path + 'mag_filtrado.txt','w') as file:
        file.write(f'Los datos de MAG filtrados para frecuencia fp = {fp*16:.3g}, fs = {fs*16:.3g}.\n')
        file.write(f'Bx  By  Bz  B.\n')
        for i in range(M):
            file.write(f'{Bx_filtrado[i]}\t{By_filtrado[i]}\t{Bz_filtrado[i]}\t{B_filtrado[i]}\t')
            file.write('\n')
