import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys

sys.path.append("..")
from funciones import fechas

"""
Hace un filtro butterworth para quitar el ruido de la señal que es de aproximadamente 180 ms.
Primero usa buttord para encontrar el orden.
Es un filtro digital con lo cual pide que las frecuencias estén normalizadas
respecto de la frec de Nyquist. En este caso Nyquist es 32Hz/2 = 16Hz.
Como las frecuencias están normalizadas, da lo mismo usar f o w.
N es el orden del filtro, en general voy a querer que N sea cercano a 10.

Comentado está otro filtro más fuerte.
"""

year, month, day, doy = fechas()
# year, month, day, doy = 2016, "03", 16, 76

path = f"../../../datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la laptop
mag = np.loadtxt(path + "MAG.asc")

hh = mag[:, 3]
mm = mag[:, 4]
ss = mag[:, 5]

t = hh + mm / 60 + ss / 3600  # hdec

M = np.size(t)  # el numero de datos

# el campo
B = np.zeros((M, 3))
B = mag[:, 6:9]

Bnorm = mag[:, -1]

# de 18.2539  a 18.255: 18:15:14 a 18:15:18

Tseg = 180e-3  # 180ms
fs = 1 / Tseg / 16  # f normalizada, da lo mismo si es omega o frec
fp = 3 / 16  # fp < fs para que sea pasabajos
N, Wn = signal.buttord(fp, fs, 3, 50)
b, a = signal.butter(N, Wn, "low")
Bx_filtrado = signal.filtfilt(b, a, B[:, 0])
By_filtrado = signal.filtfilt(b, a, B[:, 1])
Bz_filtrado = signal.filtfilt(b, a, B[:, 2])

B_filtrado = np.linalg.norm([Bx_filtrado, By_filtrado, Bz_filtrado], axis=0)

plt.plot(t, Bnorm, label="sin filtro")
plt.plot(t, B_filtrado, linewidth=0.5, label=f"fs = {fs:.3g}, fp = {fp:.3g}")
plt.legend()
plt.show(block=False)

# filtro más fuerte
# Tseg = 4
# fp = 1/Tseg /16 #fp < fs para que sea pasabajos
# fs = 1 /16 #f normalizada, da lo mismo si es omega o frec
# N, Wn = signal.buttord(fp, fs, 3, 50)
# b,a = signal.butter(N, Wn,'low')
# Bx_filtrado = signal.filtfilt(b, a, B[:,0])
# By_filtrado = signal.filtfilt(b, a, B[:,1])
# Bz_filtrado = signal.filtfilt(b, a, B[:,2])
#
# B_filtrado = np.linalg.norm([Bx_filtrado,By_filtrado,Bz_filtrado], axis=0)
#
# plt.plot(t,Bnorm, label='sin filtro')
# plt.plot(t,B_filtrado,linewidth = 1, label = f'fs = {fs:.3g}, fp = {fp:.3g}')
# plt.legend()
# plt.show(block=False)

happy = input("If happy press Y\n")
if happy == "y" or happy == "Y":
    with open(path + "mag_filtrado.txt", "w") as file:
        file.write(
            f"Los datos de MAG filtrados para frecuencia fp = {fp*16:.3g}, fs = {fs*16:.3g}.\n"
        )
        file.write(f"Bx  By  Bz  B.\n")
        for i in range(M):
            file.write(
                f"{Bx_filtrado[i]}\t{By_filtrado[i]}\t{Bz_filtrado[i]}\t{B_filtrado[i]}\t"
            )
            file.write("\n")
