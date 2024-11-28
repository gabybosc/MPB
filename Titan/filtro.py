import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys
from socket import gethostname
import os

sys.path.append("..")
from funciones import UTC_to_hdec

"""
Hace un filtro butterworth para quitar el ruido de la señal que es de aproximadamente 180 ms.
Primero usa buttord para encontrar el orden.
Es un filtro digital con lo cual pide que las frecuencias estén normalizadas
respecto de la frec de Nyquist. En este caso Nyquist es 32Hz/2 = 16Hz.
Como las frecuencias están normalizadas, da lo mismo usar f o w.
N es el orden del filtro, en general voy a querer que N sea cercano a 10.

Comentado está otro filtro más fuerte.
"""
path_hires = "../../../datos/Titan/t96_kso_hires.txt"
tiempo = np.genfromtxt(path_hires, skip_header=1, dtype="str", usecols=[1])
t = np.array([UTC_to_hdec(x) for x in tiempo])

B_total = np.genfromtxt(path_hires, skip_header=1, usecols=[-4, -3, -2, -1])
Bnorm = B_total[:, 0]
Bx, By, Bz = B_total[:, 1], B_total[:, 2], B_total[:, 3]


def filtro(Tseg):
    fs = 1 / Tseg / 16  # f normalizada, da lo mismo si es omega o frec
    fp = 0.3 * fs  # fp < fs para que sea pasabajos
    N, Wn = signal.buttord(fp, fs, 3, 50)
    b, a = signal.butter(N, Wn, "low")
    Bx_filtrado = signal.filtfilt(b, a, Bx)
    By_filtrado = signal.filtfilt(b, a, By)
    Bz_filtrado = signal.filtfilt(b, a, Bz)
    B_filtrado = np.transpose([Bx_filtrado, By_filtrado, Bz_filtrado])

    return (B_filtrado, fp, fs)


# Tseg = input("T filtro? (default 180 ms)\n")
# if not Tseg:  # si no escribe nada, toma el default
#     Tseg = 180e-3  # 180ms
# else:
#     Tseg = float(Tseg)
Tseg = 180e-3  # 180ms
B_filtrado, fp, fs = filtro(Tseg)
plt.plot(t, Bnorm, label="sin filtro")
plt.plot(
    t,
    np.linalg.norm(B_filtrado, axis=1),
    linewidth=1,
    label=f"fs = {fs:.3g}, fp = {fp:.3g}",
)
plt.legend()
plt.show()
if gethostname() == "DESKTOP-2GS0QF2":
    filt = "../../../datos/Titan/t96_kso_hires_filt.gz"
    np.savetxt(filt, B_filtrado)
