import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys
from importar_datos import importar_MAG_pds, importar_MAG_clweb

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

ti, tf = 0, 24  # tiempos()
t, B, pos = importar_MAG_pds(year, doy, ti, tf)

Bnorm = np.linalg.norm(B, axis=1)


def filtro(Tseg):
    fs = 1 / Tseg / 16  # f normalizada, da lo mismo si es omega o frec
    fp = 0.3 * fs  # fp < fs para que sea pasabajos
    N, Wn = signal.buttord(fp, fs, 3, 50)
    b, a = signal.butter(N, Wn, "low")
    Bx_filtrado = signal.filtfilt(b, a, B[:, 0])
    By_filtrado = signal.filtfilt(b, a, B[:, 1])
    Bz_filtrado = signal.filtfilt(b, a, B[:, 2])
    B_filtrado = np.transpose([Bx_filtrado, By_filtrado, Bz_filtrado])

    return (B_filtrado, fp, fs)


Tseg = input("T filtro? (default 180 ms)\n")
if not Tseg:  # si no escribe nada, toma el default
    Tseg = 180e-3  # 180ms
else:
    Tseg = float(Tseg)

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

if input("save? y/n\n") == "y":
    np.savetxt(
        f"../../../datos/VEX/filtrados/VEX_mag_filtrado_{year}{doy}.gz",
        B_filtrado,
    )
