import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys
import glob
import os

sys.path.append("..")
from funciones import fechas, t_clweb

"""
Hace un filtro butterworth para quitar el ruido de la señal que es de aproximadamente 180 ms.
Primero usa buttord para encontrar el orden.
Es un filtro digital con lo cual pide que las frecuencias estén normalizadas
respecto de la frec de Nyquist. En este caso Nyquist es 32Hz/2 = 16Hz.
Como las frecuencias están normalizadas, da lo mismo usar f o w.
N es el orden del filtro, en general voy a querer que N sea cercano a 10.

Comentado está otro filtro más fuerte.
"""

# path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"  # {hora}/"  # path a los datos desde la laptop
path = f"../../../datos/clweb/"  # desde la desktop

# Find all .asc files in the specified path
asc_files = glob.glob(os.path.join(path, "**/MAG.asc"), recursive=True)


def filtro(Tseg, B):
    fs = 1 / Tseg / 16  # f normalizada, da lo mismo si es omega o frec
    fp = 0.3 * fs  # fp < fs para que sea pasabajos
    N, Wn = signal.buttord(fp, fs, 3, 50)
    b, a = signal.butter(N, Wn, "low")
    Bx_filtrado = signal.filtfilt(b, a, B[:, 0])
    By_filtrado = signal.filtfilt(b, a, B[:, 1])
    Bz_filtrado = signal.filtfilt(b, a, B[:, 2])
    B_filtrado = np.transpose([Bx_filtrado, By_filtrado, Bz_filtrado])

    return (B_filtrado, fp, fs)


Tseg = 180e-3  # 180ms

for asc_file in asc_files:
    mag = np.loadtxt(asc_file)
    t = t_clweb(mag)

    M = np.size(t)  # el numero de datos

    # el campo
    B = mag[:, 6:9]

    Bnorm = mag[:, -1]

    B_filtrado, fp, fs = filtro(Tseg, B)

    # Extract directory path from the asc_file
    directory = os.path.dirname(asc_file)

    # Construct full path for saving MAG_filtrado.npy in the same directory
    save_path = os.path.join(directory, "MAG_filtrado.npy")

    arr = np.column_stack(
        (B_filtrado, np.transpose(np.linalg.norm(B_filtrado, axis=1)))
    )
    np.save(save_path, arr)
