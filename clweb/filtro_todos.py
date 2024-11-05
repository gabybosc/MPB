import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys

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

year, month, day, doy = fechas()
# year, month, day, doy = 2016, "03", 16, 76
# hora = input("hora\n")

# path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"  # {hora}/"  # path a los datos desde la laptop
path = f"../../../datos/clweb/{year}-{month}-{day}/"  # desde la desktop

mag = np.loadtxt(path + "MAG.asc")

t = t_clweb(mag)

M = np.size(t)  # el numero de datos

# el campo
B = np.zeros((M, 3))
B = mag[:, 6:9]

Bnorm = mag[:, -1]


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


happy = False
while not happy:
    Tseg = input("T filtro (default 180 ms)\n")
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
    plt.show(block=False)

    happy = input("If happy enter any key\n")

with open(path + "mag_filtrado.txt", "w") as file:
    file.write(
        f"Los datos de MAG filtrados para frecuencia fp = {fp * 16:.3g}, fs = {fs * 16:.3g}.\n"
    )
    file.write(f"Bx  By  Bz  B.\n")
    for i in range(M):
        file.write(
            f"{B_filtrado[i, 0]}\t{B_filtrado[i, 1]}\t{B_filtrado[i, 2]}\t{np.linalg.norm(B_filtrado, axis=1)[i]}\t"
        )
        file.write("\n")

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

year, month, day, doy = fechas()
# year, month, day, doy = 2016, "03", 16, 76
# hora = input("hora\n")

# path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"  # {hora}/"  # path a los datos desde la laptop
path = f"../../../datos/clweb/{year}-{month}-{day}/"  # desde la desktop

# Find all .asc files in the specified path
asc_files = glob.glob(os.path.join(path, "*.asc"))


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


happy = False
while not happy:
    Tseg = input("T filtro (default 180 ms)\n")
    if not Tseg:  # si no escribe nada, toma el default
        Tseg = 180e-3  # 180ms
    else:
        Tseg = float(Tseg)

    for asc_file in asc_files:
        mag = np.loadtxt(asc_file)
        t = t_clweb(mag)

        M = np.size(t)  # el numero de datos

        # el campo
        B = mag[:, 6:9]

        Bnorm = mag[:, -1]

        B_filtrado, fp, fs = filtro(Tseg, B)

        plt.plot(t, Bnorm, label="sin filtro")
        plt.plot(
            t,
            np.linalg.norm(B_filtrado, axis=1),
            linewidth=1,
            label=f"fs = {fs:.3g}, fp = {fp:.3g}",
        )
        plt.legend()
        plt.title(os.path.basename(asc_file))
        plt.show(block=False)

        with open(asc_file.replace(".asc", "_filtrado.txt"), "w") as file:
            file.write(
                f"Los datos de MAG filtrados para frecuencia fp = {fp * 16:.3g}, fs = {fs * 16:.3g}.\n"
            )
            file.write(f"Bx  By  Bz  B.\n")
            for i in range(M):
                file.write(
                    f"{B_filtrado[i, 0]}\t{B_filtrado[i, 1]}\t{B_filtrado[i, 2]}\t{np.linalg.norm(B_filtrado, axis=1)[i]}\t"
                )
                file.write("\n")

    happy = input("If happy enter any key\n")
