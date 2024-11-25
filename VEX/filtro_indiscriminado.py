import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys
from socket import gethostname
import os

from _importar_datos import importar_MAG_pds

sys.path.append("..")
from funciones import fechas, day_to_doy, donde

"""
Hace un filtro butterworth para quitar el ruido de la señal que es de aproximadamente 180 ms.
Primero usa buttord para encontrar el orden.
Es un filtro digital con lo cual pide que las frecuencias estén normalizadas
respecto de la frec de Nyquist. En este caso Nyquist es 128Hz/2 = 64Hz.
Como las frecuencias están normalizadas, da lo mismo usar f o w.
N es el orden del filtro, en general voy a querer que N sea cercano a 10.

Como es de 128Hz, va a cortarlo a 32Hz tomando un punto de cada cuatro. 
En este caso entonces, Nyquist pasa a ser 16Hz.
"""
year = 2014

# lista = np.loadtxt(f"../outputs/VEX{year}_menor65.txt", dtype=str)

# # i = int(input("indice en lista\n"))
# for i in range(len(lista)):
#     print(i)
#     l = lista[i]
#     year, month, day = l[0].split("-")
#     year, doy = day_to_doy(year, month, day)

for doy in [91, 99, 109, 116]:
    print(doy)
    ti, tf = 0, 24  # tiempos()
    tiempo, campo, posicion = importar_MAG_pds(year, str(doy).zfill(3), ti, tf)
    t = tiempo[::4]
    B = campo[::4]
    pos = posicion[::4]

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

    Tseg = 180e-3  # 180ms
    B_filtrado, fp, fs = filtro(Tseg)
    # plt.plot(t, Bnorm, label="sin filtro")
    # plt.plot(
    #     t,
    #     np.linalg.norm(B_filtrado, axis=1),
    #     linewidth=1,
    #     label=f"fs = {fs:.3g}, fp = {fp:.3g}",
    # )
    # plt.legend()
    # plt.show()

    if gethostname() == "DESKTOP-2GS0QF2":
        os.chdir(f"G:/")
        filt = f"VEX{year}/VEX_mag_filtrado_{year}{str(doy).zfill(3)}.gz"
        np.savetxt(filt, np.vstack((t, B_filtrado.T, pos.T)))
        os.chdir("C:/Users/RainbowRider/Documents/GitHub/MPB/VEX/")
