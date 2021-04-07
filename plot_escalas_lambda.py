"""
Plotea los archivos que devuelve escalas_lambda.py
"""

import matplotlib.pyplot as plt
import numpy as np
from funciones import fechas, array_datenums, donde
from funciones_plot import imshow_UTC, plot_datetime
from importar_datos import importar_mag

plt.ion()

year, month, day, doy = fechas()
hora = input("hora\n")

tiempos_txt = np.loadtxt("outputs/t1t2t3t4.txt")
for i in range(len(tiempos_txt)):
    if (
        int(year) == int(tiempos_txt[i, 0])
        and int(doy) == int(tiempos_txt[i, 1])
        and int(hora) == int(tiempos_txt[i, 2])
    ):
        tiempos = [
            tiempos_txt[i, 2],
            tiempos_txt[i, 3],
            tiempos_txt[i, 4],
            tiempos_txt[i, 5],
        ]
timestamps = array_datenums(year, month, day, tiempos)  # lo convierto a datenum


datos = np.loadtxt(f"outputs/cociente_lambdas_d{doy}_t{hora}.txt", skiprows=1)

periodo_ciclotron = datos[1:, 0]
tiempo_central = datos[1:, 1]
escalas = datos[0, 2:] * 3600
cociente = np.transpose(datos[1:, 2:])

ti = tiempo_central[0] - 0.5
tf = tiempo_central[-1] + 0.5
n = int(ti * 32 * 3600)

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

inicio = donde(t, ti)
fin = donde(t, tf)

B_norm = np.linalg.norm(B, axis=1)

inicio_MVA = donde(t, tiempo_central[0])
fin_MVA = donde(t, tiempo_central[-1])
B_MVA = B_norm[inicio_MVA:fin_MVA]
t_MVA = t[inicio_MVA:fin_MVA]

plt.figure()
imshow_UTC(year, month, day, tiempo_central, cociente, escalas, "inferno")
# plot_datetime(year, month, day,t_MVA, B_MVA, 'cyan', '-', 1, 0.5) #no sé por qué se superpone mal, tiene mal los tiempos.
for tt in timestamps[1:3]:
    plt.axvline(x=tt, color="g")  # plotea los tiempos t2 y t3
plt.title("Heatmap del cociente de lambdas en distintas escalas temporales")
plt.xlabel("Tiempo en el que está centrado (hh:mm:ss)")
plt.ylabel("Radio (s) \n |B| (nT)")

plt.figure()
plot_datetime(year, month, day, t, B_norm, "red", "-", 1, 1)
plt.ylabel("|B| (nT)")
plt.xlabel("Tiempo UTC (hh:mm:ss)")

"""
pensar criterio para la escala maxima
"""
