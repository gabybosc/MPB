"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre.
En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro.
El mapa de colores me va a dar el valor del cociente.
"""

import numpy as np
import time as time
import os
from importar_datos import importar_mag, importar_t1t2t3t4
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.widgets import MultiCursor
from datetime import datetime, timedelta

import sys

sys.path.append("..")
from funciones import Mij, fechas, donde, autovectores, array_datenums
from funciones_plot import imshow_UTC, plot_datetime, hodograma


def importar_titf(year, doy):
    path = "../outputs/VEX_times.txt"
    datos = np.loadtxt(path)
    for d in datos:
        if int(d[0]) == int(year) and int(d[1]) == int(doy):
            ti = d[2]
            tf = d[-1]
    return ti, tf


def escalas_lambda(t, B):
    tf, ti = t[-1], t[0]
    tiempo_central = np.linspace(ti, tf, int((tf - ti) * 3600))
    # el tiempo central es cada un segundo en el intervalo ti, tf

    escalas = np.linspace(1 / 3600, 60 / 3600, 60)
    # las escalas van de 1 a 60 segundos

    # precomputa los cortes de B para cada tiempo y escala
    B_cuts = {}
    for esc in escalas:
        for tc in tiempo_central:
            inicio = donde(t, tc - esc)
            fin = donde(t, tc + esc)
            B_cuts[(tc, esc)] = B[inicio: fin + 1, :]

    # Hago el MVA para cada combinación de tiempo central y escalas
    cociente = np.zeros((len(tiempo_central), len(escalas)))
    for i, tc in enumerate(tiempo_central):
        for j, esc in enumerate(escalas):
            B_cut = B_cuts[(tc, esc)]
            M_ij = np.cov(B_cut.T)
            _, lamb = autovectores(M_ij)
            cociente[i, j] = lamb[1] / lamb[2]

    # Calcular períodos de ciclotrón y diezmar
    m, q = 1.67e-27, 1.6e-19  # kg, C
    periodo_ciclotron = (
            2 * np.pi * m / (q * np.linalg.norm(B_cut, axis=1)) * 1e9
    )  # en s
    periodo_diezmado = np.interp(
        np.linspace(0, len(periodo_ciclotron) - 1, len(tiempo_central)),
        np.arange(len(periodo_ciclotron)),
        periodo_ciclotron,
    )

    # crea la matriz de salida
    matriz = np.zeros((len(tiempo_central) + 1, len(escalas) + 2))
    matriz[0, 2:] = escalas
    matriz[1:, 0] = periodo_diezmado
    matriz[1:, 1] = tiempo_central
    matriz[1:, 2:] = cociente

    # lo guarda
    with open(f"../outputs/cociente_lambdas_d{doy}_t18.txt", "w") as file:
        file.write(
            "La primera columna es el período de ciclotrón, la primera fila son las escalas, la segunda columna es el tiempo central.\n"
        )
        np.savetxt(file, matriz, delimiter="\t", fmt="%.6e")

    return B, t, escalas, cociente, tiempo_central


year, month, day, doy = 2016, "03", 16, "076"
ti, tf = 18.21, 18.25

path = f"../outputs/cociente_lambdas_d{doy}_t18.txt"

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

if os.path.isfile(path):
    datos = np.loadtxt(path, skiprows=1)

    periodo_ciclotron = datos[1:, 0]
    tiempo_central = datos[1:, 1]
    escalas = datos[0, 2:]
    cociente = datos[1:, 2:]

else:
    B, t, escalas, cociente, tiempo_central = escalas_lambda(t, B)

Bnorm = np.linalg.norm(B, axis=1)
#
# escalas_plot = escalas * 3600
# cociente = np.transpose(cociente)
#
# ti = tiempo_central[0] - 0.5
# tf = tiempo_central[-1] + 0.5
# n = int(ti * 32 * 3600)
#
# inicio = donde(t, ti)
# fin = donde(t, tf)
#
# B_cut = Bnorm[inicio:fin]
# t_cut = t[inicio:fin]
#
# inicio_MVA = donde(t, tiempo_central[0])
# fin_MVA = donde(t, tiempo_central[-1])
# B_MVA = Bnorm[inicio_MVA:fin_MVA]
# t_MVA = t[inicio_MVA:fin_MVA]
#
# xfmt = md.DateFormatter("%H:%M:%S")
#
# fig = plt.figure()
# fig.subplots_adjust(
#     top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
# )
# fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop
#
# ax1 = plt.subplot2grid((1, 2), (0, 0))
# # ax1.plot(t_cut, B_cut)
# plot_datetime(year, month, day, t_cut, B_cut, "red", "-", 1, 1)
# ax1.set_ylabel(r"|$\Delta B$|/ B")
#
# ax2 = plt.subplot2grid((1, 2), (0, 1), sharex=ax1)
# imshow_UTC(year, month, day, tiempo_central, cociente, escalas_plot, "inferno", 3)
# multi = MultiCursor(fig.canvas, (ax1, ax2), color="black", lw=1)
# plt.show()


escalas_plot = escalas * 3600
cociente = np.transpose(cociente)

ti = tiempo_central[0] - 0.5
tf = tiempo_central[-1] + 0.5
n = int(ti * 32 * 3600)

inicio = donde(t, ti)
fin = donde(t, tf)

B_cut = Bnorm[inicio:fin]
t_cut = t[inicio:fin]

inicio_MVA = donde(t, tiempo_central[0])
fin_MVA = donde(t, tiempo_central[-1])
B_MVA = Bnorm[inicio_MVA:fin_MVA]
t_MVA = t[inicio_MVA:fin_MVA]

t1, t2, t3, t4 = (
    18.2167,
    18.2204,
    18.235,
    18.2476,
)
timestamps = array_datenums(2016, 3, 16, np.array([t1, t2, t3, t4]))
# xfmt = md.DateFormatter("%H:%M:%S")

fig = plt.figure()
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((1, 2), (0, 0))
for vert in timestamps:
    ax1.axvline(x=vert, color="k")
plot_datetime(2016, 3, 16, t_cut, B_cut, "red", "-", 1, 1)
ax1.set_ylabel(r"|$\Delta B$|/ B")

ax2 = plt.subplot2grid((1, 2), (0, 1), sharex=ax1)
imshow_UTC(2016, 3, 16, tiempo_central, cociente, escalas_plot, "viridis", 3)
for vert in timestamps:
    ax2.axvline(x=vert, color="k")
multi = MultiCursor(fig.canvas, [ax1, ax2], color="black", lw=1)
plt.show(block=False)

print("Click to select time: ")
val = np.asarray(plt.ginput(1))
timestamps = array_datenums(2016, 3, 16, tiempo_central)
t_graph = md.date2num(timestamps)
idx = donde(t_graph, val[0][0])

selected_time = timestamps[idx].astype("M8[s]").astype(datetime)
# Formatea el tiempo para que solo muestre HH:MM:SS
formatted_time = selected_time.strftime("%H:%M:%S")

time_minus = (selected_time - timedelta(seconds=float(val[0][1]))).strftime("%H:%M:%S")
time_plus = (selected_time + timedelta(seconds=float(val[0][1]))).strftime("%H:%M:%S")

print("Selected values: ", formatted_time, val[0][1])
print(f"Selected time - radius: {time_minus}")
print(f"Selected time + radius: {time_plus}")

fig = plt.figure()
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
ax1 = plt.subplot2grid((1, 1), (0, 0))
imshow_UTC(2016, 3, 16, tiempo_central, cociente, escalas_plot, "viridis", 3)
ax1.axvline(x=md.date2num(selected_time), c="k")
ax1.axhline(y=timedelta(seconds=float(val[0][1])).total_seconds(), c="k")
ax1.set_title(
    f"Mapa de calor del cociente de lambdas en distintas escalas temporales\npara el día {year}-{month}-{day}"
)
ax1.set_ylabel("Radio (s)")
ax1.set_xlabel("Tiempo en el que está centrado (hh:mm:ss)")
plt.show()
#
# fig = plt.figure()
# fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop
# fig.subplots_adjust(
#     top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
# )
# ax1 = plt.subplot2grid((1, 1), (0, 0))
# imshow_UTC(year, month, day, tiempo_central, cociente, escalas_plot, "inferno", 3)
# ax1.set_title(
#     "Heatmap del cociente de lambdas en distintas escalas temporales\npara el día 28-10-2008"
# )
# ax1.set_ylabel("Radio (s)")
# ax1.set_xlabel("Tiempo en el que está centrado (hh:mm:ss)")
# plt.show()

# 28-oct-2008: central: 08:32:48, radio = 9s
