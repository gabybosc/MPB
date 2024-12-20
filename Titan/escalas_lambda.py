"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre.
En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro.
El mapa de colores me va a dar el valor del cociente.
"""
import numpy as np
import time as time
import os
import matplotlib.dates as md
import matplotlib.pyplot as plt
import sys
from matplotlib.widgets import MultiCursor
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime, timedelta

sys.path.append("..")
from funciones_plot import imshow_UTC, plot_datetime, hodograma

from funciones import donde, Mij, autovectores, error, array_datenums

"""
tdec, Bx, By, Bz, modulo B,pos x, pos y, pos z, distancia km
"""

path = "../../../datos/Titan/t96_tswis_1s.ascii"
datos = np.loadtxt(path)
tiempo = datos[:, 0]
i = donde(tiempo, 24.54)
f = donde(tiempo, 24.582)

t, B, B_norm, posicion = datos[i:f, 0], datos[i:f, 1:4], datos[i:f, 4], datos[i:f, 5:8]


def MVA(t, ti, tf, B):
    inicio = donde(t, ti)
    fin = donde(t, tf)

    B_cut = B[inicio: fin + 1, :]

    M_ij = Mij(B_cut)

    avec, lamb = autovectores(M_ij)

    cociente = lamb[1] / lamb[2]

    B1 = np.dot(B, avec[0])
    B2 = np.dot(B, avec[1])
    B3 = np.dot(B, avec[2])

    return cociente, B1, B2, B3, B_cut


def escalas_lambda(t, B):
    tf = t[-1]
    ti = t[0]
    tiempo_central = np.zeros(
        int((tf - ti) * 3600)
    )  # la cantidad de segundos entre tf y ti
    tiempo_central[0] = ti
    for i in range(len(tiempo_central) - 1):
        tiempo_central[i + 1] = (
                tiempo_central[i] + 1 / 3600
        )  # el tiempo central se va barriendo cada 5 segundos

    escalas = np.zeros(60)
    escalas[0] = 1 / 3600  # la escala más chica es de 1s
    for i in range(len(escalas) - 1):
        escalas[i + 1] = escalas[i] + 1 / 3600

    cociente = np.zeros((len(tiempo_central), len(escalas)))
    for i in range(len(tiempo_central)):
        for j in range(len(escalas)):
            ratio, B1, B2, B3, B_cut = MVA(
                t, tiempo_central[i] - escalas[j], tiempo_central[i] + escalas[j], B
            )
            cociente[i, j] = ratio

    m = 1.67e-27  # kg
    q = 1.6e-19  # C
    periodo_ciclotron = (
            2 * np.pi * m / (q * np.linalg.norm(B_cut, axis=1)) * 1e9
    )  # en s
    periodo_diezmado = np.zeros(len(tiempo_central))
    k = len(periodo_ciclotron) / len(tiempo_central)
    for i in range(len(periodo_diezmado)):
        periodo_diezmado[i] = periodo_ciclotron[int(i * k)]

    matriz = np.zeros((len(tiempo_central) + 1, len(escalas) + 2))
    matriz[0, 2:] = escalas
    matriz[1:, 0] = periodo_diezmado
    matriz[1:, 1] = tiempo_central
    matriz[1:, 2:] = cociente

    with open(f"../outputs/cociente_lambdas_titan.txt", "w") as file:
        # with open(f'outputs/cociente_lambdas_salida_d{doy}_t{hora}.txt','w') as file:
        file.write(
            "La primera columna es el período de ciclotron, la primera fila son las escalas, la segunda columna es el tiempo central.\n"
        )
        for i in range(len(matriz[:, 0])):
            for j in range(len(matriz[0, :])):
                file.write(f"{matriz[i, j]}\t")
            file.write("\n")
    # print(f'{l / len(dates) * 100}%')
    return (B, t, escalas, cociente, tiempo_central)


B, t, escalas, cociente, tiempo_central = escalas_lambda(t, B)

Bnorm = np.linalg.norm(B, axis=1)

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

# xfmt = md.DateFormatter("%H:%M:%S")

fig = plt.figure()
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((1, 2), (0, 0))
# ax1.plot(t_cut, B_cut)
plot_datetime(2013, 11, 30, t_cut, B_cut, "red", "-", 1, 1)
ax1.set_ylabel(r"|$\Delta B$|/ B")

ax2 = plt.subplot2grid((1, 2), (0, 1), sharex=ax1)
imshow_UTC(2013, 11, 30, tiempo_central, cociente, escalas_plot, "inferno", 3)
multi = MultiCursor(fig.canvas, (ax1, ax2), color="black", lw=1)

print("Click to select time: ")
val = np.asarray(plt.ginput(1))
timestamps = array_datenums(2013, 11, 30, tiempo_central)
t_graph = md.date2num(timestamps)
idx = donde(t_graph, val[0][0])

selected_time = timestamps[idx].astype("M8[s]").astype(datetime)
# Formatea el tiempo para que solo muestre HH:MM:SS
formatted_time = selected_time.strftime("%H:%M:%S")

time_minus = (selected_time - timedelta(seconds=float(val[0][1]))).strftime("%H:%M:%S")
time_plus = (selected_time + timedelta(seconds=float(val[0][1]))).strftime("%H:%M:%S")

print("Selected values: ", formatted_time, val[0][1])
print(f"Selected time - 27s: {time_minus}")
print(f"Selected time + 27s: {time_plus}")
plt.show()

# fig = plt.figure()
# fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop
# fig.subplots_adjust(
#     top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
# )
# ax1 = plt.subplot2grid((1, 1), (0, 0))
# imshow_UTC(2013, 11, 30, tiempo_central, cociente, escalas_plot, "inferno", 3)
# # ax1.axvline(x="24:33:15")
# ax1.axhline(y=10)
# ax1.set_title(
#     "Heatmap del cociente de lambdas en distintas escalas temporales\npara el día 30-11-2013"
# )
# ax1.set_ylabel("Radio (s)")
# ax1.set_xlabel("Tiempo en el que está centrado (hh:mm:ss)")
# plt.show()
#
# # 01-dic-2013: central: 24:33:15, radio = 10s
