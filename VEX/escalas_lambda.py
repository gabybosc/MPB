"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre.
En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro.
El mapa de colores me va a dar el valor del cociente.
"""

import numpy as np
import time as time
import os
from importar_datos import importar_MAG_pds
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.widgets import MultiCursor


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
    return (ti, tf)


def MVA(t, ti, tf, B):
    inicio = donde(t, ti)
    fin = donde(t, tf)

    B_cut = B[inicio : fin + 1, :]

    M_ij = Mij(B_cut)

    avec, lamb = autovectores(M_ij)

    cociente = lamb[1] / lamb[2]

    B1 = np.dot(B, avec[0])
    B2 = np.dot(B, avec[1])
    B3 = np.dot(B, avec[2])

    return (cociente, B1, B2, B3, B_cut)


def escalas_lambda(year, doy, ti, tf):
    t, B, posicion = importar_MAG_pds(year, doy, ti, tf)

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

    with open(f"../outputs/cociente_lambdas_VEX_d{doy}.txt", "w") as file:
        # with open(f'outputs/cociente_lambdas_salida_d{doy}_t{hora}.txt','w') as file:
        file.write(
            "La primera columna es el período de ciclotron, la primera fila son las escalas, la segunda columna es el tiempo central.\n"
        )
        for i in range(len(matriz[:, 0])):
            for j in range(len(matriz[0, :])):
                file.write(f"{matriz[i,j]}\t")
            file.write("\n")
    # print(f'{l / len(dates) * 100}%')
    return (B, t, escalas, cociente, tiempo_central)


year, month, day, doy = fechas()
ti, tf = importar_titf(year, doy)

path = f"../outputs/cociente_lambdas_VEX_d{doy}.txt"
if os.path.isfile(path):
    datos = np.loadtxt(path, skiprows=1)
    mag, t, B, posicion = importar_MAG_pds(year, doy, ti, tf)

    periodo_ciclotron = datos[1:, 0]
    tiempo_central = datos[1:, 1]
    escalas = datos[0, 2:]
    cociente = datos[1:, 2:]

else:
    B, t, escalas, cociente, tiempo_central = escalas_lambda(year, doy, ti, tf)


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


xfmt = md.DateFormatter("%H:%M:%S")

fig = plt.figure(2)
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((1, 2), (0, 0))
# ax1.plot(t_cut, B_cut)
plot_datetime(year, month, day, t_cut, B_cut, "red", "-", 1, 1)
ax1.set_ylabel(r"|$\Delta B$|/ B")

ax2 = plt.subplot2grid((1, 2), (0, 1), sharex=ax1)
imshow_UTC(year, month, day, tiempo_central, cociente, escalas_plot, "inferno", 3)
multi = MultiCursor(fig.canvas, (ax1, ax2), color="black", lw=1)
plt.show()
