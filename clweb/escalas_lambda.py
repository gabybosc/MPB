import numpy as np
import time as time
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
import matplotlib.dates as md
import os
from socket import gethostname

sys.path.append("..")

from funciones import find_nearest, Mij, fechas, array_datenums
from importar_datos import importar_mag
from funciones_plot import imshow_UTC, plot_datetime


plt.ion()
"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre.
En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro.
El mapa de colores me va a dar el valor del cociente.
"""
year, month, day, doy = fechas()
hora = input("Hora en HH\n")

tiempos_txt = np.loadtxt("../outputs/t1t2t3t4.txt")
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

# ti = float(input('Tiempo inicial del barrido\n'))
# tf = float(input('Tiempo final del barrido\n'))
ti = tiempos[0]
tf = tiempos[3]

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

tiempo_central = np.zeros(
    int((tf - ti) * 3600)
)  # la cantidad de segundos entre tf y ti
tiempo_central[0] = ti
for i in range(len(tiempo_central) - 1):
    tiempo_central[i + 1] = (
        tiempo_central[i] + 1 / 3600
    )  # el tiempo central se va barriendo cada 5 segundos

# print(f'tiempo final efectivo = {tiempo_central[-1]}')

escalas = np.zeros(60)
escalas[0] = 1 / 3600  # la escala más chica es de 1s
for i in range(len(escalas) - 1):
    escalas[i + 1] = escalas[i] + 1 / 3600

# print(f'escala mayor = {escalas[-1]*3600}s')

M = len(t)  # el numero de datos


program_starts = time.time()
cociente = np.zeros((len(tiempo_central), len(escalas)))
for i in range(len(tiempo_central)):
    for j in range(len(escalas)):
        inicio = np.where(t == find_nearest(t, tiempo_central[i] - escalas[j]))[0][0]
        fin = np.where(t == find_nearest(t, tiempo_central[i] + escalas[j]))[0][0]

        # ahora empieza el MVA con los datos que elegí
        B_cut = B[inicio : fin + 1, :]

        M_ij = Mij(B_cut)

        # ahora quiero los autovectores y autovalores
        [lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

        # Los ordeno de mayor a menor
        idx = lamb.argsort()[::-1]
        lamb = lamb[idx]
        x = x[:, idx]
        # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
        x1 = x[:, 0]
        x2 = x[:, 1]
        x3 = x[:, 2]
        if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
            x3 = -x3

        cociente[i, j] = lamb[1] / lamb[2]

program_ends = time.time()

print(f"El loop tardó {program_ends-program_starts:.2f} s")

m = 1.67e-27  # kg
q = 1.6e-19  # C
periodo_ciclotron = 2 * np.pi * m / (q * np.linalg.norm(B_cut, axis=1)) * 1e9  # en s
periodo_diezmado = np.zeros(len(tiempo_central))
k = len(periodo_ciclotron) / len(tiempo_central)
for i in range(len(periodo_diezmado)):
    periodo_diezmado[i] = periodo_ciclotron[int(i * k)]

matriz = np.zeros((len(tiempo_central) + 1, len(escalas) + 2))
matriz[0, 2:] = escalas
matriz[1:, 0] = periodo_diezmado
matriz[1:, 1] = tiempo_central
matriz[1:, 2:] = cociente

with open(f"../outputs/cociente_lambdas_d{doy}_t{hora}.txt", "w") as file:
    # with open(f'outputs/cociente_lambdas_salida_d{doy}_t{hora}.txt','w') as file:
    file.write(
        "La primera columna es el período de ciclotron, la primera fila son las escalas, la segunda columna es el tiempo central.\n"
    )
    for i in range(len(matriz[:, 0])):
        for j in range(len(matriz[0, :])):
            file.write(f"{matriz[i,j]}\t")
        file.write("\n")
# print(f'{l / len(dates) * 100}%')

"""
Plotea los archivos que devuelve escalas_lambda.py
"""

timestamps = array_datenums(year, month, day, tiempos)  # lo convierto a datenum

Bnorm = np.linalg.norm(B, axis=0)

datos = np.loadtxt(f"../outputs/cociente_lambdas_d{doy}_t{hora}.txt", skiprows=1)

periodo_ciclotron = datos[1:, 0]
tiempo_central = datos[1:, 1]
escalas = datos[0, 2:] * 3600
cociente = np.transpose(datos[1:, 2:])

ti = tiempo_central[0] - 0.5
tf = tiempo_central[-1] + 0.5
n = int(ti * 32 * 3600)


hh = mag[:, 3]
mm = mag[:, 4]
ss = mag[:, 5]

t = hh + mm / 60 + ss / 3600  # hdec

M = np.size(t)  # el numero de datos

inicio = np.where(t == find_nearest(t, ti))[0][0]
fin = np.where(t == find_nearest(t, tf))[0][0]

B_cut = Bnorm[inicio:fin]
t_cut = t[inicio:fin]

inicio_MVA = np.where(t == find_nearest(t, tiempo_central[0]))[0][0]
fin_MVA = np.where(t == find_nearest(t, tiempo_central[-1]))[0][0]
B_MVA = Bnorm[inicio_MVA:fin_MVA]
t_MVA = t[inicio_MVA:fin_MVA]


xfmt = md.DateFormatter("%H:%M:%S")

plt.figure(1)
imshow_UTC(year, month, day, tiempo_central, cociente, escalas, "inferno", 3)
# plot_datetime(year, month, day,t_MVA, B_MVA, 'cyan', '-', 1, 0.5) #no sé por qué se superpone mal, tiene mal los tiempos.
for tt in timestamps:
    plt.axvline(x=tt, color="g")  # plotea los tiempos t1t2t3t4
plt.title(
    f"Heatmap del cociente de lambdas en distintas escalas temporales \n  para el día {day}-{month}-{year}"
)
plt.xlabel("Tiempo en el que está centrado (hh:mm:ss)")
plt.ylabel("Radio (s) \n |B| (nT)")


plt.figure(2)
plot_datetime(year, month, day, t_cut, B_cut, "red", "-", 1, 1)
plt.ylabel("|B| (nT)")
plt.xlabel("Tiempo UTC (hh:mm:ss)")

fig = plt.figure(
    3, constrained_layout=True
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((1, 2), (0, 0))
plot_datetime(year, month, day, t_cut, B_cut, "red", "-", 1, 1)
ax1.set_ylabel(r"|$\Delta B$|/ B")

ax2 = plt.subplot2grid((1, 2), (0, 1), sharex=ax1)
imshow_UTC(year, month, day, tiempo_central, cociente, escalas, "inferno", 3)


for ax in [ax1, ax2]:
    ax.set_xlim(timestamps[0], timestamps[-1])
    ax.xaxis.set_major_formatter(xfmt)
    ax.grid()
    for tt in timestamps:
        ax.axvline(x=tt, color="g")  # plotea los tiempos t1t2t3t4

multi = MultiCursor(fig.canvas, (ax1, ax2), color="c", lw=1)

plt.show()
