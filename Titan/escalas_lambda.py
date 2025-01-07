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

from funciones import donde, autovectores, error, array_datenums, UTC_to_hdec

"""
tdec, Bx, By, Bz, modulo B,pos x, pos y, pos z, distancia km
"""

# path = "../../../datos/Titan/t96_tswis_1s.ascii"
path_t_hires = "../../../datos/Titan/t96_kso_hires.txt"  # solo para el tiempo
path_B_hires = "../../../datos/Titan/t96_kso_hires_filt.gz"  # filtrado, para B
tiempo = np.genfromtxt(path_t_hires, skip_header=1, dtype="str", usecols=[1])
t_hires = np.array([UTC_to_hdec(x) for x in tiempo])

Bfull = np.loadtxt(path_B_hires)

i = donde(t_hires, 0.54)
f = donde(t_hires, 0.582)

t = t_hires[i:f]
B = Bfull[i:f, :]
B_norm = np.linalg.norm(B, axis=1)


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
    with open(f"../outputs/cociente_lambdas_titan.txt", "w") as file:
        file.write(
            "La primera columna es el período de ciclotrón, la primera fila son las escalas, la segunda columna es el tiempo central.\n"
        )
        np.savetxt(file, matriz, delimiter="\t", fmt="%.6e")

    return B, t, escalas, cociente, tiempo_central


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

t1, t2, t3, t4 = 0.54175182, 0.55058123, 0.57651763, 0.58203602
timestamps = array_datenums(2013, 11, 30, np.array([t1, t2, t3, t4]))
# xfmt = md.DateFormatter("%H:%M:%S")

fig = plt.figure()
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((1, 2), (0, 0))
for vert in timestamps:
    ax1.axvline(x=vert, color="k")
plot_datetime(2013, 11, 30, t_cut, B_cut, "red", "-", 1, 1)
ax1.set_ylabel(r"|$\Delta B$|/ B")

ax2 = plt.subplot2grid((1, 2), (0, 1), sharex=ax1)
imshow_UTC(2013, 11, 30, tiempo_central, cociente, escalas_plot, "inferno", 3)
for vert in timestamps:
    ax2.axvline(x=vert, color="k")
multi = MultiCursor(fig.canvas, (ax1, ax2), color="black", lw=1)
plt.show(block=False)

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

fig = plt.figure()
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
ax1 = plt.subplot2grid((1, 1), (0, 0))
imshow_UTC(2013, 11, 30, tiempo_central, cociente, escalas_plot, "inferno", 3)
ax1.axvline(x=md.date2num(selected_time))
ax1.axhline(y=20)
ax1.set_title(
    "Heatmap del cociente de lambdas en distintas escalas temporales\npara el día 01-12-2013"
)
ax1.set_ylabel("Radio (s)")
ax1.set_xlabel("Tiempo en el que está centrado (hh:mm:ss)")
plt.show()
#
# # 01-dic-2013: central: 24:33:25, radio = 20s
