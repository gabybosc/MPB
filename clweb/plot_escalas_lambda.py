"""
Plotea los archivos que devuelve escalas_lambda_inbound.py.py
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import MultiCursor
import matplotlib.dates as md
import os
import sys
from socket import gethostname

sys.path.append("..")

from funciones import find_nearest, array_datenums, fechas
from funciones_plot import imshow_UTC, plot_datetime
from importar_datos import importar_t1t2t3t4

plt.ion()

year, month, day, doy = fechas()
hora = input("Hora (HH) \n")

tiempos = importar_t1t2t3t4(year, month, day, int(hora))
timestamps = array_datenums(year, month, day, tiempos)  # lo convierto a datenum

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
if len(str(int(hora))) == 1:
    t_f = "0" + str(int(hora))
else:
    t_f = int(hora)

# if gethostname() == "magneto2":
path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/{hora}/"
if not os.path.exists(path):  # si no existe, usa int(tf)
    path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"

# path = f"../../../datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la laptop
if os.path.isfile(path + "mag_filtrado.txt"):
    mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
    M = len(mag[:, 0])  # el numero de datos
    B = mag[:, :3]

    Bnorm = mag[:, -1]
    mag = np.loadtxt(path + "MAG.asc")
    Bxyz_paraperp = mag[:, 6:9]
else:
    mag = np.loadtxt(path + "MAG.asc")
    M = len(mag[:, 0])  # el numero de datos
    B = mag[:, 6:9]
    Bnorm = np.linalg.norm(B, axis=1)

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

fig = plt.figure(3)
# Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
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
