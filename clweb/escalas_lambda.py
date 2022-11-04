import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
import matplotlib.dates as md
import os
from socket import gethostname
from escalas_lambda_func import escalas_lambda, MVA

sys.path.append("..")

from funciones import (
    find_nearest,
    Mij,
    fechas,
    array_datenums,
    donde,
    tiempos,
    UTC_to_hdec,
)
from importar_datos import importar_mag, importar_t1t2t3t4
from funciones_plot import imshow_UTC, plot_datetime, hodograma

plt.ion()
"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre.
En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro.
El mapa de colores me va a dar el valor del cociente.
"""
year, month, day, doy = fechas()
hora = input("Hora en HH\n")

tt = importar_t1t2t3t4(year, month, day, int(hora))


"""
Si ya tengo el escalas_lambda, no lo hace de nuevo, sino que lo abre.
Si no, corre la función
"""
if os.path.isfile(f"../outputs/cociente_lambdas_d{doy}_t{hora}.txt"):
    datos = np.loadtxt(f"../outputs/cociente_lambdas_d{doy}_t{hora}.txt", skiprows=1)
    mag, t, B, posicion = importar_mag(year, month, day, tt[0], tt[3])

    periodo_ciclotron = datos[1:, 0]
    tiempo_central = datos[1:, 1]
    escalas = datos[0, 2:]
    cociente = datos[1:, 2:]

else:
    B, t, escalas, cociente, tiempo_central = escalas_lambda(
        year, month, day, doy, hora, tt
    )

"""
Plotea los archivos que devuelve escalas_lambda.py
"""

timestamps = array_datenums(year, month, day, tt)  # lo convierto a datenum

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

# plt.figure(1)
# imshow_UTC(year, month, day, tiempo_central, cociente, escalas_plot, "inferno", 3)
# # plot_datetime(year, month, day,t_MVA, B_MVA, 'cyan', '-', 1, 0.5) #no sé por qué se superpone mal, tiene mal los tiempos.
# for tt in timestamps:
#     plt.axvline(x=tt, color="g")  # plotea los tiempos t1t2t3t4
# plt.title(
#     f"Heatmap del cociente de lambdas en distintas escalas temporales \n  para el día {day}-{month}-{year}"
# )
# plt.xlabel("Tiempo en el que está centrado (hh:mm:ss)")
# plt.ylabel("Radio (s) \n |B| (nT)")


fig = plt.figure(2)
# Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
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


for ax in [ax1, ax2]:
    ax.set_xlim(timestamps[0], timestamps[-1])
    ax.xaxis.set_major_formatter(xfmt)
    ax.grid()
    for tt in timestamps:
        ax.axvline(x=tt, color="g")  # plotea los tiempos t1t2t3t4

multi = MultiCursor(fig.canvas, (ax1, ax2), color="c", lw=1)

plt.show()

"""
Hace el MVA para hacer el hodograma entre los tiempos que diga
"""
t_c = UTC_to_hdec(input("tiempo central HH:MM:SS\n"))
t_pm = int(input("radio de t\n")) / 3600
ratio, B1, B2, B3, B_cut = MVA(t, t_c - t_pm, t_c + t_pm, B)

hodograma(B1, B2, B3)
print("cociente de lambdas =", ratio)
plt.show()
