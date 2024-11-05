import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mpl
from cycler import cycler
import sys
from importar_datos import importar_mag, importar_t1t2t3t4

sys.path.append("..")
from funciones import (
    datenum,
    donde,
    proyecciones,
    fechas,
    tiempos,
)

np.set_printoptions(precision=4)

"""
Este script plotea mag, swea, swia y lpw en la región de interés
Es la fig principal del grl.

"""
year, month, day, doy = fechas()
ti, tf = tiempos("tiempo inicial y final en el viento solar\n")

path = f"../../../datos/clweb/{year}-{month}-{day}/"
# path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la desktop.
datos_t = np.loadtxt("../outputs/t1t2t3t4.txt")

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
t1, t2, t3, t4 = importar_t1t2t3t4(year, doy, int(ti))

t_up = t1 - 0.015
t_down = t4 + 0.015

Bnorm = np.linalg.norm(B, axis=1)


# ############ tiempos UTC
year = int(year)
month = int(month)
day = int(day)

tiempo_mag = np.array(
    [np.datetime64(datenum(year, month, day, x)) for x in t]
)  # datenum es una función mía

tm1 = donde(t, t1)
tm2 = donde(t, t2)
tm3 = donde(t, t3)
tm4 = donde(t, t4)
tm_up = donde(t, t_up)
tm_down = donde(t, t_down)

tiempo_lim = [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]

B1, B2, B3 = proyecciones(B)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)
# ####plot grl zoom proyectada en los autovectores


xup_i = posicion[tm_up, 0]
xup_f = posicion[tm1, 0]
print(f"recorrido en x upstream: {xup_f-xup_i} km")
xd_i = posicion[tm4, 0]
xd_f = posicion[tm_down, 0]
print(f"recorrido en x downstream: {xd_f-xd_i} km")

fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()

ax1.xaxis.set_major_formatter(xfmt)
ax1 = plt.subplot2grid((2, 1), (1, 0))
ax1.plot(tiempo_mag, B[:, 0], label="Bx MSO")
ax1.plot(tiempo_mag, B[:, 1], label="By MSO")
ax1.plot(tiempo_mag, B[:, 2], label="Bz MSO")
ax1.set_ylabel("B components (nT)")

ax2 = plt.subplot2grid((2, 1), (0, 0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, Bnorm)
plt.ylabel("|B| (nT)")
ax2.set_title(f"MAVEN MAG SWEA SWIA {year}-{month}-{day}")

plt.show()
