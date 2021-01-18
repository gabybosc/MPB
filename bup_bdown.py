import numpy as np
import matplotlib.pyplot as plt
from funciones import (
    fechas,
    tiempos,
    donde,
)
from funciones_plot import plot_datetime
from importar_datos import importar_mag_1s, importar_t1t2t3t4

np.set_printoptions(precision=4)

"""
Calcula diferentes valores de B_upstream y B_downstream (en diferentes lapsos)
y los grafica. Sirve para elegir qué límites tomar.
"""

year, month, day, doy = fechas()
ti, tf = tiempos()
mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)

M = len(t)  # el numero de datos

inicio = donde(t, ti)
fin = donde(t, tf)

t1, t2, t3, t4 = importar_t1t2t3t4()

Bu = np.zeros((180, 4))
for i in range(180):
    paso = t1 - 180 + 0.0001 * i
    inicio_up = donde(t, paso)
    fin_up = donde(t, t1)
    B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT
    Bu[i, 0:3] = B_upstream
    Bu[i, 3] = np.linalg.norm(B_upstream)

plt.figure(1)
plot_datetime(year, month, day, t[fin_up - 180 : fin_up], Bu[:, 0])
plot_datetime(year, month, day, t[fin_up - 180 : fin_up], Bu[:, 1], colour="C1")
plot_datetime(year, month, day, t[fin_up - 180 : fin_up], Bu[:, 2], colour="C2")
plt.xlabel("t (UTC)")
plt.ylabel("componentes Bup")

plt.figure(2)
plot_datetime(year, month, day, t[fin_up - 180 : fin_up], Bu[:, 3])
plt.xlabel("t (UTC)")
plt.ylabel("|B_up|")


Bd = np.zeros((200, 4))
for i in range(200):
    paso = t4 + 0.0001 * i
    inicio_down = donde(t, t4)
    fin_down = donde(t, paso)
    B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT
    Bd[i, 0:3] = B_downstream
    Bd[i, 3] = np.linalg.norm(B_downstream)
plt.figure(3)
plot_datetime(year, month, day, t[fin_down - 200 : fin_down], Bd[:, 0])
plot_datetime(year, month, day, t[fin_down - 200 : fin_down], Bd[:, 1], colour="C1")
plot_datetime(year, month, day, t[fin_down - 200 : fin_down], Bd[:, 2], colour="C2")
plt.xlabel("t (UTC)")
plt.ylabel("componentes Bdown")

plt.figure(4)
plot_datetime(year, month, day, t[fin_down - 200 : fin_down], Bd[:, 3])
plt.xlabel("t (UTC)")
plt.ylabel("|Bdown|")


plt.show(block=False)
