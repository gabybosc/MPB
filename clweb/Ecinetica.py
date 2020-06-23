import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_vel_swica

import sys

sys.path.append("..")
from funciones import fechas, tiempos, donde


"""
Calcula la velocidad media de los protones del viento solar de SWICA en la región
upstream proyectada sobre la normal (preguntar si está bien esta región o no)
y luego calcula Ecin.
Se puede mejorar haciendo que en los datos también esté la densidad y calculando
la media en esa región en lugar de anotandola a mano.
"""


np.set_printoptions(precision=4)

mp = 1.67e-27  # kg

year, month, day, doy = fechas()
ti, t1 = tiempos("tiempo inicial de mediciones y tiempo t1")
# t1 = 18.2167
normal = [
    float(x) for x in input("Enter the normal in N.NN N.NN N.NN format\n").split()
]
densidad = float(input("Enter the density upstream from the MPB\n"))

swia, t, vel, vel_norm = importar_vel_swica(year, month, day, ti, t1)
inicio_up = donde(t, t1 - 0.015)
fin_up = donde(t, t1)

v_para = np.dot(vel, normal)


# v_mean_sheath =
v_mean_upstream = np.mean(v_para[inicio_up:fin_up])  # km/s
print("La velocidad media en la región upstream es ", v_mean_upstream)

Ecin = 1 / 2 * mp * (v_mean_upstream * 1e3) ** 2 * (densidad * 1e6)  # J/m3

print(
    "La energía cinética en la región upstream a partir de los parámetros promedio es ",
    Ecin,
)
