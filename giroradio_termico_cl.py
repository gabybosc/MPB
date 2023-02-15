import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from funciones import fechas, donde
from funciones_plot import onpick1
from importar_datos import importar_mag_1s, importar_swia, importar_fila

np.set_printoptions(precision=4)

"""
Calcula el giroradio térmico para el 16 de marzo usando los datos de clweb
"""

# #########CONSTANTES
mp = 1.67e-27  # masa del proton en kg
# mp = 1.5e-10 #masa del proton en joules/c^2
kB = 1.38e-23  # cte de Boltzmann en J/K
q_e = 1.602e-19  # carga del electrón en C

# ########## DATOS
# year, month, day, doy = fechas()
# ti = int(input("Hora del cruce (HH)\n"))
# tf = ti + 1
# in_out = input("Inbound? (y/n)\n")

year = 2016
month = "03"
day = 16
doy = "076"
ti = 17.93333
tf = 18.1

mag, t_mag, B, posicion = importar_mag_1s(year, month, day, ti, tf)

swia = np.loadtxt("../../datos/temp_swica_ms_600-1400.asc")
t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec
temperature = swia[:, 6]

inicio = donde(t, ti)
fin = donde(t, tf)

t_swia = t[inicio:fin]
temp = temperature[inicio:fin]

# quiero diezmar el tiempo y el campo para que sean los mismos que tiene swia
idx = np.zeros(len(t_swia))
for i in range(len(idx)):
    idx[i] = donde(t_mag, t_swia[i])
idx = idx.astype(int)

t_diezmado = t_mag[idx]  # lo diezmó

B_cut = B[idx]
posicion_cut = posicion[idx]

####################

B_avg = np.empty((len(idx), 3))
B_avg_normalized = np.empty((len(idx), 3))
temp_para_xyz = np.empty((len(idx), 3))

for i in range(len(idx) - 1):
    B_avg[i, :] = np.mean(B_cut[i : i + 30, :], axis=0) * 1e-5  # lo paso a gauss
    B_avg_normalized[i, :] = B_avg[i, :] / np.linalg.norm(B_avg[i, :])  # adimensional

thermal_gyroradius = np.empty(len(temp))
B_avg_norm = np.linalg.norm(B_avg, axis=1)
for i in range(len(temp)):
    thermal_gyroradius[i] = 1.02e02 * np.sqrt(temp[i]) / B_avg_norm[i] * 1e-5  # km

print(
    f"thermal ion gyroradius mean = {np.nanmean(thermal_gyroradius[:-1], axis=0):1.3g} km"
)  # nanmean ignora los nans
