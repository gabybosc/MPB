import numpy as np
import sys
from MVA_sin_spreadsheet import MVA, ajuste, normal_coplanar
from importar_datos import importar_mag, importar_lpw
import matplotlib.pyplot as plt

sys.path.append("..")

from funciones import (
    angulo,
    ancho_mpb,
    corrientes,
    find_nearest,
    find_nearest_final,
    find_nearest_inicial,
    fechas,
    tiempos,
    donde,
)
from funciones_metodos import bootstrap

"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
"""

year, month, day, doy = fechas()
ti_MVA, tf_MVA = tiempos("Intervalo del MVA")
ti, tf = ti_MVA - 0.5, tf_MVA + 0.5  # tiempos("Región de análisis (no MVA)")

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
lpw, t_lpw, e_density, flag = importar_lpw(year, month, day, ti, tf)
x3, B_cut, t_cut, posicion_cut = MVA(year, month, day, ti_MVA, tf_MVA)
normal_ajuste, t1, t2, t3, t4 = ajuste(year, month, day, doy, ti_MVA, tf_MVA)

M = len(t)
M_cut = len(t_cut)
normal_boot, phi, delta_B3, out, out_phi = bootstrap(1000, B_cut)

# buscamos el ángulo entre las normales
angulo_mva = angulo(normal_ajuste, x3)
angulo_boot = angulo(
    normal_boot, x3
)  # Es importante que los vectoers estén normalizados!


##############
# Calculo la velocidad de la nave
v_punto = np.zeros((M_cut - 1, 3))
norma_v = np.zeros(M_cut - 1)
for i in range(len(v_punto)):
    v_punto[i, :] = (posicion_cut[i + 1, :] - posicion_cut[i]) / (1 / 32)
    # en km/s, tiene resolución de 32Hz
    norma_v[i] = np.linalg.norm(v_punto[i, :])
# la velocidad promedio
v_media = np.mean(v_punto, axis=0)

# ahora quiero ver si la nave atraviesa perpendicularmente a la MPB

angulo_v_fit = angulo(normal_ajuste, v_media)
angulo_v_mva = angulo(x3, v_media)
angulo_v_boot = angulo(normal_boot, v_media)

B_medio_vectorial = np.mean(B_cut, axis=0)
angulo_B_fit = angulo(normal_ajuste, B_medio_vectorial)
angulo_B_mva = angulo(x3, B_medio_vectorial)
angulo_B_boot = angulo(normal_boot, B_medio_vectorial)
######
# Espesor de la MPB

x_14_fit, x_23_fit = ancho_mpb(t1, t2, t3, t4, normal_ajuste, v_media)
x_14_MVA, x_23_MVA = ancho_mpb(t1, t2, t3, t4, x3, v_media)
x_14_boot, x_23_boot = ancho_mpb(t1, t2, t3, t4, normal_boot, v_media)

###########
# análisis de corrientes

inicio_up = donde(t, t1 - 0.015)
fin_up = donde(t, t1)
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(t, t4)
fin_down = donde(t, t4 + 0.015)
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

n_coplanar = normal_coplanar(B_upstream, B_downstream)

omega = angulo(B_upstream, B_downstream)

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, x_23_MVA)
J_s_fit, J_v_fit = corrientes(normal_ajuste, B_upstream, B_downstream, x_23_fit)
J_s_boot, J_v_boot = corrientes(normal_boot, B_upstream, B_downstream, x_23_boot)

fuerza_mva = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3 #en t4
fuerza_fit = np.cross(J_v_fit * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3 #en t4
fuerza_boot = np.cross(J_v_boot * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3


e_density = lpw[:, -1]
t_lpw = lpw[:, 3] + lpw[:, 4] / 60 + lpw[:, 5] / 3600

ti_lpw = np.where(t_lpw == find_nearest(t_lpw, t2))[0][0]
tf_lpw = np.where(t_lpw == find_nearest(t_lpw, t3))[0][0]
n_e = np.nanmean(e_density[ti_lpw:tf_lpw])  # hace el mean ignorando los nans #cm⁻³
n_e = n_e * 1e6  # m⁻³
if np.isnan(n_e):
    n_e = 1e7
    print("LPW no tiene datos de densidad, asumí n_e = 1E7")
q_e = 1.6e-19  # carga electron #C

E_Hall = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
E_Hall_fit = np.cross(J_v_fit * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
E_Hall_boot = np.cross(J_v_boot * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
