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
def velocidad(posicion_cut):
    M = len(posicion_cut)
    v_punto = np.zeros((M - 1, 3))
    for i in range(len(v_punto)):
        v_punto[i, :] = (posicion_cut[i + 1, :] - posicion_cut[i]) / (1 / 32)
        # en km/s, tiene resolución de 32Hz
    # la velocidad promedio
    v_media = np.mean(v_punto, axis=0)
    return v_media


def corte(t, ti, tf, vector):
    inicio = donde(t, ti)
    fin = donde(t, tf)
    if vector.ndim > 1:
        vector_cut = np.nanmean(vector[inicio:fin, :], axis=0)
    else:
        vector_cut = np.nanmean(vector[inicio:fin])
    return vector_cut


def fuerza(J, B):
    f = np.cross(J * 1e-9, B * 1e-9)  # en N/m³
    return f


v_media = velocidad(posicion_cut)
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

B_upstream = corte(t, t1 - 0.015, t1, B)
B_downstream = corte(t, t4, t4 + 0.015, B)

n_coplanar = normal_coplanar(B_upstream, B_downstream)

omega = angulo(B_upstream, B_downstream)

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, x_23_MVA)
J_s_fit, J_v_fit = corrientes(normal_ajuste, B_upstream, B_downstream, x_23_fit)
J_s_boot, J_v_boot = corrientes(normal_boot, B_upstream, B_downstream, x_23_boot)

inicio_down = donde(t, t1 - 0.015)
fuerza_mva = fuerza(J_v_MVA, B[inicio_down, :])
fuerza_fit = fuerza(J_v_fit, B[inicio_down, :])
fuerza_boot = fuerza(J_v_boot, B[inicio_down, :])


e_density = lpw[:, -1]
t_lpw = lpw[:, 3] + lpw[:, 4] / 60 + lpw[:, 5] / 3600

n_e = corte(t_lpw, t2, t3, e_density)
n_e = n_e * 1e6  # m⁻³
if np.isnan(n_e):
    n_e = 1e7
    print("LPW no tiene datos de densidad, asumí n_e = 1E7")
q_e = 1.6e-19  # carga electron #C

E_Hall = fuerza_mva / (q_e * n_e)  # V/m
E_Hall_fit = fuerza_fit / (q_e * n_e)  # V/m
E_Hall_boot = fuerza_boot / (q_e * n_e)  # V/m
