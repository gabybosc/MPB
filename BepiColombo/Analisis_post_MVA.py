import numpy as np
import matplotlib.pyplot as plt
from MVA import MVA
from leer_datos import importar_bepi

import sys

sys.path.append("../")
from funciones import (
    donde,
    corrientes,
    angulo,
    fechas,
    tiempos,
)

"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente

NO NECESITA YA ACTUALIZAR LA FECHA EN GDOCS
"""

year, month, day = 2021, "08", 10  # fechas()
ti_MVA, tf_MVA = 13.911111, 13.922222
t1, t2, t3, t4 = [13.8974, 13.9078, 13.9283, 13.9469]

t, B, posicion = importar_bepi(ti_MVA - 1, tf_MVA + 1)

(
    x3,
    normal_boot,
    normal_fit,
    inicio,
    fin,
    B_cut,
    t1,
    t2,
    t3,
    t4,
    B_medio_vectorial,
    fila,
    hoja_parametros,
    hoja_mva,
    hoja_boot,
    hoja_fit,
) = MVA(ti_MVA, tf_MVA, t, B, posicion)

#########
# buscamos el ángulo entre las normales
angulo_mva = angulo(normal_fit, x3)
angulo_boot = angulo(normal_boot, x3)


##############
# Calculo la velocidad de la nave
v_punto = np.zeros((fin - inicio, 3))
norma_v = np.zeros(fin - inicio)
posicion_cut = posicion[inicio : fin + 1, :]
t_cut = t[inicio : fin + 1] * 3600  # en segundos
for i in range(fin - inicio):
    v_punto[i, :] = (posicion[inicio + 1, :] - posicion[inicio]) / (
        t_cut[i + 1] - t_cut[i]
    )  # en km/s
    norma_v[i] = np.linalg.norm(v_punto[i, :])
# veamos que no cambia mucho punto a punto, usemos la norma
diff = max(norma_v) - min(norma_v)
# la velocidad promedio
v_media = np.array(
    [np.mean(v_punto[:, 0]), np.mean(v_punto[:, 1]), np.mean(v_punto[:, 2])]
)

# ahora quiero ver si la nave atraviesa perpendicularmente a la MPB
v_media_norm = v_media / np.linalg.norm(
    v_media
)  # la normalizo por cuestion de  cuentas nomas
angulo_v_fit = angulo(normal_fit, v_media_norm)
angulo_v_mva = angulo(x3, v_media_norm)
angulo_v_boot = angulo(normal_boot, v_media_norm)

B_intermedio = B_medio_vectorial / np.linalg.norm(
    B_medio_vectorial
)  # la normalizo por cuestion de  cuentas nomas
angulo_B_fit = angulo(normal_fit, B_intermedio)
angulo_B_mva = angulo(x3, B_intermedio)
angulo_B_boot = angulo(normal_boot, B_intermedio)

######
# Espesor de la MPB
# ahora veamos v_para
deltat_14 = (t4 - t1) * 3600
deltat_23 = (t3 - t2) * 3600

v_para = np.dot(v_media, normal_fit) * normal_fit
x_14_fit = v_para * deltat_14  # en km# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x_23_fit = v_para * deltat_23

# si ahora proyecto sobre la normal de la MVA
v_para_MVA = np.dot(v_media, x3) * x3
x_14_MVA = v_para_MVA * deltat_14  # en km
x_23_MVA = v_para_MVA * deltat_23

# si ahora proyecto sobre la normal del bootstrap
v_para_boot = np.dot(np.array([-2.487, 0.479, 2.836]), normal_boot) * normal_boot
x_14_boot = v_para_boot * deltat_14  # en km
x_23_boot = v_para_boot * deltat_23

# plot_velocidades(X1, Y1, Z1, R, normal_fit, x3, v_media, v_para, v_para_MVA)

###########
# giroradio


#########
# análisis de corrientes

inicio_up = donde(t, t1 - 0.015)  # las 18:12:00
fin_up = donde(t, t1)  # las 18:13:00
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(t, t4)  # las 18:14:51
fin_down = donde(t, t4 + 0.015)  # las 18:15:52
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT


omega = np.arccos(
    np.dot(B_upstream, B_downstream)
    / (np.linalg.norm(B_upstream) * np.linalg.norm(B_downstream))
)

mu = 4 * np.pi * 1e-7  # Tm/A

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, np.linalg.norm(x_23_MVA))
J_s_fit, J_v_fit = corrientes(
    normal_fit, B_upstream, B_downstream, np.linalg.norm(x_23_fit)
)
J_s_boot, J_v_boot = corrientes(
    normal_boot, B_upstream, B_downstream, np.linalg.norm(x_23_boot)
)


fuerza_mva = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3 #en t4
fuerza_fit = np.cross(J_v_fit * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3 #en t4
fuerza_boot = np.cross(J_v_boot * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3


ti_lpw = donde(t_lpw, t2)
tf_lpw = donde(t_lpw, t3)
n_e = np.nanmean(e_density[ti_lpw:tf_lpw])  # hace el mean ignorando los nans #cm⁻³
n_e = n_e * 1e6  # m⁻³
if np.isnan(n_e):
    n_e = 1e7
    print("LPW no tiene datos de densidad, asumí n_e = 1E7")
q_e = 1.6e-19  # carga electron #C

E_Hall = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
E_Hall_fit = np.cross(J_v_fit * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
E_Hall_boot = np.cross(J_v_boot * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m


#
# plot_FLorentz(X1, Y1, Z1, R, J_v, B_upstream, B_downstream, fuerza_mva, x3)
#
# plt.figure()
# plt.plot(t[inicio_up:fin_down], fuerza_mva)
# plt.plot(t[inicio_up:fin_down], fuerza_ajuste)
# plt.legend(['Fx', 'Fy', 'Fz'])
# plt.xlabel('Tiempo')
# plt.ylabel('Fuerza (N/m^3)')
#
plt.show(block=False)


###########
# Ahora guardamos todo en la spreadsheet


##########
# Parámetros

# hoja_parametros.update_acell(f"Z{fila}", f"{omega*180/np.pi:.3g}")
# hoja_parametros.update_acell(f"S{fila}", f"{np.linalg.norm(v_media):.3g}")
#
# cell_vel = hoja_parametros.range(f"P{fila}:R{fila}")
# for i, cell in enumerate(cell_vel):
#     cell.value = round(v_media[i], 2)
# hoja_parametros.update_cells(cell_vel)
#
# cell_Bup = hoja_parametros.range(f"T{fila}:V{fila}")
# for i, cell in enumerate(cell_Bup):
#     cell.value = round(B_upstream[i], 2)
# hoja_parametros.update_cells(cell_Bup)
#
# cell_Bdown = hoja_parametros.range(f"W{fila}:Y{fila}")
# for i, cell in enumerate(cell_Bdown):
#     cell.value = round(B_downstream[i], 2)
# hoja_parametros.update_cells(cell_Bdown)
#
#
# # La hoja del MVA
#
# hoja_mva.update_acell(f"W{fila}", f"{angulo_v_mva * 180/np.pi:.3g}")
# hoja_mva.update_acell(f"X{fila}", f"{angulo_B_mva * 180/np.pi:.3g}")
# hoja_mva.update_acell(f"Y{fila}", f"{np.linalg.norm(x_23_MVA):.3g}")
# hoja_mva.update_acell(f"Z{fila}", f"{np.linalg.norm(x_14_MVA):.3g}")
# # hoja_mva.update_acell(f'AA{fila}', f'{:.3g}')
# # hoja_mva.update_acell(f'AB{fila}', f'{:.3g}')
#
# hoja_mva.update_acell(f"AF{fila}", f"{np.linalg.norm(J_s_MVA)*1E-6:.3g}")
# hoja_mva.update_acell(f"AJ{fila}", f"{np.linalg.norm(J_v_MVA):.3g}")
# hoja_mva.update_acell(f"AK{fila}", f"{np.linalg.norm(fuerza_mva):.3g}")
# hoja_mva.update_acell(f"AO{fila}", f"{np.linalg.norm(E_Hall)*1E3:.3g}")  # mV/m
#
# cell_Js = hoja_mva.range(f"AC{fila}:AE{fila}")
# for i, cell in enumerate(cell_Js):
#     cell.value = round(J_s_MVA[i] * 1e-6, 3)
# hoja_mva.update_cells(cell_Js)
#
# cell_Jv = hoja_mva.range(f"AG{fila}:AI{fila}")
# for i, cell in enumerate(cell_Jv):
#     cell.value = round(J_v_MVA[i], 3)
# hoja_mva.update_cells(cell_Jv)
#
# cell_EH = hoja_mva.range(f"AL{fila}:AN{fila}")
# for i, cell in enumerate(cell_EH):
#     cell.value = round(E_Hall[i] * 1e3, 3)
# hoja_mva.update_cells(cell_EH)

# # La hoja del bootstrap

# hoja_boot.update_acell(f"M{fila}", f"{angulo_v_boot * 180/np.pi:.3g}")
# hoja_boot.update_acell(f"N{fila}", f"{angulo_B_boot * 180/np.pi:.3g}")
# hoja_boot.update_acell(f"O{fila}", f"{np.linalg.norm(x_23_boot):.3g}")
# hoja_boot.update_acell(f"P{fila}", f"{np.linalg.norm(x_14_boot):.3g}")
# # hoja_boot.update_acell(f'Q{fila}', f'{:.3g}')
# # hoja_boot.update_acell(f'R{fila}', f'{:.3g}')
#
# hoja_boot.update_acell(f"V{fila}", f"{np.linalg.norm(J_s_boot)*1E-6:.3g}")
# hoja_boot.update_acell(f"Z{fila}", f"{np.linalg.norm(J_v_boot):.3g}")
# hoja_boot.update_acell(f"AA{fila}", f"{np.linalg.norm(fuerza_boot):.3g}")
# hoja_boot.update_acell(f"AE{fila}", f"{np.linalg.norm(E_Hall_boot)*1E3:.3g}")
#
# cell_Js = hoja_boot.range(f"S{fila}:U{fila}")
# for i, cell in enumerate(cell_Js):
#     cell.value = round(J_s_boot[i] * 1e-6, 3)
# hoja_boot.update_cells(cell_Js)
#
# cell_Jv = hoja_boot.range(f"W{fila}:Y{fila}")
# for i, cell in enumerate(cell_Jv):
#     cell.value = round(J_v_boot[i], 3)
# hoja_boot.update_cells(cell_Jv)
#
# cell_EH = hoja_boot.range(f"AB{fila}:AD{fila}")
# for i, cell in enumerate(cell_EH):
#     cell.value = round(E_Hall_boot[i] * 1e3, 3)
# hoja_boot.update_cells(cell_EH)
#
#
# # La hoja del ajuste
#
# hoja_fit.update_acell(f"J{fila}", f"{angulo_mva * 180/np.pi:.3g}")
# hoja_fit.update_acell(f"M{fila}", f"{angulo_v_fit * 180/np.pi:.3g}")
# hoja_fit.update_acell(f"N{fila}", f"{angulo_B_fit * 180/np.pi:.3g}")
# hoja_fit.update_acell(f"O{fila}", f"{np.linalg.norm(x_23_fit):.3g}")
# hoja_fit.update_acell(f"P{fila}", f"{np.linalg.norm(x_14_fit):.3g}")
# # hoja_fit.update_acell(f'Q{fila}', f'{:.3g}')
# # hoja_fit.update_acell(f'R{fila}', f'{:.3g}')
#
# hoja_fit.update_acell(f"V{fila}", f"{np.linalg.norm(J_s_fit)*1E-6:.3g}")
# hoja_fit.update_acell(f"Z{fila}", f"{np.linalg.norm(J_v_fit):.3g}")
# hoja_fit.update_acell(f"AA{fila}", f"{np.linalg.norm(fuerza_fit):.3g}")
# hoja_fit.update_acell(f"AE{fila}", f"{np.linalg.norm(E_Hall_fit)*1E3:.3g}")
#
# cell_Js = hoja_fit.range(f"S{fila}:U{fila}")
# for i, cell in enumerate(cell_Js):
#     cell.value = round(J_s_fit[i] * 1e-6, 3)
# hoja_fit.update_cells(cell_Js)
#
# cell_Jv = hoja_fit.range(f"W{fila}:Y{fila}")
# for i, cell in enumerate(cell_Jv):
#     cell.value = round(J_v_fit[i], 3)
# hoja_fit.update_cells(cell_Jv)
#
# cell_EH = hoja_fit.range(f"AB{fila}:AD{fila}")
# for i, cell in enumerate(cell_EH):
#     cell.value = round(E_Hall_fit[i] * 1e3, 3)
# hoja_fit.update_cells(cell_EH)
