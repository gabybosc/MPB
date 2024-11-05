import numpy as np
import matplotlib.pyplot as plt
from MVA_hires import MVA, ajuste, bootstrap_completo
import sys
from importar_datos import importar_mag, importar_lpw, importar_fila

sys.path.append("..")

from funciones import (
    angulo,
    ancho_mpb,
    corrientes,
    find_nearest_final,
    find_nearest_inicial,
    fechas,
    donde,
    tiempos,
)

"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
"""

year, month, day, doy = fechas()
ti_MVA, tf_MVA = tiempos("Intervalo del MVA")
ti, tf = ti_MVA - 0.5, tf_MVA + 0.5

print(
    "Si tira error de que no encuentra el path, hay que abrir una vez el disco manualmente para que lo monte"
)

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
# lpw, t_lpw, e_density, flag = importar_lpw(year, month, day, ti, tf)
x3, B_cut, t_cut, posicion_cut, nr = MVA(year, month, day, ti_MVA, tf_MVA)
normal_ajuste, t1, t2, t3, t4 = ajuste(year, month, day, doy, ti_MVA, tf_MVA, nr)
#
normal_boot = bootstrap_completo(B_cut, len(t_cut), nr, 1000)

######
# constantes
q_e = 1.6e-19  # carga electron #C

#########
# buscamos el ángulo entre las normales
angulo_mva = angulo(normal_ajuste, x3)
angulo_boot = angulo(normal_boot, x3)


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

angulo_v_mva = angulo(x3, v_media)
angulo_v_fit = angulo(normal_ajuste, v_media)
angulo_v_boot = angulo(normal_boot, v_media)

B_medio_vectorial = np.mean(B_cut, axis=0)
angulo_B_mva = angulo(x3, B_medio_vectorial)
angulo_B_fit = angulo(normal_ajuste, B_medio_vectorial)
angulo_B_boot = angulo(normal_boot, B_medio_vectorial)

######
# Espesor de la MPB proyectado sobre las distintas normales

x_14_fit, x_23_fit = ancho_mpb(t1, t2, t3, t4, normal_ajuste, v_media)
x_14_MVA, x_23_MVA = ancho_mpb(t1, t2, t3, t4, x3, v_media)
x_14_boot, x_23_boot = ancho_mpb(t1, t2, t3, t4, normal_boot, v_media)

#########
# Análisis de corrientes

B_upstream = corte(t, t1 - 0.015, t1, B)
B_downstream = corte(t, t4, t4 + 0.015, B)

omega = angulo(B_upstream, B_downstream)

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, x_23_MVA)
J_s_fit, J_v_fit = corrientes(normal_ajuste, B_upstream, B_downstream, x_23_fit)
J_s_boot, J_v_boot = corrientes(normal_boot, B_upstream, B_downstream, x_23_boot)

inicio_down = donde(t, t1 - 0.015)
fuerza_mva = fuerza(J_v_MVA, B[inicio_down, :])
fuerza_fit = fuerza(J_v_fit, B[inicio_down, :])
fuerza_boot = fuerza(J_v_boot, B[inicio_down, :])

# e_density = lpw[:, -1]
# t_lpw = lpw[:, 3] + lpw[:, 4] / 60 + lpw[:, 5] / 3600

# n_e = corte(t_lpw, t2, t3, e_density)
# n_e = n_e * 1e6  # m⁻³
# if np.isnan(n_e):
#     n_e = 1e7
#     print("LPW no tiene datos de densidad, asumí n_e = 1E7")
# q_e = 1.6e-19  # carga electron #C

# E_Hall = fuerza_mva / (q_e * n_e)  # V/m
# E_Hall_fit = fuerza_fit / (q_e * n_e)  # V/m
# E_Hall_boot = fuerza_boot / (q_e * n_e)  # V/m

nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(
    year, month, day, int(ti_MVA)
)

date_entry = f"{year}-{month}-{day}"

##########
# Parámetros

hoja_parametros.update_acell(f"Z{nr}", f"{omega * 180 / np.pi:.3g}")
hoja_parametros.update_acell(f"S{nr}", f"{np.linalg.norm(v_media):.3g}")

cell_vel = hoja_parametros.range(f"P{nr}:R{nr}")
for i, cell in enumerate(cell_vel):
    cell.value = v_media[i]
hoja_parametros.update_cells(cell_vel)

cell_Bup = hoja_parametros.range(f"T{nr}:V{nr}")
for i, cell in enumerate(cell_Bup):
    cell.value = B_upstream[i]
hoja_parametros.update_cells(cell_Bup)

cell_Bdown = hoja_parametros.range(f"W{nr}:Y{nr}")
for i, cell in enumerate(cell_Bdown):
    cell.value = B_downstream[i]
hoja_parametros.update_cells(cell_Bdown)

# La hoja del MVA

hoja_mva.update_acell(f"W{nr}", f"{angulo_v_mva * 180 / np.pi:.3g}")
hoja_mva.update_acell(f"X{nr}", f"{angulo_B_mva * 180 / np.pi:.3g}")
hoja_mva.update_acell(f"Y{nr}", f"{np.linalg.norm(x_23_MVA):.3g}")
hoja_mva.update_acell(f"Z{nr}", f"{np.linalg.norm(x_14_MVA):.3g}")
# hoja_mva.update_acell(f'AA{nr}', f'{:.3g}')
# hoja_mva.update_acell(f'AB{nr}', f'{:.3g}')

hoja_mva.update_acell(f"AF{nr}", f"{np.linalg.norm(J_s_MVA) * 1E-6:.3g}")
hoja_mva.update_acell(f"AJ{nr}", f"{np.linalg.norm(J_v_MVA):.3g}")
hoja_mva.update_acell(f"AK{nr}", f"{np.linalg.norm(fuerza_mva):.3g}")
# hoja_mva.update_acell(f"AO{nr}", f"{np.linalg.norm(E_Hall) * 1E3:.3g}")  # mV/m

cell_Js = hoja_mva.range(f"AC{nr}:AE{nr}")
for i, cell in enumerate(cell_Js):
    cell.value = J_s_MVA[i] * 1e-6
hoja_mva.update_cells(cell_Js)

cell_Jv = hoja_mva.range(f"AG{nr}:AI{nr}")
for i, cell in enumerate(cell_Jv):
    cell.value = J_v_MVA[i]
hoja_mva.update_cells(cell_Jv)

cell_EH = hoja_mva.range(f"AL{nr}:AN{nr}")
# for i, cell in enumerate(cell_EH):
#     cell.value = E_Hall[i] * 1e3
hoja_mva.update_cells(cell_EH)

# La hoja del bootstrap

hoja_boot.update_acell(f"M{nr}", f"{angulo_v_boot * 180 / np.pi:.3g}")
hoja_boot.update_acell(f"N{nr}", f"{angulo_B_boot * 180 / np.pi:.3g}")
hoja_boot.update_acell(f"O{nr}", f"{np.linalg.norm(x_23_boot):.3g}")
hoja_boot.update_acell(f"P{nr}", f"{np.linalg.norm(x_14_boot):.3g}")
# hoja_boot.update_acell(f'Q{nr}', f'{:.3g}')
# hoja_boot.update_acell(f'R{nr}', f'{:.3g}')

hoja_boot.update_acell(f"V{nr}", f"{np.linalg.norm(J_s_boot) * 1E-6:.3g}")
hoja_boot.update_acell(f"Z{nr}", f"{np.linalg.norm(J_v_boot):.3g}")
hoja_boot.update_acell(f"AA{nr}", f"{np.linalg.norm(fuerza_boot):.3g}")
# hoja_boot.update_acell(f"AE{nr}", f"{np.linalg.norm(E_Hall_boot) * 1E3:.3g}")

cell_Js = hoja_boot.range(f"S{nr}:U{nr}")
for i, cell in enumerate(cell_Js):
    cell.value = J_s_boot[i] * 1e-6
hoja_boot.update_cells(cell_Js)

cell_Jv = hoja_boot.range(f"W{nr}:Y{nr}")
for i, cell in enumerate(cell_Jv):
    cell.value = J_v_boot[i]
hoja_boot.update_cells(cell_Jv)

cell_EH = hoja_boot.range(f"AB{nr}:AD{nr}")
# for i, cell in enumerate(cell_EH):
# cell.value = E_Hall_boot[i] * 1e3
hoja_boot.update_cells(cell_EH)

# La hoja del ajuste

hoja_fit.update_acell(f"J{nr}", f"{angulo_mva * 180 / np.pi:.3g}")
hoja_fit.update_acell(f"M{nr}", f"{angulo_v_fit * 180 / np.pi:.3g}")
hoja_fit.update_acell(f"N{nr}", f"{angulo_B_fit * 180 / np.pi:.3g}")
hoja_fit.update_acell(f"O{nr}", f"{x_23_fit:.3g}")
hoja_fit.update_acell(f"P{nr}", f"{x_14_fit:.3g}")
# hoja_fit.update_acell(f'Q{nr}', f'{:.3g}')
# hoja_fit.update_acell(f'R{nr}', f'{:.3g}')

hoja_fit.update_acell(f"V{nr}", f"{np.linalg.norm(J_s_fit) * 1E-6:.3g}")
hoja_fit.update_acell(f"Z{nr}", f"{np.linalg.norm(J_v_fit):.3g}")
hoja_fit.update_acell(f"AA{nr}", f"{np.linalg.norm(fuerza_fit):.3g}")
# hoja_fit.update_acell(f"AE{nr}", f"{np.linalg.norm(E_Hall_fit) * 1E3:.3g}")

cell_Js = hoja_fit.range(f"S{nr}:U{nr}")
for i, cell in enumerate(cell_Js):
    cell.value = J_s_fit[i] * 1e-6
hoja_fit.update_cells(cell_Js)

cell_Jv = hoja_fit.range(f"W{nr}:Y{nr}")
for i, cell in enumerate(cell_Jv):
    cell.value = J_v_fit[i]
hoja_fit.update_cells(cell_Jv)

cell_EH = hoja_fit.range(f"AB{nr}:AD{nr}")
# for i, cell in enumerate(cell_EH):
# cell.value = E_Hall_fit[i] * 1e3
hoja_fit.update_cells(cell_EH)
print("Completó.")

plt.show()
