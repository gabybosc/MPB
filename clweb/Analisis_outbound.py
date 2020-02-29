import numpy as np
import matplotlib.pyplot as plt
from MVA_hires import MVA, ajuste, acceso_spreadsheet, bootstrap_completo
import sys
from importar_datos import importar_mag, importar_lpw

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
)


"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
No estoy segura de si debería cambiar algo más que simplemente el orden como tomo t1 t2 t3 t4
"""


year, month, day, doy = fechas()
ti, tf = tiempos("Región de análisis (no MVA)")
ti_MVA, tf_MVA = tiempos("Intervalo del MVA")

print(
    "Si tira error de que no encuentra el path, hay que abrir una vez el disco manualmente para que lo monte"
)

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)
x3, B_cut, t_cut, posicion_cut, nr = MVA(year, month, day, ti_MVA, tf_MVA)
normal_ajuste, t4, t3, t2, t1 = ajuste(year, month, day, doy, ti_MVA, tf_MVA, nr)

M = len(t)
M_cut = len(t_cut)

normal_boot = bootstrap_completo(B_cut, M_cut, nr, 1000)


######
# constantes
q_e = 1.6e-19  # carga electron #C

#########
# buscamos el ángulo entre las normales
angulo_mva = angulo(normal_ajuste, x3)
angulo_boot = angulo(normal_boot, x3)

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
# Espesor de la MPB proyectado sobre las distintas normales

x_14_fit, x_23_fit = ancho_mpb(t1, t2, t3, t4, normal_ajuste, v_media)
x_14_MVA, x_23_MVA = ancho_mpb(t1, t2, t3, t4, x3, v_media)
x_14_boot, x_23_boot = ancho_mpb(t1, t2, t3, t4, normal_boot, v_media)


#########

# Análisis de corrientes

inicio_up = np.where(t == find_nearest_inicial(t, t1 - 0.015))[0][0]
fin_up = np.where(t == find_nearest_final(t, t1))[0][0]
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = np.where(t == find_nearest_inicial(t, t4))[0][0]
fin_down = np.where(t == find_nearest_final(t, t4 + 0.015))[0][0]
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

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

E_Hall = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
E_Hall_fit = np.cross(J_v_fit * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
E_Hall_boot = np.cross(J_v_boot * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m


(
    hoja_parametros,
    hoja_mva,
    hoja_boot,
    hoja_fit,
    fecha_sheet,
    hora_sheet,
) = acceso_spreadsheet()

date_entry = f"{year}-{month}-{day}"

##########
# Parámetros

hoja_parametros.update_acell(f"Z{nr}", f"{omega*180/np.pi:.3g}")
hoja_parametros.update_acell(f"S{nr}", f"{np.linalg.norm(v_media):.3g}")

cell_vel = hoja_parametros.range(f"P{nr}:R{nr}")
for i, cell in enumerate(cell_vel):
    cell.value = round(v_media[i], 2)
hoja_parametros.update_cells(cell_vel)

cell_Bup = hoja_parametros.range(f"T{nr}:V{nr}")
for i, cell in enumerate(cell_Bup):
    cell.value = round(B_upstream[i], 2)
hoja_parametros.update_cells(cell_Bup)

cell_Bdown = hoja_parametros.range(f"W{nr}:Y{nr}")
for i, cell in enumerate(cell_Bdown):
    cell.value = round(B_downstream[i], 2)
hoja_parametros.update_cells(cell_Bdown)


# La hoja del MVA

hoja_mva.update_acell(f"W{nr}", f"{angulo_v_mva * 180/np.pi:.3g}")
hoja_mva.update_acell(f"X{nr}", f"{angulo_B_mva * 180/np.pi:.3g}")
hoja_mva.update_acell(f"Y{nr}", f"{np.linalg.norm(x_23_MVA):.3g}")
hoja_mva.update_acell(f"Z{nr}", f"{np.linalg.norm(x_14_MVA):.3g}")
# hoja_mva.update_acell(f'AA{nr}', f'{:.3g}')
# hoja_mva.update_acell(f'AB{nr}', f'{:.3g}')

hoja_mva.update_acell(f"AF{nr}", f"{np.linalg.norm(J_s_MVA)*1E-6:.3g}")
hoja_mva.update_acell(f"AJ{nr}", f"{np.linalg.norm(J_v_MVA):.3g}")
hoja_mva.update_acell(f"AK{nr}", f"{np.linalg.norm(fuerza_mva):.3g}")
hoja_mva.update_acell(f"AO{nr}", f"{np.linalg.norm(E_Hall)*1E3:.3g}")  # mV/m

cell_Js = hoja_mva.range(f"AC{nr}:AE{nr}")
for i, cell in enumerate(cell_Js):
    cell.value = round(J_s_MVA[i] * 1e-6, 3)
hoja_mva.update_cells(cell_Js)

cell_Jv = hoja_mva.range(f"AG{nr}:AI{nr}")
for i, cell in enumerate(cell_Jv):
    cell.value = round(J_v_MVA[i], 3)
hoja_mva.update_cells(cell_Jv)

cell_EH = hoja_mva.range(f"AL{nr}:AN{nr}")
for i, cell in enumerate(cell_EH):
    cell.value = round(E_Hall[i] * 1e3, 3)
hoja_mva.update_cells(cell_EH)


# La hoja del bootstrap

hoja_boot.update_acell(f"M{nr}", f"{angulo_v_boot * 180/np.pi:.3g}")
hoja_boot.update_acell(f"N{nr}", f"{angulo_B_boot * 180/np.pi:.3g}")
hoja_boot.update_acell(f"O{nr}", f"{np.linalg.norm(x_23_boot):.3g}")
hoja_boot.update_acell(f"P{nr}", f"{np.linalg.norm(x_14_boot):.3g}")
# hoja_boot.update_acell(f'Q{nr}', f'{:.3g}')
# hoja_boot.update_acell(f'R{nr}', f'{:.3g}')

hoja_boot.update_acell(f"V{nr}", f"{np.linalg.norm(J_s_boot)*1E-6:.3g}")
hoja_boot.update_acell(f"Z{nr}", f"{np.linalg.norm(J_v_boot):.3g}")
hoja_boot.update_acell(f"AA{nr}", f"{np.linalg.norm(fuerza_boot):.3g}")
hoja_boot.update_acell(f"AE{nr}", f"{np.linalg.norm(E_Hall_boot)*1E3:.3g}")

cell_Js = hoja_boot.range(f"S{nr}:U{nr}")
for i, cell in enumerate(cell_Js):
    cell.value = round(J_s_boot[i] * 1e-6, 3)
hoja_boot.update_cells(cell_Js)

cell_Jv = hoja_boot.range(f"W{nr}:Y{nr}")
for i, cell in enumerate(cell_Jv):
    cell.value = round(J_v_boot[i], 3)
hoja_boot.update_cells(cell_Jv)

cell_EH = hoja_boot.range(f"AB{nr}:AD{nr}")
for i, cell in enumerate(cell_EH):
    cell.value = round(E_Hall_boot[i] * 1e3, 3)
hoja_boot.update_cells(cell_EH)


# La hoja del ajuste

hoja_fit.update_acell(f"J{nr}", f"{angulo_mva * 180/np.pi:.3g}")
hoja_fit.update_acell(f"M{nr}", f"{angulo_v_fit * 180/np.pi:.3g}")
hoja_fit.update_acell(f"N{nr}", f"{angulo_B_fit * 180/np.pi:.3g}")
hoja_fit.update_acell(f"O{nr}", f"{x_23_fit:.3g}")
hoja_fit.update_acell(f"P{nr}", f"{x_14_fit:.3g}")
# hoja_fit.update_acell(f'Q{nr}', f'{:.3g}')
# hoja_fit.update_acell(f'R{nr}', f'{:.3g}')

hoja_fit.update_acell(f"V{nr}", f"{np.linalg.norm(J_s_fit)*1E-6:.3g}")
hoja_fit.update_acell(f"Z{nr}", f"{np.linalg.norm(J_v_fit):.3g}")
hoja_fit.update_acell(f"AA{nr}", f"{np.linalg.norm(fuerza_fit):.3g}")
hoja_fit.update_acell(f"AE{nr}", f"{np.linalg.norm(E_Hall_fit)*1E3:.3g}")

cell_Js = hoja_fit.range(f"S{nr}:U{nr}")
for i, cell in enumerate(cell_Js):
    cell.value = round(J_s_fit[i] * 1e-6, 3)
hoja_fit.update_cells(cell_Js)

cell_Jv = hoja_fit.range(f"W{nr}:Y{nr}")
for i, cell in enumerate(cell_Jv):
    cell.value = round(J_v_fit[i], 3)
hoja_fit.update_cells(cell_Jv)

cell_EH = hoja_fit.range(f"AB{nr}:AD{nr}")
for i, cell in enumerate(cell_EH):
    cell.value = round(E_Hall_fit[i] * 1e3, 3)
hoja_fit.update_cells(cell_EH)
print("Completó.")

plt.show()
