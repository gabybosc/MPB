import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from importar_datos import importar_mag, importar_gdocs

import sys

sys.path.append("..")
from funciones import fechas, tiempos, donde


"""
Calcula la velocidad media de los protones del viento solar de SWICA y SWIFA
en distintas regiones proyectada sobre la normal y luego calcula la ram pressure.
La compara con la presión magnética.
"""

np.set_printoptions(precision=4)


def importar_swica(year, month, day, ti, tf):

    path = f"../../../datos/Ecinetica/{year}-{month}-{day}/"

    swica = np.loadtxt(path + "SWICA.asc")

    t = swica[:, 3] + swica[:, 4] / 60 + swica[:, 5] / 3600  # hdec

    density = swica[:, 6]

    vel_norm = swica[:, -1]
    vel_mso_xyz = swica[:, 7:10]

    return swica, t, density, vel_mso_xyz, vel_norm


def importar_swifa(year, month, day, ti, tf):

    path = f"../../../datos/Ecinetica/{year}-{month}-{day}/"

    swifa = np.loadtxt(path + "SWIFA.asc")

    t = swifa[:, 3] + swifa[:, 4] / 60 + swifa[:, 5] / 3600  # hdec

    density = swifa[:, 6]

    vel_norm = swifa[:, -1]
    vel_mso_xyz = swifa[:, 7:10]

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    vel_cut = vel_mso_xyz[inicio:fin]
    vel_norm_cut = vel_norm[inicio:fin]
    d_cut = density[inicio:fin]

    return swifa, t_cut, d_cut, vel_cut, vel_norm_cut


mp = 1.67e-27  # kg
mu_0 = np.pi * 4e-7  # T m / A
year, month, day, doy = fechas()
ti_sw, tf_sw = tiempos("tiempo inicial y final en el viento solar\n")
# ti_MS, tf_MS = tiempos("tiempo inicial y final de la magnetofunda\n")

hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()

fila = input("En qué fila del gdoc estamos?\n")
tf_MPR = float(hoja_parametros.cell(fila, 6).value)
ti_MPR = tf_MPR - 0.015
normal = [
    float(hoja_MVA.cell(fila, 16).value),
    float(hoja_MVA.cell(fila, 17).value),
    float(hoja_MVA.cell(fila, 18).value),
]

swifa, t_sw, density_sw, vel_sw, vel_norm_sw = importar_swifa(
    year, month, day, ti_sw, tf_sw
)
# swica, t_swica, density_swica, vel_swica, vel_norm_swica = importar_swica(
#     year, month, day, tf_sw, tf_MPR
# )

mag, t_mag, B_mag, posicion = importar_mag(year, month, day, ti_MPR, tf_sw)
B_norm = np.linalg.norm(B_mag, axis=1) * 1e-9  # Tesla


""" Solar wind """
v_parallel_sw = np.dot(vel_sw, normal)

v_mean_sw = np.mean(v_parallel_sw)  # km/s
density_sw_mean = np.mean(density_sw)  # 1/cm3
error_v_sw = np.std(v_parallel_sw)
error_density_sw = np.std(density_sw)
# print(
#     f"La velocidad media en el viento solar en dirección de la normal es {v_mean_sw:.3g} km/s"
# )

v_sw = ufloat(v_mean_sw, error_v_sw)
n_sw = ufloat(density_sw_mean, error_density_sw)
p_ram = mp * (v_sw * 1e3) ** 2 * (n_sw * 1e6)  # J/m3


print(
    f"La ram pressure en el SW a partir de los parámetros promedio es {p_ram*1e9:.3g} nPa"
)

# inicio_sw = donde(t_mag, ti_sw)
# fin_sw = donde(t_mag, tf_sw)
# B_mean_sw = np.mean(B_norm[inicio_sw : fin_sw + 1])
# p_mag = B_mean_sw ** 2 / (2 * mu_0)  # J/m3
#
# print(f"El B medio en el SW a partir de los parámetros promedio es {B_mean_sw:.3g} T")
#
# print(
#     f"La presión magnética en el SW a partir de los parámetros promedio es {p_mag*1e9:.3g} nPa"
# )

#
# """ Magnetosheath """
# inicio_ms = donde(t_swica, ti_MS)
# fin_ms = donde(t_swica, tf_MS)
#
# v_parallel_ms = np.dot(vel_swica[inicio_ms : fin_ms + 1], normal)
#
# v_mean_ms = np.mean(v_parallel_ms)  # km/s
# density_ms_mean = np.mean(density_swica[inicio_ms : fin_ms + 1])  # 1/cm3
# print(
#     f"La velocidad media en la magnetofunda en dirección de la normal es {v_mean_ms:.3g} km/s"
# )
#
# p_ram_ms = mp * (v_mean_ms * 1e3) ** 2 * (density_ms_mean * 1e6)  # J/m3
#
# print(
#     f"La ram pressure en la magnetofunda a partir de los parámetros promedio es {p_ram_ms:.3g} J/m3"
# )
#
# inicio_ms_mag = donde(t_mag, ti_MS)
# fin_ms_mag = donde(t_mag, tf_MS)
# B_mean_ms = np.mean(B_norm[inicio_ms_mag : fin_ms_mag + 1])
# p_mag_ms = B_mean_ms ** 2 / (2 * mu_0)  # J/m3
#
# print(
#     f"El B medio en la magnetofunda a partir de los parámetros promedio es {B_mean_ms:.3g} T"
# )
#
# print(
#     f"La presión magnética en la magnetofunda a partir de los parámetros promedio es {p_mag_ms:.3g} J/m3"
# )
#
#
# """ Upstream """
# inicio_up = donde(t_swica, tf_MS - 0.02)
# v_parallel_up = np.dot(vel_swica[inicio_up : fin_ms + 1], normal)
#
# v_mean_up = np.mean(v_parallel_up)  # km/s
# density_up_mean = np.mean(density_swica[inicio_up : fin_ms + 1])  # 1/cm3
# print(
#     f"La velocidad media en la región upstream en dirección de la normal es {v_mean_up:.3g} km/s"
# )
#
# p_ram_up = mp * (v_mean_up * 1e3) ** 2 * (density_up_mean * 1e6)  # J/m3
#
# print(
#     f"La ram pressure en la región upstream a partir de los parámetros promedio es {p_ram_up:.3g} J/m3"
# )
#
# inicio_up_mag = donde(t_mag, tf_MS - 0.02)
# B_mean_up = np.mean(B_norm[inicio_up_mag : fin_ms_mag + 1])
# p_mag_up = B_mean_up ** 2 / (2 * mu_0)  # J/m3
#
# print(
#     f"El B medio en la región upstream a partir de los parámetros promedio es {B_mean_up:.3g} T"
# )
#
# print(
#     f"La presión magnética en la región upstream a partir de los parámetros promedio es {p_mag_up:.3g} J/m3"
# )


""" MPR (ojo que acá ya swia no funciona bien ) """
# inicio_mpr = donde(t_swica, ti_MPR)
# fin_mpr = donde(t_swica, ti_MPR + 0.02)
# v_parallel_mpr = np.dot(vel_swica[inicio_mpr : fin_mpr + 1], normal)
#
# v_mean_mpr = np.mean(v_parallel_mpr)  # km/s
# density_mpr_mean = np.mean(density_swica[inicio_mpr : fin_mpr + 1])  # 1/cm3
# print(
#     f"La velocidad media en la región downstream (ojo que acá SWIA no funciona bien) en dirección de la normal es {v_mean_mpr:.3g} km/s"
# )
#
# p_ram_mpr = mp * (v_mean_mpr * 1e3) ** 2 * (density_mpr_mean * 1e6)  # J/m3
#
# print(
#     f"La ram pressure en la región downstream a partir de los parámetros promedio es {p_ram_mpr*1e9:.3g} nPa"
# )
inicio_mpr_mag = donde(t_mag, ti_MPR)
fin_mpr_mag = donde(t_mag, tf_MPR)
B_mean_mpr = np.mean(B_norm[inicio_mpr_mag : fin_mpr_mag + 1])
error_B = np.std(B_norm[inicio_mpr_mag:fin_mpr_mag])
B_MPR = ufloat(B_mean_mpr, error_B)

p_mag_mpr = B_MPR ** 2 / (2 * mu_0)  # J/m3

#
# print(
#     f"El B medio en la región downstream a partir de los parámetros promedio es {B_mean_mpr:.3g} T"
# )

print(
    f"La presión magnética en la región downstream a partir de los parámetros promedio es {p_mag_mpr*1e9:.3g} nPa"
)
