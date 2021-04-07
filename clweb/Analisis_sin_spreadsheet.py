import numpy as np
import matplotlib.pyplot as plt
import sys
from MVA_sin_spreadsheet import MVA, ajuste, normal_coplanar
from importar_datos import importar_mag, importar_lpw
import datetime as dt
import cdflib as cdf

sys.path.append("..")

from funciones import (
    angulo,
    ancho_mpb,
    corrientes,
    diezmar,
    donde,
    find_nearest,
    find_nearest_final,
    find_nearest_inicial,
    UTC_to_hdec,
    unix_to_decimal,
    fechas,
    tiempos,
)
from funciones_metodos import bootstrap

"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
"""


def importar_swia(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    # if gethostname() == "magneto2":
    #     path = f"../../../../media/gabybosc/datos/SWIA/"
    # elif gethostname() == "gabybosc":
    #     path = "../../datos/SWIA/"
    # else:
    path = f"../../../datos/SWIA/"

    swia = cdf.CDF(path + f"mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf")

    t_unix = swia.varget("time_unix")
    density = swia.varget("density")  # cm⁻³
    temperature = swia.varget("temperature_mso")  # eV
    vel_mso_xyz = swia.varget("velocity_mso")  # km/s

    t_swia = unix_to_decimal(t_unix)
    inicio = donde(t_swia, ti)
    fin = donde(t_swia, tf)

    t_cut = t_swia[inicio:fin]
    density_cut = density[inicio:fin]
    temperature_cut = temperature[inicio:fin]
    vel_mso_cut = vel_mso_xyz[inicio:fin]  # km/s

    return swia, t_cut, density_cut, temperature_cut, vel_mso_cut


year, month, day, doy = fechas()
# ti, tf = tiempos("Región de análisis (no MVA)")
# ti_MVA, tf_MVA = tiempos("Intervalo del MVA")

# year, month, day, doy = 2016, "03", 16, 76
# ti_MVA, tf_MVA = UTC_to_hdec("18:13:33"), UTC_to_hdec("18:14:06")
# ti, tf = UTC_to_hdec("17:55:00"), UTC_to_hdec("18:30:00")

# year, month, day, doy = 2015, 10, 10, #doy
# ti_MVA, tf_MVA = 12.675,12.68444444
# ti, tf = 12.4, 12.9
# year, month, day, doy = 2015, 10, 12, 285
ti_MVA, tf_MVA = 19.31666667, 19.32388889
ti, tf = 19.15, 19.45
# year, month, day, doy = 2016, 04, 05, #doy
# ti_MVA, tf_MVA = 5.271388889,5.278611111
# ti, tf = 5.15, 5.45
# year, month, day, doy = 2016, 03, 31, #doy
# ti_MVA, tf_MVA = 13.07388889,13.08055556
# ti, tf = 12.45, 13.15

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)
swia, t_swia, density, temperature, vel_mso = importar_swia(
    year, month, day, ti, tf
)
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

inicio_up = np.where(t == find_nearest_inicial(t, t1 - 0.015))[0][0]
fin_up = np.where(t == find_nearest_final(t, t1))[0][0]
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = np.where(t == find_nearest_inicial(t, t4))[0][0]
fin_down = np.where(t == find_nearest_final(t, t4 + 0.015))[0][0]
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


# plt.show()

idx = diezmar(t, t_swia)

tmag_diezmado = t[idx]
B_cut = B[idx]
posicion_cut = posicion[idx]

#magnetofunda:

ti_funda = donde(tmag_diezmado, t1-0.16)
tf_funda = donde(tmag_diezmado, t1)
####################


B_avg = np.empty((len(idx), 3))
v_maven = np.empty((len(idx), 3))
v_planet = np.empty((len(idx), 3))


for i in range(len(idx)):
    B_avg[i, :] = np.mean(B_cut[i : i + 30, :], axis=0)
    v_maven[i, :] = (
        (posicion_cut[ti_funda + 1, :] - posicion_cut[ti_funda, :])
        / (tmag_diezmado[ti_funda + 1] - tmag_diezmado[ti_funda])
        / 3600
    )  # en km/s
for i in range(len(idx)):
    v_planet[i, :] = np.nanmean(
        vel_mso[i : i + 30, :] + v_maven[i : i + 30, :], axis=0
    )  # km/s

E_convective = np.cross(-v_planet * 1e3, B_avg * 1e-9) * 1e3  # mV/m
E_convective_norm = np.linalg.norm(E_convective, axis=1)  # mV/m

print(
    f"El Ecv medio en la magnetofunda es {np.mean(E_convective[ti_funda:tf_funda] * 1e3, axis=0)} mV/m"
)
print(
    f"El |Ecv| medio en la magnetofunda es {np.mean(E_convective_norm[ti_funda:tf_funda]):1.3g} mV/m"
    )  # en mV/m
print(f"el campo de hall en la MPB es {E_Hall * 1e3} mV/m, {np.linalg.norm(E_Hall * 1e3):1.3g} mV/m")
print(f"Eh / Ecv = {np.linalg.norm(E_Hall * 1e3)/np.mean(E_convective_norm[ti_funda:tf_funda])}")
