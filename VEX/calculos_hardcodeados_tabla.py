import numpy as np
from _importar_datos import importar_MAG
import sys

sys.path.append("..")
from funciones import (
    Bpara_Bperp,
    UTC_to_hdec,
    datenum,
    donde,
    SZA,
    angulo,
    ancho_mpb,
    fechas,
    tiempos,
)

year, month, day, doy = fechas()
# ti_sw, tf_sw = tiempos()
# ti_mpr, tf_mpr = tiempos()
ti_ms = UTC_to_hdec("08:27:00")
tf_ms = UTC_to_hdec("08:32:00")
mp = 1.67e-27  # masa del proton en kg
kB = 1.38e-23  # cte de Boltzmann en J/K
q_e = 1.602e-19  # carga del electrón en C
mu0 = 4 * np.pi * 1e-7  # Tm/A

ti = ti_ms - 1
tf = tf_ms + 1
if ti < 0:
    ti = 0
if tf > 24:
    tf = 24

t, B, pos, cl, tpos = importar_MAG(year, doy, ti, tf)
if cl:
    Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)  # si son datos de clweb 1s
else:
    # para datos de PDS filtrados y diezmados
    Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)
Bnorm = np.linalg.norm(B, axis=1)

# para el 2008 es futaana, para 2014 es clweb
IMA = np.genfromtxt(f"dens_{year}{month}{day}-futaana.txt", dtype=str, skip_header=1)
t_dens = np.array(
    [(UTC_to_hdec(IMA[i, 0]) + UTC_to_hdec(IMA[i, 1])) / 2 for i in range(len(IMA))]
)
dens = np.array([float(IMA[i, 2]) for i in range(len(IMA))])
vel = np.array(
    [[float(IMA[i, 3]), float(IMA[i, 4]), float(IMA[i, 5])] for i in range(len(IMA))]
)

# IMA = np.loadtxt(
#     f"densidad_{year}{month}{day}-clweb.txt", skiprows=1
# )  # estos son los del clweb
# t_dens = np.array(IMA[:, 0] + IMA[:, 1] / 60 + IMA[:, 2] / 3600)  # hdec
# dens = IMA[:, -2]
# vel = IMA[:, -1]

"""
Longitud inercial
"""
density_ms = dens[
    donde(t_dens, ti_ms) : donde(t_dens, tf_ms)
]  # la hardcodeo porque dice 70...
ion_length = 2.28e07 / np.sqrt(np.mean(density_ms)) * 1e-5  # km
print(f"La longitud inercial de iones es {ion_length:1.3g} km")

"""
giroradio: no puedo calcularlo para 2014
"""
vel_upstream = np.mean(
    vel[donde(t_dens, ti_ms) : donde(t_dens, tf_ms), :], axis=0
)  # coincide bastante bien con lo del clweb
B_punto = (
    np.mean(B[donde(t, ti_ms) : donde(t, tf_ms)], axis=0) * 1e-9
)  # para que me de en km

dot_product = np.dot(vel_upstream, B_punto)

# Calcular el módulo del campo magnético
B_norm_squared = np.dot(B_punto, B_punto)

# Calcular la componente paralela de la velocidad
v_parallel = (dot_product / B_norm_squared) * B_punto

# Calcular la componente perpendicular de la velocidad
v_perp = vel_upstream - v_parallel

# el giroradio entonces:
gf = (q_e * np.linalg.norm(B_punto)) / mp
rg = mp * np.linalg.norm(v_perp) / (q_e * np.linalg.norm(B_punto))

print("El radio de Larmor es", rg)

# """
# presiones
# """
# dens_sw = dens[donde(t_dens, ti_sw) : donde(t_dens, tf_sw)]
# vel_sw = vel[donde(t_dens, ti_sw) : donde(t_dens, tf_sw)]
#
# # p_dyn = (
# #         mp * (np.mean(dens_sw) * 1e6) * (np.mean(np.linalg.norm(vel_sw, axis=0) * 1e3)) ** 2
# # )  # p dyn en el SW futaana
#
# p_dyn = (
#     mp * (np.mean(dens_sw) * 1e6) * (np.mean(vel_sw * 1e3)) ** 2
# )  # p dyn en el SW de clweb
#
# p_mag = (
#     (np.mean(np.linalg.norm(B[donde(t, ti_mpr) : donde(t, tf_mpr), :], axis=1) * 1e-9))
#     ** 2
#     / 2
#     / mu0
# )  # pmag en la MPR
#
# print(f"La presión dinámica en el SW es {p_dyn * 1e9:1.3g} nPa")
# print(f"La presión magnética en la MPR es {p_mag * 1e9:1.3g} nPa")
