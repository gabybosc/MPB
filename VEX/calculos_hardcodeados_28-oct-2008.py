import numpy as np
from _importar_datos import importar_MAG
import sys

sys.path.append("..")
from funciones import Bpara_Bperp, UTC_to_hdec, datenum, donde, SZA, angulo, ancho_mpb

year, month, day, doy = (2008, 10, 28, 302)
t1, t2, t3, t4 = (
    8.541556527954604,
    8.544405851015947,
    8.551476393427427,
    8.556111111111111,
)
tbs = UTC_to_hdec("08:26:40")
tmpr = UTC_to_hdec("08:34:25")
tmva = 8.549442222

ti = t1 - 0.2
tf = t4 + 0.2
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

B_mpb = B[donde(t, t1): donde(t, t4)]
B_mva = B[donde(t, t2): donde(t, tmva)]  # el mva es entre t2 y tmva

# el B medio en la mpb
print(
    "el campo medio en la mpb",
    np.mean(B_mpb, axis=0),
    np.linalg.norm(np.mean(B_mpb, axis=0)),
)

n_fit = np.array([0.587, -0.116, 0.801])
n_mva = np.array([0.517, 0.103, 0.850])
B3_mva = np.dot(B_mpb, n_mva)
B3_fit = np.dot(B_mpb, n_fit)
print("el campo medio a lo largo de la normal del mva", np.mean(B3_mva))
vel_VEX = pos[32] - pos[0]  # estos son los km en 1s (porque son 32 mediciones por seg)
print("norma de vel VEX media", np.linalg.norm(vel_VEX))

print(
    "el cociente entre el campo medio a lo largo de la normal del mva y el campo medio total",
    np.mean(B3_mva) / np.linalg.norm(np.mean(B_mpb, axis=0)),
)

theta_B_mva = angulo(np.mean(B_mpb, axis=0), n_mva) * 180 / np.pi
theta_B_fit = angulo(np.mean(B_mpb, axis=0), n_fit) * 180 / np.pi
print("el ángulo theta_B del mva es ", theta_B_mva)

theta_v_mva = angulo(vel_VEX, n_mva) * 180 / np.pi
theta_v_fit = angulo(vel_VEX, n_fit) * 180 / np.pi
print("el angulo theta_v del mva es ", theta_v_mva)

x14_mva, x23_mva = ancho_mpb(t1, t2, t3, t4, n_mva, np.linalg.norm(vel_VEX))
x14_fit, x23_fit = ancho_mpb(t1, t2, t3, t4, n_fit, np.linalg.norm(vel_VEX))
print("el ancho de la mpb del MVA son (corto y largo)", x23_mva, x14_mva)

print("el campo medio a lo largo de la normal del fit", np.mean(B3_fit))
print(
    "el cociente entre el campo medio a lo largo de la normal del fit y el campo medio total",
    np.mean(B3_fit) / np.linalg.norm(np.mean(B_mpb, axis=0)),
)
print("el ángulo theta_B del fit es ", theta_B_fit)
print("el angulo theta_v del fit es ", theta_v_fit)
print("el ancho de la mpb del MVA son (corto y largo)", x23_fit, x14_fit)

"""
Calcula el giroradio y la long inercial de iones. El giroradio lo calcula usando
v_perp al campo B y también usando la velocidad proyectada en la normal.
"""

# #########CONSTANTES
mp = 1.67e-27  # masa del proton en kg
kB = 1.38e-23  # cte de Boltzmann en J/K
q_e = 1.602e-19  # carga del electrón en C
mu0 = 4 * np.pi * 1e-7  # Tm/A

IMA = np.genfromtxt("dens_20081028-futaana.txt", dtype=str, skip_header=1)

t_IMA = [(UTC_to_hdec(IMA[i, 0]) + UTC_to_hdec(IMA[i, 1])) / 2 for i in range(len(IMA))]

dens = [float(IMA[i, 2]) for i in range(len(IMA))]
vel_protones = [
    np.array([float(IMA[i, 3]), float(IMA[i, 4]), float(IMA[i, 5])])
    for i in range(len(IMA))
]

"""
giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
where vperp is the component of the velocity perpendicular to the
direction of the magnetic field and B is the strength of the magnetic field
"""

# selecciono la región upstream usando un único punto y calculo en ese lugar.
vel_upstream = np.array([-200, -14.2, 1.14])  # coincide bastante bien con lo del clweb
B_punto = (
        np.array([-2.8294, -6.6353, -1.0376]) * 1e-9
)  # B[donde(t, t_IMA[8])] * 1e-9  # para que me de en km

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

print(f"El radio de Larmor es {rg:1.3g} km ")

"""
Proyectando en la normal: pero ya ni me acuerdo por qué hacía esto
"""
# 0.391	-0.129	0.911
v_normal = np.dot(np.mean(vel_upstream, axis=0), n_mva)

# el giroradio entonces:
rg_normal = mp * np.linalg.norm(v_normal) / (q_e * np.linalg.norm(B_punto))

print(
    f"El radio de Larmor con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
)

"""
Longitud inercial
"""
density_mean = 19  # la hardcodeo porque dice 70...
ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km
print(f"La longitud inercial de iones es {ion_length:1.3g} km")
