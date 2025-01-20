import numpy as np
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
    corrientes,
)

path = "../../../datos/Titan/t96_tswis_1s.ascii"
datos = np.loadtxt(path)
tiempo = datos[:, 0]
i = donde(tiempo, 24)
f = donde(tiempo, 25)
t1, t2, t3, t4 = 24.54175182, 24.55058123, 24.57651763, 24.58203602

t, B, Bnorm, posicion = datos[i:f, 0], datos[i:f, 1:4], datos[i:f, 4], datos[i:f, 5:8]
B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, t[0] + 0.2, t[-1] - 0.2)

tiempo_mag = np.array([np.datetime64(datenum(2013, 11, 30, x)) for x in t])
tiempo_mag_delta = np.array([np.datetime64(datenum(2013, 11, 30, x)) for x in t_plot])
tm1 = donde(t, t1)
tm2 = donde(t, t2)
tm3 = donde(t, t3)
tm4 = donde(t, t4)
tbs = 24 + UTC_to_hdec("00:24:20")
tmpr = 24 + UTC_to_hdec("00:40:30")
tmva_i = UTC_to_hdec("24:33:04")
tmva_f = UTC_to_hdec("24:33:45")

B_mpb = B[donde(t, t1): donde(t, t4)]
B_mva = B[donde(t, tmva_i): donde(t, tmva_f)]

# el B medio en la mpb
print(
    "el campo medio en la mpb",
    np.mean(B_mpb, axis=0),
    np.linalg.norm(np.mean(B_mpb, axis=0)),
)

n_mva = np.array([-0.874, -0.086, 0.478])
B3_mva = np.dot(B_mpb, n_mva)
print("el campo medio a lo largo de la normal del mva", np.mean(B3_mva))


def velocidad(posicion_cut, tpos):
    M = len(posicion_cut)
    v_punto = np.zeros((M - 1, 3))
    deltat = np.zeros(M - 1)
    for i in range(len(v_punto)):
        deltat[i] = (tpos[i + 1] - tpos[i]) * 3600  # delta t en segundos
        v_punto[i] = (posicion_cut[i + 1, :] - posicion_cut[i, :]) / deltat[i]
        # en km/s
    # la velocidad promedio
    v_media = np.mean(v_punto, axis=0)
    return v_media


vel_cassini = velocidad(posicion[tm1:tm4], t[tm1:tm4])
print("norma de vel cassini media", np.linalg.norm(vel_cassini))

print(
    "el cociente entre el campo medio a lo largo de la normal del mva y el campo medio total",
    np.mean(B3_mva) / np.linalg.norm(np.mean(B_mpb, axis=0)),
)

theta_B_mva = angulo(np.mean(B_mpb, axis=0), n_mva) * 180 / np.pi
print("el ángulo theta_B del mva es ", theta_B_mva)

theta_v_mva = angulo(vel_cassini, n_mva) * 180 / np.pi
print("el angulo theta_v del mva es ", theta_v_mva)

x14, x23 = ancho_mpb(t1, t2, t3, t4, n_mva, np.linalg.norm(vel_cassini))
print("el ancho de la mpb del MVA son (corto y largo)", x23, x14)

"""
Longitud inercial
"""
density_mean = 1.2
ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km
print(f"La longitud inercial de iones es {ion_length:1.3g} km")

Bup = np.mean(B[donde(t, t1 - 0.06): donde(t, t1), :], axis=0)
Bdown = np.mean(B[donde(t, t4): donde(t, t4 + 0.06), :], axis=0)

js, jv = corrientes(n_mva, Bup, Bdown, x23)

fuerza = np.cross(jv * 1e-9, B[donde(t, t4 + 0.06), :] * 1e-9)  # N/m^3 #en t4

print(
    f"la corriente superficial es {js} mA/m, su módulo {np.linalg.norm(js):1.3g} mA/m"
)
print(
    f"la corriente en volumen es {jv} nA/m^2 , su módulo {np.linalg.norm(jv):1.3g} mA/m^2"
)
print(f"la fuerza de Hall es {fuerza} N, su módulo {np.linalg.norm(fuerza):1.3g} N")

"""
presiones
"""
mu0 = 4 * np.pi * 1e-7  # Tm/A

p_mag = (
        (np.mean(np.linalg.norm(B[donde(t, t4): donde(t, tmpr), :], axis=1) * 1e-9)) ** 2
        / 2
        / mu0
)  # pmag en la MPR

print(f"La presión magnética en la MPR es {p_mag * 1e9:1.3g} nPa")
