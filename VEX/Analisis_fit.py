import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
from _importar_datos import (
    importar_MAG,
    importar_ELS_clweb,
    importar_t1t2t3t4,
    importar_tMVA,
)
from _fit import (
    plot_3d,
    plot_2D,
    calcular_normal,
    hallar_phi,
    rotacion,
    fit_3d,
    c_parametro,
)
import math

import sys

sys.path.append("..")
from funciones import (
    donde,
    fechas,
    tiempos,
    error,
    corrientes,
    ancho_mpb,
    Bpara_Bperp,
    find_nearest,
    SZA,
)

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

year, month, day, doy = fechas()
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day)
ti, tf = importar_tMVA(year, month, day)

t, B, posicion, cl, tpos = importar_MAG(year, doy, t1 - 1, t4 + 1)
Bnorm = np.linalg.norm(B, axis=1)
B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, t[0] + 0.2, t[-1] - 0.2)

i_mva = donde(t, ti)
f_mva = donde(t, tf)
sza = SZA(posicion, i_mva)

pos_MPB = int(0.5 * (f_mva + i_mva))

R = posicion[pos_MPB, :]
R_2d = np.array([R[0], np.sqrt(R[1] ** 2 + R[2] ** 2)])

sza_rad = SZA(posicion, pos_MPB) / 180 * np.pi
normal_2d = calcular_normal(sza_rad)
c = c_parametro(posicion, pos_MPB)

angulo_mva = np.arccos(np.clip(np.dot(n_mva2d, normal_2d), -1.0, 1.0))

print(
    f"El ángulo entre las normales 2D de MVA y del fit es {angulo_mva * 180 / np.pi:.3g}º"
)

plot_2D(pos, R_2d / 6050, normal_2d, c)
plt.title(f"{year}-{month}-{day}")
plt.show()
"""
A partir de la normal 2D, la puedo rotar y encontrar la normal 3D
Para eso, necesito hallar el ángulo phi primero
"""

phi = hallar_phi(R)[2]
normal_3d = rotacion(phi, normal_2d)
x, y, z = fit_3d(c)

plot_3d(x, y, z, R / 6050, normal_3d)

angulo_3d, theta_v, theta_Bn, x14, x23, Js, Jv = analisis(normal_3d)

print(
    f"El ángulo entre las normales 2D de MVA y del fit es {angulo_3d * 180 / np.pi:.3g}º"
)
print(f"Ancho MPB hmax = {x14:.3g}, hmin = {x23:.3g}")
print(
    f"Js = {Js} mA/m, |Js| = {np.linalg.norm(Js):.3g} mA/m \nJv = {Jv} nA/m², |Jv| = {np.linalg.norm(Jv):.3g} nA/m²"
)

v_punto = np.zeros((len(B) - 1, 3))
norma_v = np.zeros(len(B) - 1)
if cl:
    for i in range(len(v_punto)):
        v_punto[i, :] = posicion[i + 1, :] - posicion[i] / 10
        # en km/s, tiene resolución de 128Hz
        norma_v[i] = np.linalg.norm(v_punto[i, :])
else:
    for i in range(len(v_punto)):
        v_punto[i, :] = (posicion[i + 1, :] - posicion[i]) / (1 / 32)
        # en km/s, tiene resolución de 128Hz
        norma_v[i] = np.linalg.norm(v_punto[i, :])
# la velocidad promedio
v_media = np.mean(v_punto, axis=0)

x14, x23 = ancho_mpb(t1, t2, t3, t4, norm, v_media)
print(f"Ancho MPB hmax = {x14:.3g} km, hmin = {x23:.3g} km")

inicio_up = donde(t, t1 - 0.015)
fin_up = donde(t, t1)
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(t, t4)
fin_down = donde(t, t4 + 0.015)
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

J_s_MVA, J_v_MVA = corrientes(norm, B_upstream, B_downstream, x23)

fuerza_mva = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3 #en t4
print(
    f"Js = {J_s_MVA} mA/m, |Js| = {np.linalg.norm(J_s_MVA):.3g} mA/m \nJv = {J_v_MVA} nA/m², |Jv| = {np.linalg.norm(J_v_MVA):.3g} nA/m²"
)

plt.show()
E_Hall = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m
