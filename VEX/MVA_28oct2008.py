import numpy as np
import matplotlib as mpl
from cycler import cycler
from _importar_datos import importar_MAG, importar_fila

# import matplotlib.pyplot as plt
# from _fit_venus import plot_orbita, fit_Xu, fit_R
# import math
from _update_parametros import (
    hoja_param,
    hoja_MVA_update,
    hoja_MVA_analisis,
    hoja_t1t2t3t4,
    hoja_bootstrap_p1,
    hoja_bootstrap_p2,
)
from MVA import MVA

import sys

sys.path.append("..")
from funciones import (
    donde,
    corrientes,
    ancho_mpb,
    SZA,
    angulo,
)
from funciones_metodos import bootstrap, plot_bootstrap

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

"""
hardcodeado para el 28 de octubre y que calcule todo lo que quiero
"""


def analisis(normal):
    angulo_vs_mva = angulo(normal, x3)
    B_medio_normalizado = B_medio_vectorial / np.linalg.norm(
        B_medio_vectorial
    )  # la normalizo por cuestion de cuentas nomas
    v_media_norm = v_media / np.linalg.norm(v_media)
    angulo_v = angulo(normal, v_media_norm)
    angulo_B = angulo(normal, B_medio_normalizado)
    x_14, x_23 = ancho_mpb(t1, t2, t3, t4, normal, v_media)
    J_s, J_v = corrientes(normal, B_upstream, B_downstream, x_23)

    return angulo_vs_mva, angulo_v, angulo_B, x_14, x_23, J_s, J_v


def vmedia():
    v_punto = np.zeros((len(B) - 1, 3))
    norma_v = np.zeros(len(B) - 1)
    if cl:
        for i in range(len(v_punto)):
            v_punto[i, :] = pos[i + 1, :] - pos[i] / 10
            # en km/s, tiene resolución de 128Hz
            norma_v[i] = np.linalg.norm(v_punto[i, :])
    else:
        for i in range(len(v_punto)):
            v_punto[i, :] = (pos[i + 1, :] - pos[i]) / (1 / 32)
            # en km/s, tiene resolución de 128Hz
            norma_v[i] = np.linalg.norm(v_punto[i, :])
    # la velocidad promedio
    return np.mean(v_punto, axis=0)


def bootstrap_completo(B, N=1000):
    normal, phi, delta_B3, out, out_phi = bootstrap(1000, B[inicio_MVA:fin_MVA])

    muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

    B3 = np.dot(B, normal)

    B_medio_vectorial = np.mean(B, axis=0)
    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    if sigma31 > sigma32:
        error = sigma31
    else:
        error = sigma32

    print("Fin del bootstrap. ")
    nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(year, month, day)
    M = fin_MVA - inicio_MVA
    hoja_bootstrap_p1(hoja_boot, nr, M, N, normal, error, B3, B_norm_medio, sigmaB)

    return normal


year, month, day, doy = (2008, 10, 28, 302)
t1, t2, t3, t4 = (
    8.541556527954604,
    8.544405851015947,
    8.551476393427427,
    8.556111111111111,
)
ti, tf = (t2, 8.549442222)
q_e = 1.6e-19  # carga electron #C

t, B, pos, cl, tpos = importar_MAG(year, doy, ti - 1, tf + 1)
Bnorm = np.linalg.norm(B, axis=1)

inicio_MVA = donde(t, ti)
fin_MVA = donde(t, tf)

sza = SZA(pos, inicio_MVA)
B_medio_vectorial, x3 = MVA(B[inicio_MVA:fin_MVA])

v_media = vmedia()

inicio_up = donde(t, t1 - 0.015)
fin_up = donde(t, t1)
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(t, t4)
fin_down = donde(t, t4 + 0.015)
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

omega = np.arccos(
    np.dot(B_upstream, B_downstream)
    / (np.linalg.norm(B_upstream) * np.linalg.norm(B_downstream))
)

"""
MVA
"""

angulo_cero, angulo_v_mva, angulo_B_mva, x14, x23, J_s_MVA, J_v_MVA = analisis(x3)
print("RESULTADOS MVA\n")
print(f"Ancho MPB hmax = {x14:.3g} km, hmin = {x23:.3g} km")
print(
    f"Js = {J_s_MVA} mA/m, |Js| = {np.linalg.norm(J_s_MVA):.4g} mA/m \nJv = {J_v_MVA} nA/m², |Jv| = {np.linalg.norm(J_v_MVA):.4g} nA/m²"
)

print(
    f"theta_Bn = {angulo_B_mva * 180 / np.pi:.3g}, theta_vn = {angulo_v_mva * 180 / np.pi:.3g}"
)
# E_Hall = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m

# nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(year, month, day)
# hoja_param(hoja_parametros, nr, sza, v_media, omega, B_upstream, B_downstream)

"""
Bootstrap
"""
normal_boot = bootstrap_completo(B)
(
    angulo_vs_mva,
    angulo_v_boot,
    angulo_B_boot,
    x14_boot,
    x23_boot,
    J_s_boot,
    J_v_boot,
) = analisis(normal_boot)
print("RESULTADOS BOOTSTRAP\n")
print(f"normal bootstrap: {normal_boot}")
print(
    f"ángulo que forma la normal del boostrap respecto a la del MVA: {angulo_vs_mva * 180 / np.pi}º"
)
print(f"Ancho MPB hmax = {x14_boot:.3g} km, hmin = {x23_boot:.3g} km")
print(
    f"Js = {J_s_boot} mA/m, |Js| = {np.linalg.norm(J_s_boot):.4g} mA/m \nJv = {J_v_boot} nA/m², |Jv| = {np.linalg.norm(J_v_boot):.4g} nA/m²"
)

print(
    f"theta_Bn = {angulo_B_boot * 180 / np.pi:.3g}, theta_vn = {angulo_v_boot * 180 / np.pi:.3g}"
)

nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(year, month, day)
hoja_bootstrap_p2(
    hoja_boot,
    nr,
    angulo_v_boot,
    angulo_B_boot,
    x23_boot,
    x14_boot,
    J_s_boot,
    J_v_boot,
    0,
    [0, 0, 0],
)
