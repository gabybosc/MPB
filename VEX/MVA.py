import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

# import datetime as dt
# import math

from _importar_datos import (
    importar_MAG,
    importar_fila,
    importar_t1t2t3t4,
    importar_tMVA,
    importar_BupBdown,
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
from _update_parametros import (
    hoja_param,
    hoja_MVA_update,
    hoja_MVA_analisis,
    hoja_t1t2t3t4,
    hoja_bootstrap_p1,
    hoja_bootstrap_p2,
    hoja_fit,
)

import sys

sys.path.append("..")
from funciones import (
    donde,
    fechas,
    tiempos,
    Bpara_Bperp,
    Mij,
    autovectores,
    error,
    find_nearest,
    angulo,
    ancho_mpb,
    corrientes,
    UTC_to_hdec,
    SZA,
    altitude,
)
from funciones_metodos import bootstrap, plot_bootstrap
from funciones_plot import hodograma

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

"""
hardcodeado para el 28 de octubre y que calcule todo lo que quiero
"""


def velocidad(posicion_cut, tpos):
    M = len(posicion_cut)
    v_punto = np.zeros((M - 1, 3))
    deltat = np.zeros(M - 1)
    if np.mean(posicion_cut) < 10:  # es decir, pos está en RV en vez de km
        posicion_cut = posicion_cut * 6050
    for i in range(len(v_punto)):
        deltat[i] = (tpos[i + 1] - tpos[i]) * 3600  # delta t en segundos
        v_punto[i] = (posicion_cut[i + 1, :] - posicion_cut[i, :]) / deltat[i]
        # en km/s
    # la velocidad promedio
    v_media = np.mean(v_punto, axis=0)
    return v_media


def analisis(normal, B_medio_vectorial, v_media):
    angulo_vs_mva = angulo(normal, x3)
    B_medio_normalizado = B_medio_vectorial / np.linalg.norm(
        B_medio_vectorial
    )  # la normalizo por cuestion de cuentas nomas
    angulo_v = angulo(normal, np.linalg.norm(v_media))
    angulo_B = angulo(normal, B_medio_normalizado)
    x_14, x_23 = ancho_mpb(t1, t2, t3, t4, normal, v_media)
    J_s, J_v = corrientes(normal, B_upstream, B_downstream, x_23)

    return angulo_vs_mva, angulo_v, angulo_B, x_14, x_23, J_s, J_v


def vmedia(pos):
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


def MVA(t, B, posicion):
    M = len(t)

    n_p = int(M / 2)

    M_ij = Mij(B)

    avec, lamb = autovectores(M_ij)

    print("la normal del MVA es ", avec[2])

    # las proyecciones
    B1 = np.dot(B, avec[0])
    B2 = np.dot(B, avec[1])
    B3 = np.dot(B, avec[2])

    # el B medio
    B_medio_vectorial = np.mean(B, axis=0)
    altitud = np.linalg.norm(posicion, axis=1) - 3390  # km
    altitud_media = np.mean(altitud)

    # SZA = angulo(posicion[n_p, :], [1, 0, 0]) * 180 / np.pi
    # print(f"altitud = {altitud_media}, SZA = {SZA}")

    print(f"l1 = {lamb[0]}, l2 = {lamb[1]}, l3 = {lamb[2]}")
    print("cociente de lambdas = ", lamb[1] / lamb[2])
    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    print(f"El B medio es {B_medio_vectorial}, su norma es {B_norm_medio}\n")
    print(f"La componente normal media del B es {np.mean(B3)}\n")
    print(
        r"El cociente <B_3>/|<B>|$ es",
        f"{np.mean(B3) / B_norm_medio}",
    )
    hodograma(B1, B2, B3)

    print(
        f"El ángulo entre el vector de campo magnético medio y la normal es {angulo(B_medio_vectorial, avec[2]) * 180 / np.pi} "
    )

    # el error
    phi, delta_B3 = error(lamb, B, avec[2])
    phi = 0
    delta_B3 = 0
    hoja_MVA_update(hoja_mva, nr, lamb, avec[2], phi, B3, delta_B3, B_norm_medio)
    print("MVA terminado")
    return avec[2], B, B_medio_vectorial


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


year, month, day, doy = fechas()
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day)
ti, tf = importar_tMVA(year, month, day)
q_e = 1.6e-19  # carga electron #C
nr, hoja_parametros, hoja_mva, hoja_boot, hoja_ajuste = importar_fila(year, month, day)

t, B, posicion, cl, tpos = importar_MAG(year, doy, ti - 1, tf + 1)
Bnorm = np.linalg.norm(B, axis=1)

inicio_MVA = donde(t, ti)
fin_MVA = donde(t, tf)

# sza = SZA(posicion, inicio_MVA)
x3, B_cut, B_medio_vectorial = MVA(
    t[inicio_MVA:fin_MVA], B[inicio_MVA:fin_MVA], posicion[inicio_MVA:fin_MVA]
)
v_media = velocidad(
    posicion[donde(tpos, ti): donde(tpos, tf)], tpos[donde(tpos, ti): donde(tpos, tf)]
)

B_upstream, B_downstream = importar_BupBdown(year, month, day)

omega = np.arccos(
    np.dot(B_upstream, B_downstream)
    / (np.linalg.norm(B_upstream) * np.linalg.norm(B_downstream))
)

"""
MVA
"""

angulo_cero, angulo_v_mva, angulo_B_mva, x14, x23, J_s_MVA, J_v_MVA = analisis(
    x3, B_medio_vectorial, v_media
)
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
hoja_MVA_analisis(hoja_mva, nr, ti, tf, x14, x23, J_s_MVA, J_v_MVA)

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
) = analisis(normal_boot, B_medio_vectorial)
print("\nRESULTADOS BOOTSTRAP\n")
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

# nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(year, month, day)
# hoja_bootstrap_p2(
#     hoja_boot,
#     nr,
#     angulo_v_boot,
#     angulo_B_boot,
#     x23_boot,
#     x14_boot,
#     J_s_boot,
#     J_v_boot,
#     0,
#     [0, 0, 0],
# )


"""
El fit
"""
print("\nRESULTADOS FIT\n")

# year, doy = 2008, 302
# date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)
# month = date_orbit.strftime("%m")
# day = date_orbit.strftime("%d")
#
# ti_MVA, tf_MVA = 8.5444425, 8.549442222
# t1, t2, t3, t4 = [8.541556528, 8.544405851, 8.551476393, 8.556111111]
# t, B, posicion, cl, tpos = importar_MAG(year, doy, t1 - 0.5, t4 + 0.5)
# Bnorm = np.linalg.norm(B, axis=1)

i_MVA = donde(tpos, ti)
f_MVA = donde(tpos, tf)
pos_MPB = int(0.5 * (f_MVA + i_MVA))

R = posicion[pos_MPB, :]
R_2d = np.array([R[0], np.sqrt(R[1] ** 2 + R[2] ** 2)])

sza_rad = SZA(posicion, pos_MPB) / 180 * np.pi
normal_2d = calcular_normal(sza_rad)
c = c_parametro(posicion, pos_MPB)

n_mva2d = np.array([x3[0], np.sqrt(x3[1] ** 2 + x3[2] ** 2)])
n_mva = x3

angulo_mva = np.arccos(np.clip(np.dot(n_mva2d, normal_2d), -1.0, 1.0))

print(
    f"El ángulo entre las normales 2D de MVA y del fit es {angulo_mva * 180 / np.pi:.3g}º"
)

plot_2D(posicion, R_2d / 6050, normal_2d, c)
plt.title(f"{year}-{month}-{day}")
plt.show()

"""
A partir de la normal 2D, la puedo rotar y encontrar la normal 3D
Para eso, necesito hallar el ángulo phi primero
"""

phi = hallar_phi(R)[2]
normal_3d = rotacion(phi, normal_2d)
x, y, z = fit_3d(c)

plot_3d(x, y, z, R / 6050, normal_3d, x3)

angulo_3d, theta_v, theta_Bn, x14, x23, Js, Jv = analisis(normal_3d, B_medio_vectorial)

print(
    f"El ángulo entre las normales 2D de MVA y del fit es {angulo_3d * 180 / np.pi:.3g}º"
)
print(f"Ancho MPB hmax = {x14:.3g}, hmin = {x23:.3g}")
print(
    f"Js = {Js} mA/m, |Js| = {np.linalg.norm(Js):.3g} mA/m \nJv = {Jv} nA/m², |Jv| = {np.linalg.norm(Jv):.3g} nA/m²"
)

hoja_fit(
    hoja_ajuste,
    nr,
    np.array([0.11, -0.22, c]),
    normal_3d,
    angulo_3d,
    theta_Bn,
    theta_v,
    x14,
    x23,
    Js,
    Jv,
)
