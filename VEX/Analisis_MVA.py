import numpy as np
import sys

from _importar_datos import importar_MAG, importar_t1t2t3t4, importar_tMVA
from _old_fit_venus import plot_orbita, fit_Xu, fit_R
import math

sys.path.append("..")
from funciones_plot import onpick1, hodograma
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

"""
tdec, Bx, By, Bz, modulo B,pos x, pos y, pos z, distancia km
"""
year, month, day, doy = fechas()
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day)
ti, tf = importar_tMVA(year, month, day)

t, B, posicion, cl, tpos = importar_MAG(year, doy, t1 - 1, t4 + 1)
Bnorm = np.linalg.norm(B, axis=1)
B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, t[0] + 0.2, t[-1] - 0.2)

i_mva = donde(t, ti)
f_mva = donde(t, tf)


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

    SZA = angulo(posicion[n_p, :], [1, 0, 0]) * 180 / np.pi
    print(f"altitud = {altitud_media}, SZA = {SZA}")

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
    print("MVA terminado")
    return avec[2], B, t, posicion


def corte(t, ti, tf, vector):
    inicio = donde(t, ti)
    fin = donde(t, tf)
    if vector.ndim > 1:
        vector_cut = np.nanmean(vector[inicio:fin, :], axis=0)
    else:
        vector_cut = np.nanmean(vector[inicio:fin])
    return vector_cut


x3, B_cut, t_cut, posicion_cut = MVA(
    t[i_mva:f_mva], B[i_mva:f_mva], posicion[i_mva:f_mva]
)

B_upstream = corte(t, t1 - 0.015, t1, B)
B_downstream = corte(t, t4, t4 + 0.015, B)

omega = angulo(B_upstream, B_downstream)
print(f"omega = {omega * 180 / np.pi}")

v_punto = np.zeros((len(posicion_cut) - 1, 3))
norma_v = np.zeros(len(posicion_cut) - 1)
for i in range(len(v_punto)):
    v_punto[i, :] = posicion_cut[i + 1, :] - posicion_cut[i]
    # en km/s, tiene resolución de 1Hz
    norma_v[i] = np.linalg.norm(v_punto[i, :])

# la velocidad promedio
v_media = np.mean(v_punto, axis=0)
print("v media", v_media)
print(f"El ángulo entre v media y la normal es {angulo(v_media, x3) * 180 / np.pi}")

x14, x23 = ancho_mpb(t1, t2, t3, t4, x3, v_media)
print(f"El ancho x14 es {x14}km y el ancho x23={x23} km")

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, x23)

print(
    f"J sup = {J_s_MVA}, |Js| = {np.linalg.norm(J_s_MVA)}, J vol {J_v_MVA}, |Jv| = {np.linalg.norm(J_v_MVA)}"
)
inicio_down = donde(t, t1 - 0.015)
# fuerza_mva = fuerza(J_v_MVA, B[inicio_down, :])

# altitude(posicion, 2574)
print("SZA para t1, t2, t3, t4:\n")
for ts in [t1, t2, t3, t4]:
    print(ts, SZA(posicion, donde(t, ts)))
