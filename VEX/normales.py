import numpy as np

from _importar_datos import (
    importar_MAG,
    importar_t1t2t3t4,
    importar_tMVA,
)

from _fit import (
    calcular_normal,
    hallar_phi,
    rotacion,
)

import sys

sys.path.append("..")
from funciones import (
    donde,
    fechas,
    Mij,
    autovectores,
    SZA,
)
from funciones_metodos import bootstrap

np.set_printoptions(precision=4)

"""
calcula las normales de MVA y del fit, nada más.
"""


def MVA(B):
    M_ij = Mij(B)

    avec, lamb = autovectores(M_ij)

    print("la normal del MVA es ", avec[2])
    return avec[2]


year, month, day, doy = fechas()
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day)
ti, tf = importar_tMVA(year, month, day)
q_e = 1.6e-19  # carga electron #C

t, B, posicion, cl, tpos = importar_MAG(year, doy, ti - 1, tf + 1)
Bnorm = np.linalg.norm(B, axis=1)

inicio_MVA = donde(t, ti)
fin_MVA = donde(t, tf)
B_cut = B[inicio_MVA:fin_MVA]

n_mva = MVA(B_cut)

normal_boot, ang, delta_B3, out, out_phi = bootstrap(1000, B_cut)

print("la normal del bootstrap es ", normal_boot)

"""
El fit
"""

i_MVA = donde(tpos, ti)
f_MVA = donde(tpos, tf)
pos_MPB = int(0.5 * (f_MVA + i_MVA))

R = posicion[pos_MPB, :]
R_2d = np.array([R[0], np.sqrt(R[1] ** 2 + R[2] ** 2)])

sza_rad = SZA(posicion, pos_MPB) / 180 * np.pi
normal_2d = calcular_normal(sza_rad)

"""
A partir de la normal 2D, la puedo rotar y encontrar la normal 3D
Para eso, necesito hallar el ángulo phi primero
"""

phi = hallar_phi(R)[2]
normal_3d = rotacion(phi, normal_2d)

print("la normal del fit es ", normal_3d)

angulo_mva = np.arccos(np.clip(np.dot(normal_boot, normal_3d), -1.0, 1.0))

print(
    f"El ángulo entre las normales 3D de MVA y del fit es {angulo_mva * 180 / np.pi:.3g}º"
)
