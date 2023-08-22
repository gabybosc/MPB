import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
from leer_datos import importar_bepi


import sys

sys.path.append("..")
from funciones import (
    donde,
    Mij,
    corrientes,
    ancho_mpb,
    Bpara_Bperp,
    SZA,
    UTC_to_hdec,
)
from funciones_plot import hodograma

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)


def MVA(B):
    # ya los importa cortados a los datos, entonces no hace falta que haga el cut yo
    M_ij = Mij(B)

    # ahora quiero los autovectores y autovalores
    [lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

    # Los ordeno de mayor a menor
    idx = lamb.argsort()[::-1]
    lamb = lamb[idx]
    x = x[:, idx]
    # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
    x1 = x[:, 0]
    x2 = x[:, 1]
    x3 = x[:, 2]

    if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
        x3 = -x3
    if any(np.cross(x1, x2) - x3) > 0.01:
        print("Cambio el signo de x1 para que los av formen terna derecha")
        x1 = -x1

    print("la normal del MVA es ", x3)

    # las proyecciones
    B1 = np.dot(B, x1)
    B2 = np.dot(B, x2)
    B3 = np.dot(B, x3)

    # el B medio
    B_medio_vectorial = np.mean(B, axis=0)

    print("cociente de lambdas = ", lamb[1] / lamb[2])
    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    print(f"El B medio es {B_norm_medio}")
    hodograma(B1, B2, B3)

    return x3


ti = 13.5
tf = 14.1

t, B, pos = importar_bepi(ti, tf)
Bnorm = np.linalg.norm(B, axis=1)

Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)

ti_mva, tf_mva = 13.8933, 13.9094
inicio_MVA = donde(t, ti_mva)
fin_MVA = donde(t, tf_mva)

sza = SZA(pos, inicio_MVA)
x3 = MVA(B[inicio_MVA:fin_MVA])

t1, t2, t3, t4 = 13.8807, 13.8933, 13.9094, 13.929

v_punto = np.zeros((len(B) - 1, 3))
norma_v = np.zeros(len(B) - 1)
for i in range(len(v_punto)):
    v_punto[i, :] = (pos[i + 1, :] - pos[i]) / (1 / 32)
    # en km/s, tiene resolución de 32Hz
    norma_v[i] = np.linalg.norm(v_punto[i, :])
# la velocidad promedio
v_media = np.mean(v_punto, axis=0)


x14, x23 = ancho_mpb(t1, t2, t3, t4, x3, v_media)
print(f"Ancho MPB hmax = {x14:.3g}, hmin = {x23:.3g}")

inicio_up = donde(t, t1 - 0.015)
fin_up = donde(t, t1)
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(t, t4)
fin_down = donde(t, t4 + 0.015)
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, x23)

fuerza_mva = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3 #en t4
print(
    f"Js = {J_s_MVA} mA/m, |Js| = {np.linalg.norm(J_s_MVA):.3g} mA/m \nJv = {J_v_MVA} nA/m², |Jv| = {np.linalg.norm(J_v_MVA):.3g} nA/m²"
)

plt.show()
# E_Hall = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m

# buenas órbitas: SZA no tan alto, el campo en SW no es Bx
# 21 nov 2007
# 14 abr 2007
