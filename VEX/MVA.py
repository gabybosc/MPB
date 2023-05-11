import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import MultiCursor
from cycler import cycler
from importar_datos import importar_MAG_pds, importar_fila
from fit_venus import plot_orbita, fit_Xu
from multiplot import multi_plot_MAG_only
from update_parametros import (
    hoja_MVA_update,
    hoja_MVA_analisis,
    hoja_t1t2t3t4,
)
import sys

sys.path.append("..")
from funciones import (
    donde,
    fechas,
    Mij,
    error,
    corrientes,
    ancho_mpb,
    Bpara_Bperp,
    SZA,
)
from funciones_plot import hodograma

np.set_printoptions(precision=4)

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)


def MVA(B, year, month, day):
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

    av = np.concatenate([x1, x2, x3])

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

    # el error
    phi, delta_B3 = error(lamb, B, x)
    if phi[2, 1] > phi[2, 0]:
        error_normal = phi[2, 1] * 180 / np.pi
    else:
        error_normal = phi[2, 0] * 180 / np.pi

    print("MVA terminado")

    nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(year, month, day)

    # #######update la hoja de MVA
    hoja_MVA_update(hoja_mva, nr, lamb, av, error_normal, B3, delta_B3, B_norm_medio)
    return x3


year, month, day, doy = fechas()
ti, tf = 0, 24  # tiempos()
t, B, pos = importar_MAG_pds(year, doy, ti, tf)
Bnorm = np.linalg.norm(B, axis=1)

Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)

val = multi_plot_MAG_only(t, tpara, B, Bnorm, Bpara, Bperp, 2)
nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(year, month, day)

inicio_MVA = donde(t, min(val))
fin_MVA = donde(t, max(val))

sza = SZA(pos, inicio_MVA)
x3 = MVA(B[inicio_MVA:fin_MVA], year, month, day)

plt.show()
xx, yz = fit_Xu()
pos_RV = pos / 6050
orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)

plt.figure()
plot_orbita(pos_RV, orbita, xx, yz)

plt.quiver(
    pos_RV[inicio_MVA, 0],
    np.sqrt(pos_RV[inicio_MVA, 1] ** 2 + pos_RV[inicio_MVA, 2] ** 2),
    x3[0],
    np.sqrt(x3[1] ** 2 + x3[2] ** 2),
)


dd = np.loadtxt("../outputs/VEX_times.txt", usecols=1)
times = np.loadtxt("../outputs/VEX_times.txt")
idx = donde(dd, int(doy))
t1 = times[idx][2]
t2 = times[idx][3]
t3 = times[idx][4]
t4 = times[idx][5]
hoja_t1t2t3t4(hoja_parametros, nr, t1, t2, t3, t4)

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
hoja_MVA_analisis(
    hoja_mva, nr, t[inicio_MVA], t[fin_MVA], x14, x23, J_s_MVA, J_v_MVA, fuerza_mva
)
# E_Hall = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9) / (q_e * n_e)  # V/m

# buenas órbitas: SZA no tan alto, el campo en SW no es Bx
# 21 nov 2007
# 14 abr 2007
