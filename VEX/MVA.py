import numpy as np
import matplotlib.pyplot as plt
import sys
from importar_datos import importar_VEX_mag_AMDA

sys.path.append("..")
from funciones import donde, fechas, tiempos, Mij, error, corrientes
from funciones_plot import hodograma


np.set_printoptions(precision=4)


def MVA(year, month, day, ti_MVA, tf_MVA):

    t, B = importar_VEX_mag_AMDA(year, month, day, ti_MVA, tf_MVA)
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

    # el error
    phi, delta_B3 = error(lamb, B, x)
    print("MVA terminado")
    return x3, B, t


# year, month, day, doy = fechas()
# ti, tf = tiempos()s
year, month, day = 2007, 11, 21
# ti, tf = 2.26012, 2.26799
ti = 2.26095222
tf = 2.2667852777777777
x3, B, t = MVA(year, month, day, ti, tf)
plt.show()

t1 = ti
t4 = tf
inicio_up = donde(t, t1 - 0.015)
fin_up = donde(t, t1)
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(t, t4)
fin_down = donde(t, t4 + 0.015)
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, x_23_MVA)

fuerza_mva = np.cross(J_v_MVA * 1e-9, B[inicio_down, :] * 1e-9)  # N/m^3 #en t4
# buenas órbitas: SZA no tan alto, el campo en SW no es Bx
# 21 nov 2007
# 14 abr 2007
