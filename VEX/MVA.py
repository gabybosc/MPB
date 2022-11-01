import numpy as np
import matplotlib.pyplot as plt
import sys
from glob import glob

sys.path.append("..")
from funciones import (
    donde,
    fechas,
    tiempos,
    Mij,
    error,
)
from funciones_plot import hodograma


np.set_printoptions(precision=4)


def importar_VEX_mag_AMDA(year, month, day, ti, tf):
    path = glob(f"../../../datos/VEX MAG/{day}{month}{year}/*.txt")
    B = np.genfromtxt(path[0], usecols=[1, 2, 3])
    tt = np.genfromtxt(path[0], usecols=0, dtype="str")
    fecha = np.array([x.split("T") for x in tt])
    hora = np.array([x.split(":") for x in fecha[:, 1]])
    hh = np.array([int(x) for x in hora[:, 0]])
    mm = np.array([int(x) for x in hora[:, 1]])
    ss = np.array([float(x) for x in hora[:, 2]])
    t = hh + mm / 60 + ss / 3600  # hdec
    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    return t_cut, B_cut


year, month, day, doy = fechas()
ti, tf = tiempos()


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


x3, B, t = MVA(year, month, day, ti, tf)
plt.show()

# buenas órbitas: SZA no tan alto, el campo en SW no es Bx
# 21 nov 2007
# 14 abr 2007
