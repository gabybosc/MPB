import numpy as np
import sys

sys.path.append("..")
from funciones import (
    Mij,
)
from funciones_plot import hodograma


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

    print(f"la normal del MVA es {x3}")

    # las proyecciones
    B1 = np.dot(B, x1)
    B2 = np.dot(B, x2)
    B3 = np.dot(B, x3)

    # el B medio
    B_medio_vectorial = np.mean(B, axis=0)

    print(f"cociente de lambdas = {lamb[1] / lamb[2]:.4g}")
    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    print(f"El B medio es {B_medio_vectorial}, con norma {B_norm_medio:.3g}")
    hodograma(B1, B2, B3)

    return B_medio_vectorial, x3
