import numpy as np
import sys

sys.path.append("..")

from funciones import (
    Mij,
    fechas,
    donde,
)
from importar_datos import importar_mag, importar_t1t2t3t4


"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre.
En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro.
El mapa de colores me va a dar el valor del cociente.
"""


def MVA(t, ti, tf, B):
    inicio = donde(t, ti)
    fin = donde(t, tf)

    B_cut = B[inicio : fin + 1, :]

    M_ij = Mij(B_cut)

    [lamb, x] = np.linalg.eigh(M_ij)
    idx = lamb.argsort()[::-1]
    lamb = lamb[idx]
    x = x[:, idx]
    # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
    x1 = x[:, 0]
    x2 = x[:, 1]
    x3 = x[:, 2]
    if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
        x3 = -x3

    cociente = lamb[1] / lamb[2]

    B1 = np.dot(B_cut, x1)
    B2 = np.dot(B_cut, x2)
    B3 = np.dot(B_cut, x3)

    return (cociente, B1, B2, B3, B_cut)


def escalas_lambda(year, month, day, doy, hora, tt):
    ti = tt[0]
    tf = tt[3]

    mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

    tiempo_central = np.zeros(
        int((tf - ti) * 3600)
    )  # la cantidad de segundos entre tf y ti
    tiempo_central[0] = ti
    for i in range(len(tiempo_central) - 1):
        tiempo_central[i + 1] = (
            tiempo_central[i] + 1 / 3600
        )  # el tiempo central se va barriendo cada 5 segundos

    escalas = np.zeros(60)
    escalas[0] = 1 / 3600  # la escala más chica es de 1s
    for i in range(len(escalas) - 1):
        escalas[i + 1] = escalas[i] + 1 / 3600

    cociente = np.zeros((len(tiempo_central), len(escalas)))
    for i in range(len(tiempo_central)):
        for j in range(len(escalas)):
            ratio, B1, B2, B3, B_cut = MVA(
                t, tiempo_central[i] - escalas[j], tiempo_central[i] + escalas[j], B
            )
            cociente[i, j] = ratio

    m = 1.67e-27  # kg
    q = 1.6e-19  # C
    periodo_ciclotron = (
        2 * np.pi * m / (q * np.linalg.norm(B_cut, axis=1)) * 1e9
    )  # en s
    periodo_diezmado = np.zeros(len(tiempo_central))
    k = len(periodo_ciclotron) / len(tiempo_central)
    for i in range(len(periodo_diezmado)):
        periodo_diezmado[i] = periodo_ciclotron[int(i * k)]

    matriz = np.zeros((len(tiempo_central) + 1, len(escalas) + 2))
    matriz[0, 2:] = escalas
    matriz[1:, 0] = periodo_diezmado
    matriz[1:, 1] = tiempo_central
    matriz[1:, 2:] = cociente

    with open(f"../outputs/cociente_lambdas_d{doy}_t{hora}.txt", "w") as file:
        # with open(f'outputs/cociente_lambdas_salida_d{doy}_t{hora}.txt','w') as file:
        file.write(
            "La primera columna es el período de ciclotron, la primera fila son las escalas, la segunda columna es el tiempo central.\n"
        )
        for i in range(len(matriz[:, 0])):
            for j in range(len(matriz[0, :])):
                file.write(f"{matriz[i,j]}\t")
            file.write("\n")
    # print(f'{l / len(dates) * 100}%')
    return (B, t, escalas, cociente, tiempo_central)
