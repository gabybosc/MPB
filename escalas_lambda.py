"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre.
En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro.
El mapa de colores me va a dar el valor del cociente.
"""

import numpy as np
import time as time
from funciones import find_nearest, Mij, fechas
from importar_datos import importar_t1t2t3t4, importar_mag

# dates = np.loadtxt('outputs/t1t2t3t4.txt')

# for l in range(len(dates)):


year, month, day, doy = fechas()
hora = input("Hora en HH\n")

tiempos = importar_t1t2t3t4(year, month, day, int(hora))

# ti = float(input('Tiempo inicial del barrido\n'))
# tf = float(input('Tiempo final del barrido\n'))
ti = tiempos[0]
tf = tiempos[3]

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

tiempo_central = np.zeros(
    int((tf - ti) * 3600)
)  # la cantidad de segundos entre tf y ti
tiempo_central[0] = ti
for i in range(len(tiempo_central) - 1):
    tiempo_central[i + 1] = (
            tiempo_central[i] + 1 / 3600
    )  # el tiempo central se va barriendo cada 5 segundos

# print(f'tiempo final efectivo = {tiempo_central[-1]}')

escalas = np.zeros(60)
escalas[0] = 1 / 3600  # la escala más chica es de 1s
for i in range(len(escalas) - 1):
    escalas[i + 1] = escalas[i] + 1 / 3600

# print(f'escala mayor = {escalas[-1]*3600}s')


program_starts = time.time()
cociente = np.zeros((len(tiempo_central), len(escalas)))
for i in range(len(tiempo_central)):
    for j in range(len(escalas)):
        inicio = np.where(t == find_nearest(t, tiempo_central[i] - escalas[j]))[0][0]
        fin = np.where(t == find_nearest(t, tiempo_central[i] + escalas[j]))[0][0]

        # ahora empieza el MVA con los datos que elegí
        B_cut = B[inicio: fin + 1, :]

        M_ij = Mij(B_cut)

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

        cociente[i, j] = lamb[1] / lamb[2]

program_ends = time.time()

print(f"El loop tardó {program_ends - program_starts:.2f} s")

m = 1.67e-27  # kg
q = 1.6e-19  # C
periodo_ciclotron = 2 * np.pi * m / (q * np.linalg.norm(B_cut, axis=1)) * 1e9  # en s
periodo_diezmado = np.zeros(len(tiempo_central))
k = len(periodo_ciclotron) / len(tiempo_central)
for i in range(len(periodo_diezmado)):
    periodo_diezmado[i] = periodo_ciclotron[int(i * k)]

matriz = np.zeros((len(tiempo_central) + 1, len(escalas) + 2))
matriz[0, 2:] = escalas
matriz[1:, 0] = periodo_diezmado
matriz[1:, 1] = tiempo_central
matriz[1:, 2:] = cociente

with open(f"outputs/cociente_lambdas_d{doy}_t{hora}.txt", "w") as file:
    # with open(f'outputs/cociente_lambdas_salida_d{doy}_t{hora}.txt','w') as file:
    file.write(
        "La primera columna es el período de ciclotron, la primera fila son las escalas, la segunda columna es el tiempo central.\n"
    )
    for i in range(len(matriz[:, 0])):
        for j in range(len(matriz[0, :])):
            file.write(f"{matriz[i, j]}\t")
        file.write("\n")
# print(f'{l / len(dates) * 100}%')
