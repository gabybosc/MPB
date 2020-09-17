import numpy as np
import time as time
import sys
import gif
import matplotlib.pyplot as plt
from importar_datos import importar_mag

sys.path.append("..")

from funciones import find_nearest, Mij, fechas, angulo, hdec_to_UTC

"""
Hace el MVA en un barrido en ventanas (nesting) para ver cómo varía la normal.
Entonces, para cada tiempo t entre t1 y t4 hace el MVA y calcula la normal en
radios de 15s a 45s. Luego, calcula el ángulo entre las diferentes normales para
cada t, lo que me devuelve una matriz con ceros en la antidiagonal y que es
persimétrica (simétrica en la antidiagonal). Hace el heatmap de estos ángulos
para cada t y finalmente me devuelve todo en un gif con un cuadro cada 2 segundos.
"""

year, month, day, doy = 2016, "03", 16, "076"  # fechas()
hora = 18  # input("Hora en HH\n")

tiempos_txt = np.loadtxt("../outputs/t1t2t3t4.txt")
for i in range(len(tiempos_txt)):
    if (
        int(year) == int(tiempos_txt[i, 0])
        and int(doy) == int(tiempos_txt[i, 1])
        and int(hora) == int(tiempos_txt[i, 2])
    ):
        tiempos = [
            tiempos_txt[i, 2],
            tiempos_txt[i, 3],
            tiempos_txt[i, 4],
            tiempos_txt[i, 5],
        ]

ti = tiempos[1]
tf = tiempos[2]

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

tiempo_central = np.zeros(int((tf - ti) * 3600))
tiempo_central[0] = ti
for i in range(len(tiempo_central) - 1):
    tiempo_central[i + 1] = tiempo_central[i] + 1 / 3600


"""
Me va a dar las normales en un rango de 30s a 90s en diámetro, con 1s de diferencia.
Empieza en 30s porque nos nos interesa ver la normal adentro de la hoja de
corriente (es la que ya tenemos), pero sí ver cómo varía a medida que aumentamos
la cantidad de puntos
"""
escalas = np.zeros(30)
escalas[0] = 15 / 3600
for i in range(len(escalas) - 1):
    escalas[i + 1] = escalas[i] + 1 / 3600  # va de 15s a 45s

program_starts = time.time()
normales = np.zeros((len(tiempo_central), len(escalas), 3))
# un array en 3 dimensiones
for i in range(len(tiempo_central)):
    for j in range(len(escalas)):
        inicio = np.where(t == find_nearest(t, tiempo_central[i] - escalas[j]))[0][0]
        fin = np.where(t == find_nearest(t, tiempo_central[i] + escalas[j]))[0][0]

        B_cut = B[inicio : fin + 1, :]

        M_ij = Mij(B_cut)

        [lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

        idx = lamb.argsort()[::-1]
        lamb = lamb[idx]
        x = x[:, idx]
        # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
        x1 = x[:, 0]
        x2 = x[:, 1]
        x3 = x[:, 2]
        if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
            x3 = -x3

        normales[i, j] = x3

program_ends = time.time()

print(f"El loop tardó {program_ends-program_starts:.2f} s")

"""
Una vez que está normales lleno, elijo un tiempo y las normales en todo el "radio"
que le corresponde. Hago un gif con el ángulo entre las normales para cada radio
centrado en diferentes tiempos (cada cuadro es un tiempo central, cada heatmap
es una matriz de los ángulos). Es decir, calculo el ángulo a tiempo fijo para
cada intervalo.
Finalmente, un histograma con el número de veces que aparece cada ángulo de 0 a
20, sin contar los ceros de la diagonal.
"""

minimo = escalas[0] * 3600
maximo = escalas[-1] * 3600


@gif.frame
def plot(k, angle, minimo=minimo, maximo=maximo):
    plt.imshow(
        angle,
        aspect="auto",
        origin="lower",
        extent=(minimo, maximo, minimo, maximo),
        cmap="viridis",
        vmin=0,
        vmax=10,
    )
    cbar = plt.colorbar()
    cbar.set_label("Angle between normals (º)", rotation=270)
    plt.title(f"Heatmap centered in t = {k}")
    plt.xlabel("Radius (s)")
    plt.ylabel("Radius (s)")


# frames = []
n = 1
b = [i for i in range(20)]
histograma = 0
for k in range(int(len(tiempo_central) / n)):
    NN = int(k * n)
    # tiempo_central[NN]
    normales_t1 = normales[NN]
    angle = np.zeros((len(normales_t1), len(normales_t1)))
    thdec = hdec_to_UTC(tiempo_central[NN])
    for i in range(len(normales_t1)):
        for j in range(len(normales_t1)):
            angle[i, j] = angulo(normales_t1[i], normales_t1[j]) * 180 / np.pi
    # frame = plot(thdec, angle)
    # frames.append(frame)
    hist, bins = np.histogram(angle, bins=b)
    hist[0] = hist[0] - 29  # saco los 29 ceros de la diagonal
    histograma += hist

t = 0
s = 0
for i in range(len(histograma)):
    s += histograma[i] * ((bins[i] + bins[i + 1]) / 2)
mean = s / np.sum(histograma)
for i in range(len(histograma)):
    t += histograma[i] * (bins[i] - mean) ** 2
std = np.sqrt(t / np.sum(histograma))
# gif.save(frames, "heatmap.gif", duration=500)

center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, histograma)
plt.axvline(x=std, ls="--", color="orange", alpha=0.7, label=r"$\sigma$")
plt.axvline(x=2 * std, ls="--", color="red", alpha=0.7, label=r"2 $\sigma$")
plt.xlabel("Angle (º)")
plt.ylabel("Counts")
plt.title(
    "Angle between normals centered at the same time, different interval\n(between t1 and t4)"
)
plt.legend()
plt.show()

"""
Histograma a tiempo fijo variando el intervalo (lo mismo que antes).
Ángulo entre estas normales y la normal del MVA.
"""
for bines in [10, 20]:
    b = [i for i in range(bines)]
    normal = [0.920, -0.302, 0.251]
    histograma = 0
    for k in range(len(tiempo_central)):
        normales_t1 = normales[k, :, :]
        angle = np.zeros(len(normales_t1))
        for i in range(len(normales_t1)):
            angle[i] = angulo(normal, normales_t1[i]) * 180 / np.pi
        hist, bins = np.histogram(angle, bins=b)
        histograma += hist

    center = (bins[:-1] + bins[1:]) / 2
    plt.figure()
    plt.bar(center, histograma)
    plt.xlabel("Angle (º)")
    plt.ylabel("Counts")
    plt.title("Fixed time, varying interval.\nAngle with respect to MVA normal.")
plt.show()

"""
Histograma a intervalo fijo pero variando el tiempo. Ángulo entre estas normales
y la normal del MVA.
"""
b = [i for i in range(10)]
normal = [0.920, -0.302, 0.251]
histograma = 0
for k in range(len(escalas)):
    normales_intervalo = normales[:, k, :]
    angle = np.zeros(len(normales_intervalo))
    for i in range(len(normales_intervalo)):
        angle[i] = angulo(normal, normales_intervalo[i]) * 180 / np.pi
    hist, bins = np.histogram(angle, bins=10)
    histograma += hist

center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, histograma)
plt.xlabel("Angle (º)")
plt.ylabel("Counts")
plt.title("Fixed interval, varying time.\nAngle with respect to MVA normal.")
plt.show()

"""
Para hacer screenshots de tiempos particulares sólo tengo que elegir el valor de NN
y descomentar lo que sigue
"""
# NN = 28
# tiempo_central[NN]
# normales_t1 = normales[NN]
# angle = np.zeros((len(normales_t1), len(normales_t1)))
# thdec = hdec_to_UTC(tiempo_central[NN])
# for i in range(len(normales_t1)):
#     for j in range(len(normales_t1)):
#         angle[i, j] = angulo(normales_t1[i], normales_t1[j]) * 180 / np.pi
#
# plt.imshow(
#     angle,
#     aspect="auto",
#     origin="lower",
#     extent=(minimo, maximo, minimo, maximo),
#     cmap="viridis",
#     vmin=0,
#     vmax=10,
# )
# cbar = plt.colorbar()
# cbar.set_label("Angle between normals", rotation=270)
# plt.title(f"Heatmap centered in t = {thdec}")
# plt.xlabel("Radius (s)")
# plt.ylabel("Radius (s)")
# plt.show()
