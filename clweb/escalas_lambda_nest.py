import numpy as np
import time as time
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as md

sys.path.append("..")

from funciones import array_datenums, Mij, fechas, angulo, donde
from importar_datos import importar_mag
from funciones_plot import imshow_UTC, plot_datetime, scatter_datetime

"""
Calcula el cociente de lambdas y el ángulo entre la normal del MVA y la normal
obtenida centrándome en cada t para un radio fijo r = 15s. Superpone ambos.
"""
# year, month, day, doy = fechas()
# hora = input("Hora en HH\n")

year, month, day, doy = 2016, "03", 16, "076"
hora = 18

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

# ti = float(input('Tiempo inicial del barrido\n'))
# tf = float(input('Tiempo final del barrido\n'))
ti = tiempos[0]
tf = tiempos[3]
timestamps = array_datenums(year, month, day, tiempos)  # lo convierto a datenum


mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

tiempo_central = np.zeros(
    int((tf - ti) * 3600 * 32)
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

M = len(t)  # el numero de datos


program_starts = time.time()
cociente = np.zeros((len(tiempo_central), len(escalas)))
normales = np.zeros((len(tiempo_central), len(escalas), 3))

for i in range(len(tiempo_central)):
    for j in range(len(escalas)):
        inicio = donde(t, tiempo_central[i] - escalas[j])
        fin = donde(t, tiempo_central[i] + escalas[j])

        # ahora empieza el MVA con los datos que elegí
        B_cut = B[inicio : fin + 1, :]

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
        normales[i, j] = x3

program_ends = time.time()

print(f"El loop tardó {program_ends-program_starts:.2f} s")

"""
Quiero el ángulo contra la normal del MVA para cada tiempo t con escala r = 15s
"""
normal = [0.920, -0.302, 0.251]
rr = 17  # radio
normales_radio = normales[:, rr - 1, :]
angle = np.zeros(len(normales_radio))
for i in range(len(normales_radio)):
    angle[i] = angulo(normal, normales_radio[i]) * 180 / np.pi

tiempos_MVA = [18.2258, 18.235]
tMVA1 = donde(tiempo_central, tiempos_MVA[0])
tMVA2 = donde(tiempo_central, tiempos_MVA[1])
angulo_medio = np.mean(angle[tMVA1:tMVA2])
print(f"el ángulo medio es {angulo_medio}")

times = array_datenums(year, month, day, tiempo_central)
t_graph = md.date2num(times)
plt.figure()
plt.subplots_adjust(bottom=0.2)
plt.xticks(rotation=25)
ax = plt.gca()
xfmt = md.DateFormatter("%H:%M:%S")
ax.xaxis.set_major_formatter(xfmt)
plt.imshow(
    np.transpose(cociente),
    aspect="auto",
    origin="lower",
    extent=(t_graph[0], t_graph[-1], escalas[0] * 3600, escalas[-1] * 3600),
    cmap="viridis",
    vmin=0,
    vmax=20,
)
cbar = plt.colorbar()
scatter_datetime(year, month, day, tiempo_central, angle, "red")
for tt in [times[tMVA1], times[tMVA2]]:
    plt.axvline(x=tt, color="m")  # plotea los tiempos t1t2t3t4
plt.axhline(y=rr, color="k")
plt.title(rf"Heatmap of $\lambda2/\lambda3$ for the {day}-{month}-{year} crossing")
cbar.set_label(r"$\lambda2/\lambda3$", rotation=90)
plt.xlabel("Central Time (hh:mm:ss)")
plt.ylabel("Interval radius (s)")
sec_y = ax.secondary_yaxis("right", color="r")
sec_y.set_ylabel("Angle with respect to the MVA normal (º)", color="red")
plt.show()
