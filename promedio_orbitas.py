import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag
from funciones import donde

mu0 = 4 * np.pi * 1e-7


mag, t, B, posicion = importar_mag(2016, "03", 16, 0.1, 23.9)
Bcut = B[400000:]
Bnorm = np.linalg.norm(Bcut, axis=1)
pos = np.linalg.norm(posicion[400000:], axis=1)
# puntos = int(len(Bnorm) / 7)

plt.plot(Bnorm[:512000])
plt.show()

puntos = 512000 + 5908
plt.figure()
for i in range(6):
    plt.plot(pos[int(i * puntos) : int((i + 1) * puntos)])
    # están perfectamente superpuestas las posiciones

plt.figure()
for i in range(6):
    plt.plot(Bnorm[int(i * puntos) : int((i + 1) * puntos)])
    plt.legend(["1", "2", "3", "4", "5"])
plt.show()

lista = np.array(
    [Bnorm[int(i * puntos) : int((i + 1) * puntos - 200000)] for i in range(4)]
)
prom = np.mean(lista, axis=0)

for i in range(6):
    plt.plot(Bnorm[int(i * puntos) : int((i + 1) * puntos)], linewidth=0.5)
plt.plot(prom, linewidth=2, label="prom")
plt.legend()
plt.show()

# ahora lo mismo pero para las componentes
comps = np.array(
    [Bcut[int(i * puntos) : int((i + 1) * puntos - 200000), :] for i in range(4)]
)
medio = np.mean(comps, axis=0)

plt.figure()
for i in range(6):
    plt.plot(Bcut[int(i * puntos) : int((i + 1) * puntos), 0], linewidth=0.5)

plt.plot(t[:317908], medio)
plt.xlabel("t")
plt.ylabel("componentes de B promediadas")
plt.show()

datos_On = np.loadtxt(
    "../../../../media/gabybosc/datos/Chuanfei/sat_trajectory_HallOn.sat", skiprows=2
)
datos = datos_On[
    1000:1250
]  # le saco los primeros mil puntos porque son el lado de atrás que no me interesa

x = datos[:, 8]
y = datos[:, 9]
z = datos[:, 10]
B = datos[:, 15:18]
b1 = datos[:, 40:43]

velocidad_plasma = datos[:, 12:15]
velocidad_H = datos[:, 21:24]
J = datos[:, -3:]

P_plasma = datos[:, 19]

posicion = datos[:, 8:11]
presion = {
    "e-": datos[:, 18],
    "H+": datos[:, 24],
    "O+": datos[:, 34],
    "O2+": datos[:, 29],
    "CO2+": datos[:, 39],
}
densidad = {
    "e-": datos[:, 11],
    "H+": datos[:, 20],
    "O+": datos[:, 30],
    "O2+": datos[:, 25],
    "CO2+": datos[:, 35],
}

P_heavy = presion["O+"] + presion["O2+"] + presion["CO2+"]
P_B = np.linalg.norm(B, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram = 1.67e-6 * densidad["H+"] * velocidad_H[:, 0] ** 2  # nPa
P_total = P_heavy + P_B + presion["H+"] + P_ram + presion["e-"]

beta = P_plasma / P_B

density_mean = [
    np.mean(densidad["H+"][i : i + 50]) for i in range(len(densidad["H+"]) - 50)
]
ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km
ion_length = np.append(ion_length, ion_length[-51:-1])

t1t2t3t4 = [18.2167, 18.2204, 18.235, 18.2476]
i1 = donde(t, t1t2t3t4[0])
i2 = donde(t, t1t2t3t4[1])
i3 = donde(t, t1t2t3t4[2])
i4 = donde(t, t1t2t3t4[3])

limites_simu = [
    0.8978406562499996,
    1.0030125729166666,
    1.3053818333333331,
    1.3974072604166663,
]
# las posiciones de entrada y salida de la MPB en el eje x en RM (de la simu)
x1, x2, x3, x4 = limites_simu[0], limites_simu[1], limites_simu[2], limites_simu[3]

pos_MPB_simu = np.mean(limites_simu)
