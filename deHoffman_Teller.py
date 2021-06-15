"""Its aim is to find the frame velocity vector (HT velocity), V_HT , that best
agrees with the set of measured values of magnetic field, B, and electric field,
E, or magnetic field, B, and plasma bulk velocity, v, in cases where the convection
electric field, -v × B, can be used as a proxy for E."""

import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag_1s, importar_swia
from funciones import fechas, tiempos, diezmar

float_formatter = "{:.3g}".format
np.set_printoptions(formatter={"float_kind": float_formatter})

mu0 = 4e-7 * np.pi  # T m/A


def delta_ij(i, j):
    if i == j:
        delta = 1
    else:
        delta = 0

    return delta


def Kij(B, m):
    """Calcula la matriz Kij para el vector m-ésimo de campo B."""
    P_ij = np.zeros((3, 3))
    for i in range(3):  # para las tres coordenadas
        for j in range(3):
            P_ij[i, j] = (
                delta_ij(i, j) - B[m, i] * B[m, j] / np.linalg.norm(B[m, :]) ** 2
            )

    K_ij = np.linalg.norm(B[m, :]) ** 2 * P_ij
    return K_ij


year, month, day, doy = 2016, "03", 16, "076"  # fechas()
ti, tf = 18.2258, 18.235  # tiempos()
mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
swia, t_swia, density, temp, v = importar_swia(year, month, day, ti, tf)

# si hay buenos datos de v_electrones usamos ve en lugar de v


normal = [0.920, -0.302, 0.251]
# B y v no tienen la misma cantidad de puntos, tengo que diezmar (como siempre)

idx = diezmar(t, t_swia)
B_cut = B[idx]

K0 = np.mean(
    [Kij(B_cut, m) for m in range(len(B_cut))], axis=0
)  # tiene que ser una matriz
KK = np.mean([np.dot(Kij(B_cut, m), v[m, :]) for m in range(len(B_cut))], axis=0)

v_HT = np.dot(np.linalg.inv(K0), KK)
vn = np.dot(normal, v_HT)

# v alfvén
v_a = np.array(
    [(mu0 * density[i]) ** (-1 / 2) * B_cut[i, :] * 1e-9 for i in range(len(B_cut))]
)

v_prima = v - v_HT
print(f"v_HT = {v_HT}, v a lo largo de la normal = {vn:.3g}")

Ecv = np.cross(-v, B_cut)
EHT = np.cross(-v_HT, B_cut)

plt.figure()
plt.scatter(Ecv[:, 0], EHT[:, 0])
plt.scatter(Ecv[:, 1], EHT[:, 1])
plt.scatter(Ecv[:, 2], EHT[:, 2])
plt.title("Ecv vs EHT")
plt.xlabel("Ecv")
plt.ylabel("EHT")

# plt.figure()
# plt.scatter(v_a[:, 0], v_prima[:, 0])
# plt.scatter(v_a[:, 1], v_prima[:, 1])
# plt.scatter(v_a[:, 2], v_prima[:, 2])
# plt.title("walen")
# plt.xlabel("v alfven")
# plt.ylabel("v-vHT")
plt.show()
