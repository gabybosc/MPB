"""Its aim is to find the frame velocity vector (HT velocity), V_HT , that best
agrees with the set of measured values of magnetic field, B, and electric field,
E, or magnetic field, B, and plasma bulk velocity, v, in cases where the convection
electric field, -v × B, can be used as a proxy for E."""

import numpy as np

t = [
    0.00,
    4.38,
    8.75,
    13.12,
    17.49,
    21.87,
    26.24,
    30.61,
    34.98,
    39.36,
    43.73,
    48.10,
    52.47,
    56.85,
    61.22,
    65.59,
]
Bx = [
    -13.6,
    -14.8,
    -13.4,
    -14.0,
    -7.1,
    -0.9,
    -10.0,
    -6.1,
    1.2,
    -3.4,
    -0.9,
    -1.0,
    11.0,
    19.1,
    24.9,
    29.2,
]
By = [
    -24.7,
    -24.9,
    -17.2,
    -25.0,
    -4.5,
    -5.0,
    -0.4,
    -4.8,
    1.6,
    -3.9,
    1.2,
    -1.5,
    13.2,
    34.4,
    50.1,
    47.1,
]
Bz = [
    54.6,
    58.7,
    62.4,
    43.8,
    33.5,
    44.4,
    44.6,
    21.1,
    21.0,
    4.1,
    5.0,
    12.3,
    29.7,
    20.1,
    1.9,
    -10.6,
]
vx = [
    -111.0,
    -102.0,
    -111.0,
    -135.0,
    -128.0,
    -82.0,
    -139.0,
    -143.0,
    -132.0,
    -117.0,
    -112.0,
    -98.0,
    -90.0,
    -83.0,
    -82.0,
    -93.0,
]
vy = [
    -211.0,
    -213.0,
    -196.0,
    -229.0,
    -252.0,
    -237.0,
    -228.0,
    -241.0,
    -226.0,
    -217.0,
    -210.0,
    -212.0,
    -186.0,
    -168.0,
    -129.0,
    -123.0,
]
vz = [
    57.0,
    41.0,
    34.0,
    54.0,
    54.0,
    51.0,
    77.0,
    57.0,
    80.0,
    79.0,
    93.0,
    92.0,
    104.0,
    121.0,
    88.0,
    53.0,
]
N = [
    12.18,
    9.22,
    9.92,
    18.08,
    20.39,
    15.00,
    20.19,
    23.53,
    24.31,
    25.91,
    26.17,
    24.49,
    22.20,
    22.86,
    17.56,
    17.86,
]


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


B = np.transpose([Bx, By, Bz])
v = np.transpose([vx, vy, vz])

K0 = np.mean([Kij(B, m) for m in range(len(B))], axis=0)  # tiene que ser una matriz
KK = np.mean([np.dot(Kij(B, m), v[m, :]) for m in range(len(B))], axis=0)

v_HT = np.dot(np.linalg.inv(K0), KK)

print(v_HT)
