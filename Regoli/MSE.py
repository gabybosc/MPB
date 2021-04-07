import numpy as np
import matplotlib.pyplot as plt
import sys

path = "../../../datos/simulacion_leonardo/"
datos = np.load(
    path + "recorte_1RM.npy"
)  # todos los datos ordenados con y,z menores a 1 RM

sys.path.append("..")
from funciones import donde, angulo


"""
x y z Rho Ux Uy Uz Bx By Bz (0 a 9)
P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz (10 a 19)
O2pP OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz (20 a 29)
CO2pP b1x b1y b1z e jx jy jz (30 a 37)
"""

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
g = 3.7  # Mars surface gravity, m/s²

# vectores
posicion = datos[:, :3]
B = datos[:, 7:10]
J = datos[:, -3:]
b = datos[:, 31:34]

velocidad_plasma = datos[:, 4:7]  # km/s
velocidad_i = datos[:, 12:15]

# densidades
rho = datos[:, 3]
HpRho = datos[:, 11]
O2pRho = datos[:, 16]
OpRho = datos[:, 21]
CO2pRho = datos[:, 26]

# Presiones
Ptermica = datos[:, 10]  # nPa
HP = datos[:, 15]
O2P = datos[:, 20]
OP = datos[:, 25]
CO2P = datos[:, 30]

P_heavy = OP + O2P + CO2P
P_B = np.linalg.norm(B, axis=1) ** 2 * 1e-9 / (2 * mu0)
Pe = Ptermica - CO2P - HP - OP - O2P
P_ram = 1.67e-6 * HpRho * velocidad_i[:, 0] ** 2  # nPa

beta = Ptermica / P_B

beta_str = (Ptermica + P_ram) / P_B

rho_heavies = OpRho + O2pRho + CO2pRho
density_ratio = HpRho / rho_heavies
mass_ratio = HpRho / (OpRho * 16 + O2pRho * 32 + CO2pRho * 44)

upstream = donde(posicion[:, 0], 2.5)
v_sw = velocidad_i[upstream, :]
B_sw = B[upstream, :]

x_MSE = -v_sw / np.linalg.norm(v_sw)
z_MSE = np.cross(-v_sw, B_sw)
z_MSE = z_MSE / np.linalg.norm(z_MSE)
y_MSE = np.cross(z_MSE, x_MSE)
y_MSE = y_MSE / np.linalg.norm(y_MSE)
theta = angulo(z_MSE, [0, 0, 1])
# phi = angulo()

"""ahora que tengo el ángulo de la rotación entre MSO/MSE simplemente
tengo que rotar la grilla con tres matrices (yaw-pitch-roll)
| ca sa 0| | cb 0 sb| |1  0  0 |
|-sa ca 0| | 0  1  0| |0 cc -sc|
| 0  0  1| |-sb 0 cb| |0 sc  cc|

La matriz final queda:
|ca*cb   ca*sb*sc-sa*cc     ca*sb*cc+sa*sc|
|sa*cb   sa*sb*sc+ca*cc     sa*sb*cc-ca*sc|
| -sb          cb*sc             cb*cc     |
https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions

La matriz es para vectores en vertical
"""
beta = -np.arcsin(x_MSE[2])
alpha = np.arcsin(x_MSE[1] / np.cos(beta))
gamma = np.arccos(z_MSE[2] / np.cos(beta))

ca = np.cos(alpha)
sa = np.sin(alpha)
cb = np.cos(beta)
sb = np.sin(beta)
cg = np.cos(gamma)
sg = np.sin(gamma)
yaw = np.array([[ca, sa, 0], [-sa, ca, 0], [0, 0, 1]])
pitch = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])
roll = np.array([[1, 0, 0], [0, cg, -sg], [0, sg, cg]])

rot_matrix = np.linalg.multi_dot([yaw, pitch, roll])
