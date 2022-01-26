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
e_SI = 1.6e-19  # C

# vectores
posicion = datos[:, :3]
B = datos[:, 7:10]
J = datos[:, -3:]  # uA/m²
b = datos[:, 31:34]

velocidad_plasma = datos[:, 4:7]  # km/s
velocidad_i = datos[:, 12:15]

# densidades
densidad = datos[:, 3]
densidad_H = datos[:, 11]
densidad_O2 = datos[:, 16]
densidad_O = datos[:, 21]
densidad_CO2 = datos[:, 26]

# Presiones
Ptermica = datos[:, 10]  # nPa
presion_H = datos[:, 15]
presion_O2 = datos[:, 20]
presion_O = datos[:, 25]
presion_CO2 = datos[:, 30]

P_heavy = presion_O + presion_O2 + presion_CO2
P_mag = np.linalg.norm(B, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_e = Ptermica - P_heavy - presion_H
P_dyn = 1.67e-6 * densidad_H * velocidad_i[:, 0] ** 2  # nPa

n_SI = densidad_H * 1e6  # 1/m3
delta_v = np.array([J[i, :] * 1e-6 / (e_SI * n_SI[i]) for i in range(len(J))])  # m/s
velocidad_e = velocidad_i - delta_v * 1e-3

beta = Ptermica / P_mag

beta_str = (Ptermica + P_dyn) / P_mag

densidad_heavies = densidad_O + densidad_O2 + densidad_CO2
density_ratio = densidad_H / densidad_heavies
mass_ratio = densidad_H / (densidad_O * 16 + densidad_O2 * 32 + densidad_CO2 * 44)

x_MSE = np.array([1, 0, 0])  # esto se puede después cambiar a que sea v_sw
upstream = donde(posicion[:, 0], 2.5)
B_sw = B[upstream, :]
z_MSE = np.cross(x_MSE, B_sw)
z_MSE = z_MSE / np.linalg.norm(z_MSE)
y_MSE = np.cross(z_MSE, x_MSE)
y_MSE = y_MSE / np.linalg.norm(y_MSE)
theta = angulo(y_MSE, [0, 1, 0])

"""ahora que tengo el ángulo de la rotación entre MSO/MSE simplemente
tengo que rotar la grilla con la matriz
1 0 0
0 c s
0 -s c
https://en.wikipedia.org/wiki/Rotation_of_axes
"""

matriz = np.array(
    [[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]]
)

# Descomentar la siguiente parte para chequear que MSE es lo que debería
# MSO = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
# MSE = np.dot(matriz, MSO)

posicion_MSE = np.zeros((len(datos), 3))
B_MSE = np.zeros((len(datos), 3))
J_MSE = np.zeros((len(datos), 3))
b1_MSE = np.zeros((len(datos), 3))
v_MSE = np.zeros((len(datos), 3))
for i in range(len(posicion)):
    posicion_MSE[i, :] = np.dot(matriz, posicion[i, :])
    B_MSE[i, :] = np.dot(matriz, B[i, :])
    J_MSE[i, :] = np.dot(matriz, J[i, :])
    b1_MSE[i, :] = np.dot(matriz, b[i, :])
    v_MSE[i, :] = np.dot(matriz, velocidad_plasma[i, :])


v_SI = v_MSE * 1e3  # m/s
B_SI = b1_MSE * 1e-9  # T
n_SI = densidad_H * 1e6  # 1/m3
J_SI = J_MSE * 1e-6  # A/m2

Ecv_MSE = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B_MSE))])  # V/m
Ehall_MSE = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B_MSE))]
)

# corte en z_MSE=0

a = [np.abs(posicion_MSE[i, 2]) <= 0.05 for i in range(len(posicion_MSE))]
idx = [i for i, x in enumerate(a) if x]

x = posicion_MSE[idx, 0]
y = posicion_MSE[idx, 1]

coordenadas = {0: "x", 1: "y", 2: "z"}

Ehall_norma = np.linalg.norm(Ehall_MSE[idx, :], axis=1)
Ecv_norma = np.linalg.norm(Ecv_MSE[idx, :], axis=1)

# Campo B y Corrientes


fig, axs = plt.subplots(2, 3)

for i in range(3):
    sc = axs[0, i].scatter(
        x, y, c=B_MSE[idx, i], vmin=-75, vmax=75, s=35, cmap="coolwarm",
    )
    axs[0, i].set_xlim([0, 2])
    axs[0, i].set_ylim([-1, 1])
    axs[0, i].set_title(f"Z=0, B{coordenadas[i]} ")
    plt.setp(axs[0, i].get_xticklabels(), visible=False)
    if i > 0:
        plt.setp(axs[0, i].get_yticklabels(), visible=False)
    else:
        axs[0, i].set_ylabel("Y MSE (RM)")

    sc2 = axs[1, i].scatter(
        x, y, c=J_MSE[idx, i], vmin=-0.25, vmax=0.25, s=35, cmap="coolwarm",
    )
    axs[1, i].set_xlabel("X MSE (RM)")
    axs[1, i].set_title(f"Z=0, J{coordenadas[i]}")
    axs[1, i].set_xlim([0, 2])
    axs[1, i].set_ylim([-1, 1])
    if i > 0:
        plt.setp(axs[1, i].get_yticklabels(), visible=False)
    else:
        axs[1, i].set_ylabel("Y MSE (RM)")

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.55, 0.05, 0.38])  # [left, bottom, width, height]
fig.colorbar(sc, cax=cbar_ax)

fig.subplots_adjust(right=0.8)
cbar_ax2 = fig.add_axes([0.85, 0.05, 0.05, 0.38])
fig.colorbar(sc2, cax=cbar_ax2)


plt.show(block=False)

# Campo eléctrico

fig, axs = plt.subplots(2, 3)
for i in range(3):
    sc = axs[0, i].scatter(
        x, y, c=Ecv_MSE[idx, i] * 1e3, s=35, vmin=-2.5, vmax=2.5, cmap="coolwarm",
    )
    axs[0, i].set_xlim([0, 2])
    axs[0, i].set_ylim([-1, 1])
    axs[0, i].set_title(f"Z=0, Ecv_{coordenadas[i]}")
    plt.setp(axs[0, i].get_xticklabels(), visible=False)
    if i > 0:
        plt.setp(axs[0, i].get_yticklabels(), visible=False)
    else:
        axs[0, i].set_ylabel("Y MSE (RM)")

    sc2 = axs[1, i].scatter(
        x, y, c=Ehall_MSE[idx, i] * 1e3, s=35, vmin=-5, vmax=5, cmap="coolwarm",
    )
    axs[1, i].set_xlabel("X MSE (RM)")
    axs[1, i].set_title(f"Z=0, Ehall_{coordenadas[i]}")
    axs[1, i].set_xlim([0, 2])
    axs[1, i].set_ylim([-1, 1])
    if i > 0:
        plt.setp(axs[1, i].get_yticklabels(), visible=False)
    else:
        axs[1, i].set_ylabel("Y MSE (RM)")

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.55, 0.05, 0.38])  # [left, bottom, width, height]
fig.colorbar(sc, cax=cbar_ax)

fig.subplots_adjust(right=0.8)
cbar_ax2 = fig.add_axes([0.85, 0.05, 0.05, 0.38])
fig.colorbar(sc2, cax=cbar_ax2)

plt.show(block=False)

# Norma de Ecv/Ehall

plt.figure()
plt.scatter(
    x, y, c=np.log10(Ecv_norma / Ehall_norma), vmin=-9, vmax=9, s=35, cmap="coolwarm",
)
plt.xlim([0, 2])
plt.ylim([-1, 1])
plt.title("log(|Ecv|/|Ehall|) en Z=0")
plt.xlabel("X MSE (Rm)")
plt.ylabel("Y MSE (Rm)")
plt.colorbar()
plt.show(block=False)

# Presión electrónica

plt.figure()
plt.scatter(x, y, c=P_e[idx])
plt.title("P_e en Z=0")
plt.xlabel("X MSE (Rm)")
plt.ylabel("Y MSE (Rm)")
plt.colorbar()
plt.show()
