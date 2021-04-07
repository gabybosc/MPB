import numpy as np
import matplotlib.pyplot as plt

# import sys


path = "../../../datos/simulacion_leonardo/"
reordenados = np.load(
    path + "recorte_zcero.npy"
)  # todos los datos ordenados con y,z menores a 1 RM

# sys.path.append("..")
# from funciones import donde


"""
x y z Rho Ux Uy Uz Bx By Bz (0 a 9)
P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz (10 a 19)
O2pP OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz (20 a 29)
CO2pP b1x b1y b1z e jx jy jz (30 a 37)
"""

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
g = 3.7  # Mars surface gravity, m/s²

x = reordenados[:, 0]
y = reordenados[:, 1]
z = reordenados[:, 2]
rho = reordenados[:, 3]
B = reordenados[:, 7:10]
J = reordenados[:, -3:]
b = reordenados[:, 31:34]

velocidad_plasma = reordenados[:, 4:7]  # km/s
velocidad_i = reordenados[:, 12:15]

HpRho = reordenados[:, 11]
O2pRho = reordenados[:, 16]
OpRho = reordenados[:, 21]
CO2pRho = reordenados[:, 26]

# Presiones
Ptermica = reordenados[:, 10]  # nPa
HP = reordenados[:, 15]
O2P = reordenados[:, 20]
OP = reordenados[:, 25]
CO2P = reordenados[:, 30]

P_heavy = OP + O2P + CO2P
P_B = np.linalg.norm(B, axis=1) ** 2 * 1e-9 / (2 * mu0)
Pe = Ptermica - CO2P - HP - OP - O2P
P_ram = 1.67e-6 * HpRho * velocidad_i[:, 0] ** 2  # nPa

beta = Ptermica / P_B

beta_str = (Ptermica + P_ram) / P_B

rho_heavies = OpRho + O2pRho + CO2pRho
density_ratio = HpRho / rho_heavies
mass_ratio = HpRho / (OpRho * 16 + O2pRho * 32 + CO2pRho * 44)


# MPB
L = 0.96
x0 = 0.78
e = 0.9

R = [1.082, -0.064, 0.515]
theta = np.linspace(0, np.pi * 2, 100)

r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(theta))

X1 = x0 + r1 * np.cos(theta)
Y1 = r1 * np.sin(theta)

# BS

L_BS = 2.04
x0_BS = 0.64
e_BS = 1.03

r0_BS = 1.64 - x0_BS
theta0_BS = np.arccos(r0_BS)

L0_BS = np.linalg.norm(r0_BS) * (1 + e_BS * np.cos(theta0_BS))
r1_BS = L0_BS / (1 + e_BS * np.cos(theta))

X_BS = x0_BS + r1_BS * np.cos(theta)
Y_BS = r1_BS * np.sin(theta)

# Campo E
e_SI = 1.6e-19  # C

v_SI = velocidad_plasma * 1e3  # m/s
B_SI = B * 1e-9  # T
n_SI = HpRho * 1e6  # 1/m3
J_SI = J * 1e-6  # A/m2

Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)

# Campo B y Corrientes

fig, axs = plt.subplots(2, 3)
axs[0, 0].scatter(
    x, y, c=B[:, 0], vmin=-20, vmax=20, s=35, cmap="coolwarm",
)
axs[0, 0].plot(X_BS, Y_BS, c="C4", label="BS")
axs[0, 0].plot(X1, Y1, label="MPB", c="C2")
axs[0, 0].set_ylabel("Y (RM)")
axs[0, 0].set_title("Z=0, Bx")
axs[0, 0].set_xlim([0, 2])
axs[0, 0].set_ylim([-1, 1])
plt.setp(axs[0, 0].get_xticklabels(), visible=False)

sc = axs[0, 1].scatter(x, y, c=B[:, 1], vmin=-20, vmax=20, s=35, cmap="coolwarm",)
axs[0, 1].plot(
    X_BS, Y_BS, c="C4",
)
axs[0, 1].plot(X1, Y1, c="C2")
axs[0, 1].set_xlim([0, 2])
axs[0, 1].set_ylim([-1, 1])
axs[0, 1].set_title("Z=0, By")
plt.setp(axs[0, 1].get_yticklabels(), visible=False)
plt.setp(axs[0, 1].get_xticklabels(), visible=False)

sc = axs[0, 2].scatter(x, y, c=B[:, 2], vmin=-20, vmax=20, s=35, cmap="coolwarm",)
axs[0, 2].plot(X_BS, Y_BS, c="C4")
axs[0, 2].plot(X1, Y1, c="C2")
axs[0, 2].set_xlim([0, 2])
axs[0, 2].set_ylim([-1, 1])
axs[0, 2].set_title("Z=0, Bz")
plt.setp(axs[0, 2].get_yticklabels(), visible=False)
plt.setp(axs[0, 2].get_xticklabels(), visible=False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.55, 0.05, 0.38])  # [left, bottom, width, height]
fig.colorbar(sc, cax=cbar_ax)

axs[1, 0].scatter(
    x, y, c=J[:, 0], vmin=-0.05, vmax=0.05, s=35, cmap="coolwarm",
)
axs[1, 0].plot(X_BS, Y_BS, c="C4")
axs[1, 0].plot(X1, Y1, c="C2")
axs[1, 0].set_xlabel("X (RM)")
axs[1, 0].set_ylabel("Y (RM)")
axs[1, 0].set_title("Z=0, Jx")
axs[1, 0].set_xlim([0, 2])
axs[1, 0].set_ylim([-1, 1])

sc = axs[1, 1].scatter(x, y, c=J[:, 1], vmin=-0.05, vmax=0.05, s=35, cmap="coolwarm",)
axs[1, 1].plot(X_BS, Y_BS, c="C4")
axs[1, 1].plot(X1, Y1, c="C2")
axs[1, 1].set_xlim([0, 2])
axs[1, 1].set_ylim([-1, 1])
axs[1, 1].set_xlabel("X (RM)")
axs[1, 1].set_title("Z=0, Jy")
plt.setp(axs[1, 1].get_yticklabels(), visible=False)

sc = axs[1, 2].scatter(x, y, c=J[:, 2], vmin=-0.05, vmax=0.05, s=35, cmap="coolwarm",)
axs[1, 2].plot(X_BS, Y_BS, label="BS", c="C4")
axs[1, 2].plot(X1, Y1, label="MPB", c="C2")
axs[1, 2].set_xlim([0, 2])
axs[1, 2].set_ylim([-1, 1])
axs[1, 2].set_xlabel("X (RM)")
axs[1, 2].set_title("Z=0, Jz")
axs[1, 2].legend()
plt.setp(axs[1, 2].get_yticklabels(), visible=False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.05, 0.05, 0.38])
fig.colorbar(sc, cax=cbar_ax)


plt.show(block=False)

# Campo Ecv y Ehall

fig, axs = plt.subplots(2, 3)
axs[0, 0].scatter(
    x, y, c=Ecv[:, 0] * 1e3, vmin=-1, vmax=1, s=35, cmap="coolwarm",
)
axs[0, 0].plot(X_BS, Y_BS, c="C4", label="BS")
axs[0, 0].plot(X1, Y1, label="MPB", c="C2")
axs[0, 0].set_ylabel("Y (RM)")
axs[0, 0].set_title("Z=0, Ecv_x")
axs[0, 0].set_xlim([0, 2])
axs[0, 0].set_ylim([-1, 1])
plt.setp(axs[0, 0].get_xticklabels(), visible=False)

axs[0, 1].scatter(
    x, y, c=Ecv[:, 1] * 1e3, vmin=-1, vmax=1, s=35, cmap="coolwarm",
)
axs[0, 1].plot(
    X_BS, Y_BS, c="C4",
)
axs[0, 1].plot(X1, Y1, c="C2")
axs[0, 1].set_xlim([0, 2])
axs[0, 1].set_ylim([-1, 1])
axs[0, 1].set_title("Z=0, Ecv_y")
plt.setp(axs[0, 1].get_yticklabels(), visible=False)
plt.setp(axs[0, 1].get_xticklabels(), visible=False)

sc = axs[0, 2].scatter(x, y, c=Ecv[:, 2] * 1e3, vmin=-1, vmax=1, s=35, cmap="coolwarm",)
axs[0, 2].plot(X_BS, Y_BS, c="C4")
axs[0, 2].plot(X1, Y1, c="C2")
axs[0, 2].set_xlim([0, 2])
axs[0, 2].set_ylim([-1, 1])
axs[0, 2].set_title("Z=0, Ecv_z")
plt.setp(axs[0, 2].get_yticklabels(), visible=False)
plt.setp(axs[0, 2].get_xticklabels(), visible=False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.55, 0.05, 0.38])  # [left, bottom, width, height]
fig.colorbar(sc, cax=cbar_ax)

axs[1, 0].scatter(
    x, y, c=Ehall[:, 0] * 1e3, vmin=-1, vmax=1, s=35, cmap="coolwarm",
)
axs[1, 0].plot(X_BS, Y_BS, c="C4")
axs[1, 0].plot(X1, Y1, c="C2")
axs[1, 0].set_xlabel("X (RM)")
axs[1, 0].set_ylabel("Y (RM)")
axs[1, 0].set_title("Z=0, Ehall_x")
axs[1, 0].set_xlim([0, 2])
axs[1, 0].set_ylim([-1, 1])

axs[1, 1].scatter(
    x, y, c=Ehall[:, 1] * 1e3, vmin=-1, vmax=1, s=35, cmap="coolwarm",
)
axs[1, 1].plot(X_BS, Y_BS, c="C4")
axs[1, 1].plot(X1, Y1, c="C2")
axs[1, 1].set_xlim([0, 2])
axs[1, 1].set_ylim([-1, 1])
axs[1, 1].set_xlabel("X (RM)")
axs[1, 1].set_title("Z=0, Ehall_y")
plt.setp(axs[1, 1].get_yticklabels(), visible=False)

sc = axs[1, 2].scatter(
    x, y, c=Ehall[:, 2] * 1e3, vmin=-1, vmax=1, s=35, cmap="coolwarm",
)
axs[1, 2].plot(X_BS, Y_BS, label="BS", c="C4")
axs[1, 2].plot(X1, Y1, label="MPB", c="C2")
axs[1, 2].set_xlim([0, 2])
axs[1, 2].set_ylim([-1, 1])
axs[1, 2].set_xlabel("X (RM)")
axs[1, 2].set_title("Z=0, Ehall_z")
axs[1, 2].legend()
plt.setp(axs[1, 2].get_yticklabels(), visible=False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.05, 0.05, 0.38])
fig.colorbar(sc, cax=cbar_ax)


plt.show(block=False)

# qué es el b1
fig, axs = plt.subplots(3, 3)
axs[0, 0].scatter(
    x, y, c=B[:, 0], vmin=-20, vmax=20, s=35, cmap="coolwarm",
)
axs[0, 0].plot(X_BS, Y_BS, c="C4", label="BS")
axs[0, 0].plot(X1, Y1, label="MPB", c="C2")
axs[0, 0].set_ylabel("Y (RM)")
axs[0, 0].set_title("Z=0, Bx", fontsize=8)
axs[0, 0].set_xlim([0, 2])
axs[0, 0].set_ylim([-1, 1])
plt.setp(axs[0, 0].get_xticklabels(), visible=False)

sc = axs[0, 1].scatter(x, y, c=B[:, 1], vmin=-20, vmax=20, s=35, cmap="coolwarm",)
axs[0, 1].plot(
    X_BS, Y_BS, c="C4",
)
axs[0, 1].plot(X1, Y1, c="C2")
axs[0, 1].set_xlim([0, 2])
axs[0, 1].set_ylim([-1, 1])
axs[0, 1].set_title("Z=0, By", fontsize=8)
plt.setp(axs[0, 1].get_yticklabels(), visible=False)
plt.setp(axs[0, 1].get_xticklabels(), visible=False)

sc = axs[0, 2].scatter(x, y, c=B[:, 2], vmin=-20, vmax=20, s=35, cmap="coolwarm",)
axs[0, 2].plot(X_BS, Y_BS, c="C4")
axs[0, 2].plot(X1, Y1, c="C2")
axs[0, 2].set_xlim([0, 2])
axs[0, 2].set_ylim([-1, 1])
axs[0, 2].set_title("Z=0, Bz", fontsize=8)
plt.setp(axs[0, 2].get_yticklabels(), visible=False)
plt.setp(axs[0, 2].get_xticklabels(), visible=False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])  # [left, bottom, width, height]
fig.colorbar(sc, cax=cbar_ax)

axs[1, 0].scatter(x, y, c=b[:, 0], vmin=-20, vmax=20, s=35, cmap="coolwarm")
axs[1, 0].plot(X_BS, Y_BS, c="C4")
axs[1, 0].plot(X1, Y1, c="C2")
axs[1, 0].set_ylabel("Y (RM)")
axs[1, 0].set_title("Z=0, b1x", fontsize=8)
axs[1, 0].set_xlim([0, 2])
axs[1, 0].set_ylim([-1, 1])
plt.setp(axs[1, 0].get_xticklabels(), visible=False)

sc = axs[1, 1].scatter(x, y, c=b[:, 1], vmin=-20, vmax=20, s=35, cmap="coolwarm",)
axs[1, 1].plot(X_BS, Y_BS, c="C4")
axs[1, 1].plot(X1, Y1, c="C2")
axs[1, 1].set_xlim([0, 2])
axs[1, 1].set_ylim([-1, 1])
axs[1, 1].set_title("Z=0, b1y", fontsize=8)
plt.setp(axs[1, 1].get_yticklabels(), visible=False)
plt.setp(axs[1, 1].get_xticklabels(), visible=False)

sc = axs[1, 2].scatter(x, y, c=b[:, 2], vmin=-20, vmax=20, s=35, cmap="coolwarm",)
axs[1, 2].plot(X_BS, Y_BS, label="BS", c="C4")
axs[1, 2].plot(X1, Y1, label="MPB", c="C2")
axs[1, 2].set_xlim([0, 2])
axs[1, 2].set_ylim([-1, 1])
axs[1, 2].set_title("Z=0, b1z", fontsize=8)
plt.setp(axs[1, 2].get_yticklabels(), visible=False)
plt.setp(axs[1, 2].get_xticklabels(), visible=False)
# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.05, 0.05, 0.38])
# fig.colorbar(sc, cax=cbar_ax)

axs[2, 0].scatter(x, y, c=B[:, 0] - b[:, 0], vmin=-20, vmax=20, s=35, cmap="coolwarm")
axs[2, 0].plot(X_BS, Y_BS, c="C4")
axs[2, 0].plot(X1, Y1, c="C2")
axs[2, 0].set_xlabel("X (RM)")
axs[2, 0].set_ylabel("Y (RM)")
axs[2, 0].set_title("Z=0, Bx-b1x", fontsize=8)
axs[2, 0].set_xlim([0, 2])
axs[2, 0].set_ylim([-1, 1])

sc = axs[2, 1].scatter(
    x, y, c=B[:, 1] - b[:, 1], vmin=-20, vmax=20, s=35, cmap="coolwarm",
)
axs[2, 1].plot(X_BS, Y_BS, c="C4")
axs[2, 1].plot(X1, Y1, c="C2")
axs[2, 1].set_xlim([0, 2])
axs[2, 1].set_ylim([-1, 1])
axs[2, 1].set_xlabel("X (RM)")
axs[2, 1].set_title("Z=0, By-b1y", fontsize=8)
plt.setp(axs[2, 1].get_yticklabels(), visible=False)

sc = axs[2, 2].scatter(
    x, y, c=B[:, 2] - b[:, 2], vmin=-20, vmax=20, s=35, cmap="coolwarm",
)
axs[2, 2].plot(X_BS, Y_BS, label="BS", c="C4")
axs[2, 2].plot(X1, Y1, label="MPB", c="C2")
axs[2, 2].set_xlim([0, 2])
axs[2, 2].set_ylim([-1, 1])
axs[2, 2].set_xlabel("X (RM)")
axs[2, 2].set_title("Z=0, Bz-b1z", fontsize=8)
plt.setp(axs[2, 2].get_yticklabels(), visible=False)

plt.show(block=False)
