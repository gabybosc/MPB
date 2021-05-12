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

mach = P_ram / Ptermica

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

# longitud inercial de protones en z=0

density_mean = [np.mean(HpRho[i : i + 50]) for i in range(len(HpRho) - 50)]

ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km

plt.figure()
plt.plot(x[:-50], ion_length)
plt.xlabel("X (RM)")
plt.ylabel("proton inertial length")
plt.title("Z = 0")
plt.show(block=False)
#  Heatmap de presiones

fig, axs = plt.subplots(1, 3)
axs[0].scatter(
    x, y, c=np.log10(P_B * 1e-9), vmin=-12, vmax=-9, s=35, cmap="inferno",
)
axs[0].plot(X_BS, Y_BS, label="BS")
axs[0].plot(X1, Y1, label="MPB", c="C2")
axs[0].set_xlabel("X (RM)")
axs[0].set_ylabel("Y (RM)")
axs[0].set_title("Z=0, log10(Pmag)")
axs[0].set_xlim([0, 2])
axs[0].set_ylim([-1, 1])

sc = axs[1].scatter(
    x, y, c=np.log10(P_ram * 1e-9), vmin=-12, vmax=-9, s=35, cmap="inferno",
)
axs[1].plot(X_BS, Y_BS, label="BS")
axs[1].plot(X1, Y1, label="MPB", c="C2")
axs[1].set_xlim([0, 2])
axs[1].set_ylim([-1, 1])
axs[1].set_xlabel("X (RM)")
# axs[1].set_ylabel("Y (RM)")
axs[1].set_title("Z=0, log10(Pdyn)")
axs[1].legend()
plt.setp(axs[1].get_yticklabels(), visible=False)

sc = axs[2].scatter(
    x, y, c=np.log10(Ptermica * 1e-9), vmin=-12, vmax=-9, s=35, cmap="inferno",
)
axs[2].plot(X_BS, Y_BS, label="BS")
axs[2].plot(X1, Y1, label="MPB", c="C2")
axs[2].set_xlim([0, 2])
axs[2].set_ylim([-1, 1])
axs[2].set_xlabel("X (RM)")
axs[2].set_title("Z=0, log10(Pth)")
axs[2].legend()
plt.setp(axs[2].get_yticklabels(), visible=False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(sc, cax=cbar_ax)
# plt.colorbar(sc)
plt.show(block=False)


# Beta

fig, axs = plt.subplots(1, 2)
axs[0].scatter(
    x, y, c=np.log10(beta), vmin=-4, vmax=4, s=35, cmap="coolwarm",
)
axs[0].plot(X_BS, Y_BS, label="BS")
axs[0].plot(X1, Y1, label="MPB", c="C2")
axs[0].set_xlabel("X (RM)")
axs[0].set_ylabel("Y (RM)")
axs[0].set_title("Z=0, log10(beta)")
axs[0].set_xlim([0, 2])
axs[0].set_ylim([-1, 1])

sc = axs[1].scatter(x, y, c=np.log10(beta_str), vmin=-4, vmax=4, s=35, cmap="coolwarm",)
axs[1].plot(X_BS, Y_BS, label="BS")
axs[1].plot(X1, Y1, label="MPB", c="C2")
axs[1].set_xlim([0, 2])
axs[1].set_ylim([-1, 1])
axs[1].set_xlabel("X (RM)")
# axs[1].set_ylabel("Y (RM)")
axs[1].set_title("Z=0, log10(beta*)")
axs[1].legend()
plt.setp(axs[1].get_yticklabels(), visible=False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(sc, cax=cbar_ax)
# plt.colorbar(sc)
plt.show(block=False)

# Densidades

fig, axs = plt.subplots(1, 2)
axs[0].scatter(
    x, y, c=np.log10(HpRho), vmin=-1, vmax=2, s=35, cmap="inferno",
)
axs[0].plot(X_BS, Y_BS, label="BS")
axs[0].plot(X1, Y1, label="MPB", c="C2")
axs[0].set_xlabel("X (RM)")
axs[0].set_ylabel("Y (RM)")
axs[0].set_title("Z=0, log10(Rho_H)")
axs[0].set_xlim([0, 2])
axs[0].set_ylim([-1, 1])

axs[1].scatter(
    x, y, c=np.log10(rho_heavies), vmin=-1, vmax=2, s=35, cmap="inferno",
)
axs[1].plot(X_BS, Y_BS, label="BS")
axs[1].plot(X1, Y1, label="MPB", c="C2")
axs[1].set_xlim([0, 2])
axs[1].set_ylim([-1, 1])
axs[1].set_xlabel("X (RM)")
axs[1].set_title("Z=0, log10(Rho_heavies)")
axs[1].legend()
plt.setp(axs[1].get_yticklabels(), visible=False)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(sc, cax=cbar_ax)
# plt.colorbar(sc)
plt.show(block=False)


# cociente de densidades


fig, axs = plt.subplots(1, 2)
axs[0].scatter(
    x, y, c=np.log10(density_ratio), vmin=-4, vmax=4, s=35, cmap="coolwarm",
)
axs[0].plot(X_BS, Y_BS, label="BS")
axs[0].plot(X1, Y1, label="MPB", c="C2")
axs[0].set_xlabel("X (RM)")
axs[0].set_ylabel("Y (RM)")
axs[0].set_title("Z=0, log10(cociente densidad masa)")
axs[0].set_xlim([0, 2])
axs[0].set_ylim([-1, 1])

sc = axs[1].scatter(
    x, y, c=np.log10(mass_ratio), vmin=-4, vmax=4, s=35, cmap="coolwarm",
)
axs[1].plot(X_BS, Y_BS, label="BS")
axs[1].plot(X1, Y1, label="MPB", c="C2")
axs[1].set_xlim([0, 2])
axs[1].set_ylim([-1, 1])
axs[1].set_xlabel("X (RM)")
# axs[1].set_ylabel("Y (RM)")
axs[1].set_title("Z=0, log10(cociente densidad numérica)")
axs[1].legend()
plt.setp(axs[1].get_yticklabels(), visible=False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(sc, cax=cbar_ax)
# plt.colorbar(sc)
plt.show(block=False)


plt.figure()
# plt.scatter(x, y, c=np.log10(mach), vmin=-1, vmax=1, cmap="coolwarm")
plt.scatter(x, y, c=mach, cmap="coolwarm")
plt.plot(X_BS, Y_BS, label="BS")
plt.plot(X1, Y1, label="MPB", c="C2")
plt.xlabel("X (RM)")
plt.ylabel("Y (RM)")
plt.title("Z=0, Pdyn/Pth")
plt.xlim([0, 2])
plt.ylim([-1, 1])
plt.colorbar()
plt.show(block=False)
