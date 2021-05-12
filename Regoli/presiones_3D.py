import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


path = "../../../datos/simulacion_leonardo/"
datos = np.load(path + "cubo_Daniel.npy")  # 1 < x < 2 , -1.5 < y,z < 1.5

sys.path.append("..")
from funciones import donde


"""
Rho Ux Uy Uz Bx By Bz P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x b1y b1z e jx jy jz
"""

"""
Plots de las diferentes presiones en 3D
"""

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
g = 3.7  # Mars surface gravity, m/s²

# reordenados = []
# for i in range(len(datos) - 1):
#     if datos[i, 0] != datos[i + 1, 0]:
#         reordenados.append(datos[i, :])
#
# reordenados = np.array(reordenados)

reordenados = datos

x = reordenados[:, 0]
pos = reordenados[:, :3]
rho = reordenados[:, 3]
B = reordenados[:, 7:10]
J = reordenados[:, -3:]
velocidad = reordenados[:, 4:7]  # km/s
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

beta1 = np.array([beta[i] for i in range(len(beta)) if 0.9 < beta[i] < 1.1])
superficie_beta = np.array([pos[i, :] for i in range(len(beta)) if 0.9 < beta[i] < 1.1])

idx = []
for i in range(len(P_B)):
    if P_B[i] > 0.5:
        idx.append(i)

fig = plt.figure()
ax = fig.gca(projection="3d")
p = ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=Ptermica, marker=".", cmap="hot", s=1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.colorbar(p, ax=ax)
plt.title("P termica")

fig = plt.figure()
ax = fig.gca(projection="3d")
p = ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=P_ram, marker=".", cmap="hot", s=1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.colorbar(p, ax=ax)
plt.title("P ram")


fig = plt.figure()
ax = fig.gca(projection="3d")
p = ax.scatter(
    pos[idx, 0],
    pos[idx, 1],
    pos[idx, 2],
    c=P_B[idx],
    marker=".",
    cmap="YlOrRd",
    s=1,
    vmax=5,
)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.colorbar(p, ax=ax)
plt.title("P magnética")

fig = plt.figure()
ax = fig.gca(projection="3d")
ax.scatter(
    superficie_beta[:, 0], superficie_beta[:, 1], superficie_beta[:, 2],
)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("sup beta ~ 1")
plt.show()
