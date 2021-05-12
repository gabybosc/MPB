import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

sys.path.append("..")

from funciones import donde, angulo

path = "../../../datos/simulacion_leonardo/"
posicion = np.load(path + "pos_mhd.npy")
variables = np.load(path + "variables.npy")

"""
Rho Ux Uy Uz Bx By Bz P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x b1y b1z e jx jy jz
"""

ceros = [i for i in range(len(posicion)) if posicion[i, 1] == 0]

pos_cut = np.delete(posicion, ceros, axis=0)
var_cut = np.delete(variables, ceros, axis=0)

val = np.concatenate((pos_cut, var_cut), axis=1)

X = np.array([val[i, :] for i in range(len(val)) if 0 < val[i, 0] < 3])
Y = np.array([X[i, :] for i in range(len(X)) if np.abs(X[i, 1]) <= 3])
Z = np.array([Y[i, :] for i in range(len(Y)) if np.abs(Y[i, 2]) <= 3])
# Z es el cut final y el que voy a usar el resto del tiempo

# quiero reordenar los datos de forma creciente en x
reordenados = np.array(sorted(Z, key=lambda f: f[0]))

pos = reordenados[:, :3]

inicio_MPB = donde(pos[:, 0], 1.2)
fin_MPB = donde(pos[:, 0], 1.36)
inicio_BS = donde(pos[:, 0], 1.67)
fin_BS = donde(pos[:, 0], 1.72)

MPB = donde(pos[:, 0], 1.26)  # 1.26 es donde es máx la J en x.
BS = donde(pos[:, 0], 1.7)

B = reordenados[:, 7:10]
J = reordenados[:, -3:] * 1000

print("Angulo J, B en la MPB:", angulo(J[MPB, :], B[MPB, :]) * 180 / np.pi)
print("Angulo J, B en el BS:", angulo(J[BS, :], B[BS, :]) * 180 / np.pi)


MPB = np.array([Z[i, :] for i in range(len(Z)) if 1.2 <= Z[i, 0] <= 1.36])
pos = MPB[:, :3]

B = MPB[:, 7:10]
J = MPB[:, -3:] * 1000

plt.plot(J[:, 1], J[:, 2], ".")
plt.plot(B[:, 1], B[:, 2], ".")
plt.show()

#
# fig = plt.figure()
# ax = fig.gca(projection="3d")
# ax.quiver(
#     0, 0, 0, B[MPB, 0], B[MPB, 1], B[MPB, 2], length=0.1, normalize=True, label="B"
# )
# ax.quiver(
#     0,
#     0,
#     0,
#     J[MPB, 0],
#     J[MPB, 1],
#     J[MPB, 2],
#     color="C1",
#     length=0.1,
#     normalize=True,
#     label="J",
# )
# ax.quiver(
#     0,
#     0,
#     0,
#     JxB[MPB, 0],
#     JxB[MPB, 1],
#     JxB[MPB, 2],
#     color="C2",
#     length=0.1,
#     normalize=True,
#     label="JxB",
# )
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# plt.legend()
# plt.title("MPB")
#
#
# fig = plt.figure()
# ax = fig.gca(projection="3d")
# ax.quiver(0, 0, 0, B[BS, 0], B[BS, 1], B[BS, 2], length=0.1, normalize=True, label="B")
# ax.quiver(
#     0,
#     0,
#     0,
#     J[BS, 0],
#     J[BS, 1],
#     J[BS, 2],
#     color="C1",
#     length=0.1,
#     normalize=True,
#     label="J",
# )
# ax.quiver(
#     0,
#     0,
#     0,
#     JxB[BS, 0],
#     JxB[BS, 1],
#     JxB[BS, 2],
#     color="C2",
#     length=0.1,
#     normalize=True,
#     label="JxB",
# )
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# plt.legend()
# plt.title("BS")

# plt.show()
