import numpy as np
import matplotlib.pyplot as plt

"""
Quiero comparar J con rotB y con qn(vi-ve)
"""

path = "../../../datos/simulacion_leonardo/"
reordenados = np.load(
    path + "ordenado_cut_0dot05.npy"
)  # todos los datos ordenados con y,z menores a 0.05

x = reordenados[:, 0]  # Rm
pos = reordenados[:, :3]
rho = reordenados[:, 3]  # 1/cm³
v_e = reordenados[:, 4:7]  # km/s
B = reordenados[:, 7:10]  # nT
Ptot = reordenados[:, 10]  # nPa
HpRho = reordenados[:, 11]
v_i = reordenados[:, 12:15]
HP = reordenados[:, 15]
O2pRho = reordenados[:, 16]
O2P = reordenados[:, 20]
OpRho = reordenados[:, 21]
OP = reordenados[:, 25]
CO2pRho = reordenados[:, 26]
CO2P = reordenados[:, 30]
J = reordenados[:, -3:]  # uA/m²

e = 1.6e-19  # C
mu_0 = 4e-7 * np.pi  # T m / A = nT m / nA

j_ampere = (
    1e9 * e * np.array([HpRho[i] * (v_i[i, :] - v_e[i, :]) for i in range(len(v_i))])
)  # A/m²


def dif_finitas(f, dx):
    diff = np.array([f[i + 1] - f[i] for i in range(len(f) - 1)])
    # agregamos el que falta al final
    diff = np.append(diff, f[-1] - f[-2])
    sol = diff / dx
    return sol


dy = np.abs(reordenados[4, 1] - reordenados[2, 1])
dyBz = np.gradient(B[:, 2], dy * 3390) / 3390000  # nT / m

dz = np.abs(reordenados[4, 2] - reordenados[0, 2])
dzBy = np.gradient(B[:, 1], dz * 3390) / 3390000

j_curl = (dyBz - dzBy) / mu_0  # nA / m²

plt.plot(x, J[:, 0] * 1e3, label="J simu")
plt.plot(x, j_ampere[:, 0] * 1e9, label="J ampère")
plt.plot(x, j_curl, label="J rotor")
plt.xlabel("x (RM)")
plt.ylabel("Jx (nA/m²)")
plt.legend()
plt.show()
