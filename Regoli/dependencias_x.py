import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

from funciones import donde, Mij
from funciones_plot import hodograma


"""
Rho Ux Uy Uz Bx By Bz P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x b1y b1z e jx jy jz
"""

"""
Quiero plotear el campo y ver que me aparezca una MPB y un BS.
Primero, voy a buscar las posiciones y=0, z=0 y con eso plotear para esos índices
|B| vs x.
Veo que efectivamente donde tengo la MPB aparece un salto en J.
"""

path = "../../../datos/simulacion_leonardo/"
reordenados = np.load(
    path + "ordenado_cut_0dot05.npy"
)  # todos los datos ordenados con y,z menores a 0.05
limite = 0.05

x = reordenados[:, 0]
pos = reordenados[:, :3]
rho = reordenados[:, 3]
v_plasma = reordenados[:, 4:7]
B = reordenados[:, 7:10]
Ptot = reordenados[:, 10]
HpRho = reordenados[:, 11]
v_i = reordenados[:, 12:15]
HP = reordenados[:, 15]
O2pRho = reordenados[:, 16]
O2P = reordenados[:, 20]
OpRho = reordenados[:, 21]
OP = reordenados[:, 25]
CO2pRho = reordenados[:, 26]
CO2P = reordenados[:, 30]
J = reordenados[:, -3:] * 1000

inicio_MPB = donde(pos[:, 0], 1.2)
fin_MPB = donde(pos[:, 0], 1.36)
inicio_BS = donde(pos[:, 0], 1.67)
fin_BS = donde(pos[:, 0], 1.72)

B_MPB = reordenados[inicio_MPB:fin_MPB, 7:10]
J_MPB = reordenados[inicio_MPB:fin_MPB, -3:] * 1000

# MVA:

M_ij = Mij(B_MPB)
[lamb, av] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica
idx = lamb.argsort()[::-1]
lamb = lamb[idx]
av = av[:, idx]
x1 = av[:, 0]
x2 = av[:, 1]
x3 = av[:, 2]

av = np.concatenate([x1, x2, x3])

if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
    x3 = -x3
if any(np.cross(x1, x2) - x3) > 0.01:
    print("Cambio el signo de x1 para que los av formen terna derecha")
    x1 = -x1

# las proyecciones
B1 = np.dot(B_MPB, x1)
B2 = np.dot(B_MPB, x2)
B3 = np.dot(B_MPB, x3)

hodograma(B1, B2, B3)

xi = 1.2
xf = 1.36
inicio_up = donde(x, xi - 126)
fin_up = donde(x, xi)
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(x, xf)
fin_down = donde(x, xf + 146)
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT


omega = (
    np.arccos(
        np.dot(B_upstream, B_downstream)
        / (np.linalg.norm(B_upstream) * np.linalg.norm(B_downstream))
    )
    * 180
    / np.pi
)

print(f"Bup = {B_upstream} \n Bdown = {B_downstream}\nomega = {omega}")

# nu_in = (rho[inicio_down] * velocidad[inicio_down, 0])
# nu_out = (rho[fin_down] * velocidad[fin_down, 0])
#
# print(f"Conservación del caudal: {nu_in}= {nu_out}")

q_e = 1.6e-19  # C
m_p = 1.6e-27  # kg

# plt.plot(j)
# plt.plot(J)
# plt.legend(["jx", "jy", "jz", "Jx", "Jy", "Jz"])
# plt.show()

#
# plt.plot(J[:, 1], J[:, 2], label="J")
# plt.plot(B[:, 1], B[:, 2], label="B")
# plt.legend()
# plt.show()


plt.figure()
plt.plot(x, HpRho, label="Densidad H+")
plt.plot(x, OpRho, label="Densidad O+")
plt.plot(x, O2pRho, label="Densidad O2+")
plt.plot(x, CO2pRho, label="Densidad CO2+")
plt.axvline(x=1.25, c="black", ls="--", label="MPB")
plt.axvline(x=1.7, c="m", ls="--", label="BS")
plt.title(f"Densidad de partículas para y,z < {limite}")
plt.legend()
plt.xlabel("x (RM)")
plt.ylim((0, 200))

plt.figure()
plt.plot(x, np.linalg.norm(B, axis=1))
plt.plot(x, B)
plt.title(f"Bx,By,Bz (nT) vs x (RM) para y,z < {limite}")
plt.legend(["|B|", "Bx", "By", "Bz"])
plt.xlim((1, 2))
plt.ylim((-10, 100))

plt.figure()
plt.plot(x, J)
plt.plot(x, np.linalg.norm(J, axis=1), label="|J| (nA/m2)")
plt.title(f"Jx,Jy,Jz (nA/m2) vs x (RM) para y,z < {limite}")
plt.legend(["Jx", "Jy", "Jz", "|J|"])
plt.xlim((1, 2))
plt.ylim((-100, 100))

plt.figure()
plt.plot(x, Ptot - CO2P - HP - OP - O2P, label="Presion electronica?")
plt.plot(x, CO2P, label="Presion CO2+")
plt.plot(x, HP, label="Presion H+")
plt.plot(x, OP, label="Presion O+")
plt.plot(x, O2P, label="Presion O2+")
plt.axvline(x=1.25, c="black", ls="--", label="MPB")
plt.axvline(x=1.7, c="m", ls="--", label="BS")
plt.title(f"Presión para y,z < {limite}")
plt.legend()
plt.xlabel("x (RM)")
# plt.ylim((0, 10))

plt.figure()
plt.plot(x, v_i)

plt.show()
