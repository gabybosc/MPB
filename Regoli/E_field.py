import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

from funciones import donde

path = "../../../datos/simulacion_leonardo/"
# variables = np.loadtxt(path + "variables.txt")
# np.save(path+'variables.npy', variables)

"""
Rho Ux Uy Uz Bx By Bz P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x b1y b1z e jx jy jz
"""

"""
Vamos a escribir cada término de la ley de Ohm generalizada.
Necesito los valores de v, B, j y p en una región más amplia que el eje x nomás
"""

reordenados = np.load(path + "ordenado_cut_0dot05.npy")

x = reordenados[:, 0]
pos = reordenados[:, :3]
rho = reordenados[:, 3]
v_e = reordenados[:, 4:7]
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
J = reordenados[:, -3:]

inicio_MPB = donde(pos[:, 0], 1.2)
fin_MPB = donde(pos[:, 0], 1.36)
inicio_BS = donde(pos[:, 0], 1.67)
fin_BS = donde(pos[:, 0], 1.72)

#
# xi = 1.2
# xf = 1.36
# inicio_up = donde(x, xi - 126)
# fin_up = donde(x, xi)
# B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT
#
# inicio_down = donde(x, xf)
# fin_down = donde(x, xf + 146)
# B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

e_SI = 1.6e-19  # C

v_SI = v_i * 1e3  # m/s
B_SI = B * 1e-9  # T
n_SI = HpRho * 1e6  # 1/m3
J_SI = J * 1e-6  # A/m2
P_SI = Ptot * 1e-9  # Pa

Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)
Ep = -1 / (3390 * e_SI * n_SI) * np.gradient(P_SI, x[10] - x[7] * 3390)

J_mean = np.mean(J[inicio_MPB:fin_MPB], axis=0)
B_mean = np.mean(B[inicio_MPB:fin_MPB], axis=0)

print(
    f"corriente simu j = {J_mean*1e3} nA/m², corriente MAVEN j = (-34.773,-233.373,-153.737) nA/m²"
)
print(f"campo simu B = {B_mean} nT, campo MAVEN B = (11.71,26.29,-19.55) nT")

plt.figure()
plt.plot(x[inicio_MPB:fin_BS], Ecv[inicio_MPB:fin_BS, :] * 1e3)
plt.plot(x[inicio_MPB:fin_BS], np.linalg.norm(Ecv[inicio_MPB:fin_BS, :], axis=1) * 1e3)
plt.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
plt.legend(["x", "y", "z", "|B|", "MPB"])
plt.xlabel("x (R_M)")
plt.ylabel("E_cv (mV/m)")

plt.figure()
plt.plot(x[inicio_MPB:fin_BS], Ehall[inicio_MPB:fin_BS, :] * 1e3)
plt.legend(["x", "y", "z"])
plt.xlabel("x (R_M)")
plt.ylabel("E_hall (mV/m)")

plt.figure()
plt.plot(x[inicio_MPB:fin_BS], Ep[inicio_MPB:fin_BS])
plt.xlabel("x (R_M)")
plt.ylabel("E_p (mV/m)")

plt.figure()
plt.plot(x[inicio_MPB:fin_BS], Ecv[inicio_MPB:fin_BS, 0] * 1e3)
plt.plot(x[inicio_MPB:fin_BS], Ehall[inicio_MPB:fin_BS, 0] * 1e3)
plt.plot(x[inicio_MPB:fin_BS], Ep[inicio_MPB:fin_BS] * 1e3)
plt.legend(["cv", "hall", "p"])
plt.xlabel("x (R_M)")
plt.ylabel("E (mV/m)")
plt.show()

plt.show()
