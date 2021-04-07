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

reordenados_full = np.load(path + "ordenado_cut_0dot05.npy")

# Lo "diezmamos" para que no repita valores de x
x_full = reordenados_full[:, 0]
idx = [i for i in range(len(x_full) - 1) if x_full[i] != x_full[i + 1]]

reordenados = reordenados_full[idx]

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

# idx = [i for i in range(len(x) - 1) if x[i] != x[i + 1]]
# P_diezmado = P_SI[idx]
# x_diezmado = x[idx] * 3390e3
# n_diezmado = n_SI[idx]
# Ep = -1 / (e_SI * n_diezmado) * np.gradient(P_diezmado, x_diezmado)


Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)
# Ep = -1 / (e_SI * n_SI) * np.gradient(P_SI, (x[10] - x[7]) * 3390e3)
Ep = -1 / (e_SI * n_SI) * np.gradient(P_SI, x * 3390e3)

J_mean = np.mean(J[inicio_MPB:fin_MPB], axis=0)
B_mean = np.mean(B[inicio_MPB:fin_MPB], axis=0)

print(
    f"corriente simu j = {J_mean*1e3} nA/m², corriente MAVEN j = (-34.773,-233.373,-153.737) nA/m²"
)
print(f"campo simu B = {B_mean} nT, campo MAVEN B = (11.71,26.29,-19.55) nT")

plt.figure()
ax1 = plt.subplot2grid((3, 1), (0, 0))
ax2 = plt.subplot2grid((3, 1), (1, 0))
ax3 = plt.subplot2grid((3, 1), (2, 0))
# plt.plot(x[inicio_MPB:fin_BS], np.linalg.norm(Ecv[inicio_MPB:fin_BS, :], axis=1) * 1e3)

ax1.plot(x, Ecv * 1e3)
ax1.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
ax1.axvline(x=1.708, c="black", ls="--")
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("E_cv (mV/m)")
ax1.set_xlim([1.1, 1.9])

ax2.plot(x, B)
ax2.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
ax2.axvline(x=1.708, c="black", ls="--")
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_xlim([1.1, 1.9])

ax3.plot(x, v_i)
ax3.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
ax3.axvline(x=1.708, c="black", ls="--")
ax3.legend(["x", "y", "z", "BS", "MPB"])
ax3.set_ylabel("v_i (km/s)")
ax3.set_xlim([1.1, 1.9])
ax3.set_xlabel("x (RM)")

plt.show()

plt.figure()
ax1 = plt.subplot2grid((3, 1), (0, 0))
ax2 = plt.subplot2grid((3, 1), (1, 0))
ax3 = plt.subplot2grid((3, 1), (2, 0))
# plt.plot(x[inicio_MPB:fin_BS], np.linalg.norm(Ecv[inicio_MPB:fin_BS, :], axis=1) * 1e3)

ax1.plot(x, Ehall * 1e3)
ax1.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
ax1.axvline(x=1.708, c="black", ls="--")
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("E_hall (mV/m)")
ax1.set_xlim([1.1, 1.9])
ax1.set_ylim([-1, 5])

ax2.plot(x, B)
ax2.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
ax2.axvline(x=1.708, c="black", ls="--")
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_xlim([1.1, 1.9])
ax2.set_ylim(ymax=51)

ax3.plot(x, J * 1e3)
ax3.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
ax3.axvline(x=1.708, c="black", ls="--")
ax3.legend(["x", "y", "z", "BS", "MPB"])
ax3.set_ylabel("J (uA/m²)")
ax3.set_xlim([1.1, 1.9])
ax3.set_ylim([-0.07, 0.07])
ax3.set_xlabel("x (RM)")

plt.show()


plt.figure()
plt.plot(x[inicio_MPB:fin_BS], Ep[inicio_MPB:fin_BS] * 1e3)
plt.xlabel("x (R_M)")
plt.ylabel("E_p (mV/m)")

plt.figure()
plt.plot(x, Ecv[:, 0] * 1e3)
plt.plot(x, Ehall[:, 0] * 1e3)
plt.plot(x, Ep * 1e3)
plt.axvspan(x[inicio_MPB], x[fin_MPB], color="red", alpha=0.2)
plt.axvline(x=1.708, c="black", ls="--")
plt.legend(["E cv", "E hall", "Ep", "BS", "MPB"])
plt.xlim([1.1, 1.9])
plt.ylim([-5, 5])
plt.title("Términos del campo eléctrico sobre el eje x")
plt.xlabel("x (R_M)")
plt.ylabel("E (mV/m)")
plt.show()
