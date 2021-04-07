import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

from funciones import donde

path = "../../../datos/simulacion_leonardo/"
# variables = np.loadtxt(path + "variables.txt")
# np.save(path+'variables.npy', variables)

"""it year mo dy hr mn sc msc X Y (0 a 9)
Z Rho Ux Uy Uz Bx By Bz P HpRho (10 a 19)
HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP OpRho (20 a 29)
OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x (30 a 39)
b1y b1z e jx jy jz (40 a 45)
"""


"""
Vamos a escribir cada término de la ley de Ohm generalizada.
Necesito los valores de v, B, j y p en una región más amplia que el eje x nomás
"""

datos_enteros = np.loadtxt(path + "simu_tray.sat", skiprows=8641)

# Lo "diezmamos" para que no repita valores de x
# x_full = reordenados_full[:, 0]
# idx = [i for i in range(len(x_full) - 1) if x_full[i] != x_full[i + 1]]

pi = 0
pf = 600
datos = datos_enteros[pi:pf]

x = datos[:, 8]
y = datos[:, 9]
z = datos[:, 10]
B = datos[:, 15:18]

velocidad = datos[:, 12:15]
v_i = datos[:, 20:23]
J = datos[:, -3:]

rho = datos[:, 11]
HpRho = datos[:, 19]
O2pRho = datos[:, 24]
OpRho = datos[:, 29]
CO2pRho = datos[:, 34]

Ptermica = datos[:, 18]
presion_H = datos[:, 23]
presion_O2 = datos[:, 28]
presion_O = datos[:, 33]
presion_CO2 = datos[:, 38]

x2 = 1.3234990984508788  # sacados de trayectoria.py
x3 = 1.2017981543969254

inicio_MPB = donde(x, x2)
fin_MPB = donde(x, x3)

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
P_SI = Ptermica * 1e-9  # Pa

Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)
Epx = (
    -1 / (e_SI * n_SI) * np.gradient(P_SI, x * 3390e3)
)  # ojo que esto es una derivada direccional más que un gradiente

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
ax1.axvspan(x2, x3, color="red", alpha=0.2)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("E_cv (mV/m)")
ax1.set_xlim([1, 1.5])

ax2.plot(x, B)
ax2.axvspan(x2, x3, color="red", alpha=0.2)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_xlim([1, 1.5])

ax3.plot(x, v_i)
ax3.axvspan(x2, x3, color="red", alpha=0.2)
ax3.legend(["x", "y", "z", "MPB"], loc="lower left")
ax3.set_ylabel("v_i (km/s)")
ax3.set_xlim([1, 1.5])
ax3.set_ylim(ymin=-200)
ax3.set_xlabel("x (RM)")

plt.show()

plt.figure()
ax1 = plt.subplot2grid((3, 1), (0, 0))
ax2 = plt.subplot2grid((3, 1), (1, 0))
ax3 = plt.subplot2grid((3, 1), (2, 0))
# plt.plot(x[inicio_MPB:fin_BS], np.linalg.norm(Ecv[inicio_MPB:fin_BS, :], axis=1) * 1e3)

ax1.plot(x, Ehall * 1e3)
ax1.axvspan(x2, x3, color="red", alpha=0.2)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("E_hall (mV/m)")
ax1.set_xlim([1, 1.5])
ax1.set_ylim([-1, 5])

ax2.plot(x, B)
ax2.axvspan(x2, x3, color="red", alpha=0.2)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_xlim([1, 1.5])
# ax2.set_ylim(ymax=51)

ax3.plot(x, J * 1e3)
ax3.axvspan(x2, x3, color="red", alpha=0.2)
ax3.legend(["x", "y", "z", "MPB"])
ax3.set_ylabel("J (uA/m²)")
ax3.set_xlim([1, 1.5])
ax3.set_ylim(ymin=-100)
ax3.set_xlabel("x (RM)")

plt.show()


plt.figure()
ax1 = plt.subplot2grid((1, 1), (0, 0))

ax1.plot(x, Epx * 1e3, label="x")
ax1.set_ylabel("E_p (mV/m)")
ax1.set_ylim([-2.5, 10])
ax1.legend()
ax1.axvspan(x2, x3, color="red", alpha=0.2, label="MPB")
ax1.set_xlim([1, 1.5])
ax1.set_xlabel("x (R_M)")

plt.figure()
ax2 = plt.subplot2grid((3, 1), (0, 0))
ax3 = plt.subplot2grid((3, 1), (1, 0))
ax4 = plt.subplot2grid((3, 1), (2, 0))
ax2.plot(x, P_SI * 1e9)
ax2.set_ylabel("Presión (nPa)")
ax2.axvspan(x2, x3, color="red", alpha=0.2, label="MPB")
ax2.set_xlim([1, 1.5])
ax2.set_ylim([0, 0.8])
ax2.set_xlabel("x (R_M)")

ax3.plot(
    y[inicio_MPB - 10 : fin_MPB + 10],
    P_SI[inicio_MPB - 10 : fin_MPB + 10] * 1e9,
    c="C1",
)
ax3.set_ylabel("Presión (nPa)")
ax3.set_xlabel("y (RM)")
ax3.axvspan(y[inicio_MPB], y[fin_MPB], color="red", alpha=0.2, label="MPB")
ax3.set_xlim([-0.2, 0.3])
ax3.set_ylim([0, 0.8])

ax4.plot(
    z[inicio_MPB - 10 : fin_MPB + 10],
    P_SI[inicio_MPB - 10 : fin_MPB + 10] * 1e9,
    c="C2",
)
ax4.set_ylabel("Presión (nPa)")
ax4.set_xlabel("z (RM)")
ax4.set_xlim([0, 0.5])
ax4.set_ylim([0, 0.8])
ax4.axvspan(z[inicio_MPB], z[fin_MPB], color="red", alpha=0.2, label="MPB")
plt.show()

plt.figure()
plt.plot(x, np.linalg.norm(Ecv, axis=1) * 1e3)
plt.plot(x, np.linalg.norm(Ehall, axis=1) * 1e3)
plt.plot(x, Epx * 1e3)
plt.axvspan(x2, x3, color="red", alpha=0.2)
plt.legend(["|E cv|", "|E hall|", "Epx", "MPB"])
plt.xlim([1, 1.5])
plt.ylim([0, 13])
plt.title("Términos del campo eléctrico sobre la órbita")
plt.xlabel("x (R_M)")
plt.ylabel("E (mV/m)")
plt.show()
