import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

path = "../../../datos/simulacion_leonardo/"
# variables = np.loadtxt(path + "variables.txt")
# np.save(path+'variables.npy', variables)
# posicion = np.load(path + "pos_mhd.npy")
# variables = np.load(path + "variables.npy")
reordenados = np.load(
    path + "ordenado_cut_0dot05.npy"
)  # todos los datos ordenados con y,z menores a 0.05


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

"""
Por como están hechos los datos, cuando una de las posiciones es cero, todas las
variables se anulan. Voy a borrar esos datos entonces.
"""

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
g = 3.7  # Mars surface gravity, m/s²

x = reordenados[:, 0]
pos = reordenados[:, :3]
rho = reordenados[:, 3]
B = reordenados[:, 7:10]
J = reordenados[:, -3:] * 1000
HpRho = reordenados[:, 11]
velocidad = reordenados[:, 4:7]  # km/s

# Presiones
Ptermica = reordenados[:, 10]  # nPa
HP = reordenados[:, 15]
O2P = reordenados[:, 20]
OP = reordenados[:, 25]
CO2P = reordenados[:, 30]

P_heavy = OP + O2P + CO2P
P_B = np.linalg.norm(B, axis=1) ** 2 * 1e-9 / (2 * mu0)
Pe = Ptermica - CO2P - HP - OP - O2P
P_ram = 1.67e-6 * HpRho * velocidad[:, 0] ** 2  # nPa

dxx = np.abs(x[4] - x[0])

grav = np.trapz(-HpRho * 1e6 * mp * g, x * 3390e3) * 1e9
# grav = np.trapz(-HpRho * 1e6 * mp * g, dx=dxx * 3390e3)

# plt.plot(x, -HpRho * 1e6 * mp * g)
# plt.show()

P_grav = grav * np.ones(len(grav))
P_total = P_heavy + P_B + Pe + P_ram + HP + P_grav
#
# inicio_MPB = donde(pos[:, 0], 1.2)
# fin_MPB = donde(pos[:, 0], 1.36)
# inicio_BS = donde(pos[:, 0], 1.67)
# fin_BS = donde(pos[:, 0], 1.72)
#
# B_MPB = reordenados[inicio_MPB:fin_MPB, 7:10]
# J_MPB = reordenados[inicio_MPB:fin_MPB, -3:] * 1000
#
#
# plt.figure()
# plt.plot(x, Pe, label="Presion electronica?")
# plt.plot(x, CO2P, label="Presion CO2+")
# plt.plot(x, HP, label="Presion H+")
# plt.plot(x, OP, label="Presion O+")
# plt.plot(x, O2P, label="Presion O2+")
# plt.axvline(x=1.25, c="black", ls="--", label="MPB")
# plt.axvline(x=1.7, c="m", ls="--", label="BS")
# plt.title(f"Presión para y,z < {limite}")
# plt.legend()
# plt.xlabel("x (RM)")
#
#
plt.figure()
plt.scatter(x, Pe, label="Presion electronica?")
plt.scatter(x, P_heavy, label="Presion heavy")
plt.scatter(x, P_B, label="Presion B")
plt.scatter(x, P_ram, label="Presion ram")
plt.scatter(x, P_grav, label="Presion grav")
plt.scatter(x, P_total, label="presión total")
plt.scatter(x, HP, label="Presion H+", marker=".")
plt.axvline(x=1.25, c="black", ls="--", label="MPB")
plt.axvline(x=1.7, c="m", ls="--", label="BS")
plt.title(f"Presión para y,z < 0.05")
plt.legend()
plt.xlabel("x (RM)")
plt.ylabel("P (nPa)")
plt.ylim(0, 1)


plt.show()
