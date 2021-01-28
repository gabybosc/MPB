import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

path = "../../../datos/simulacion_leonardo/"
datos = np.load(
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

reordenados = []
for i in range(len(datos) - 1):
    if datos[i, 0] != datos[i + 1, 0]:
        reordenados.append(datos[i, :])

reordenados = np.array(reordenados)

x = reordenados[:, 0]
pos = reordenados[:, :3]
rho = reordenados[:, 3]
B = reordenados[:, 7:10]
J = reordenados[:, -3:]
velocidad = reordenados[:, 4:7]  # km/s

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
P_ram = 1.67e-6 * HpRho * velocidad[:, 0] ** 2  # nPa

grav = []
grav_O2 = []
grav_O = []
grav_CO2 = []
for i in range(len(HpRho) - 2):
    grav.append(np.trapz(HpRho[i : i + 2] * 1e6 * mp * g, x[i : i + 2] * 3390e3))
    grav_O2.append(np.trapz(O2pRho[i : i + 2] * 1e6 * mp * g, x[i : i + 2] * 3390e3))
    grav_O.append(np.trapz(OpRho[i : i + 2] * 1e6 * mp * g, x[i : i + 2] * 3390e3))
    grav_CO2.append(np.trapz(CO2pRho[i : i + 2] * 1e6 * mp * g, x[i : i + 2] * 3390e3))

np.trapz(HpRho[i : i + 2] * 1e6 * mp * g, x[i : i + 2] * 3390e3)


def trapecio(f, x):
    sol = np.zeros(len(x))
    for i in range(len(x) - 1):
        sol[i] = (x[i + 1] - x[i]) * (f[i] + f[i + 1]) / 2
    sol[-1] = sol[-2]
    return sol


trapecio(HpRho * 1e6 * mp * g, x * 3390e3)

# grav.append(grav[-3:])
# grav_O.append(grav_O[-2:])
# grav_O2.append(grav_O2[-2:])
# grav_CO2.append(grav_CO2[-2:])

g_total = (
    np.array(grav) + np.array(grav_O) + np.array(grav_O2) + np.array(grav_CO2)
) * 1e9

g_total = np.append(g_total, g_total[-2:])

plt.plot(x, g_total)
plt.show()

P_total = P_heavy + P_B + Pe + P_ram + HP + g_total
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
# plt.scatter(x, P_grav, label="Presion grav")
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

from scipy.optimize import curve_fit
from funciones import donde


def sigmoid(x, L, x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return y


i = donde(x, 1.18)
f = donde(x, 1.7)
xx = x[i:f]

plt.plot(x[i:f], HpRho[i:f])
plt.plot(x[i:f], O2pRho[i:f])
plt.plot(x[i:f], OpRho[i:f])
plt.plot(x[i:f], CO2pRho[i:f])
plt.show()


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


densidad = O2pRho[i:f]
guess = [max(densidad), 1, min(densidad)]  # this is an mandatory initial guess

popt, pcov = curve_fit(func, xx, densidad, guess, method="dogbox")
yp = func(x, *popt)

plt.figure()
plt.plot(x, O2pRho, "ko", label="Original Noised Data")
plt.plot(x, yp, "r-", label="Fitted Curve")
plt.legend()
plt.show()


r = HpRho[i:f]
p0 = [max(r), np.median(xx), 1, min(r)]  # this is an mandatory initial guess

popt, pcov = curve_fit(sigmoid, xx, r, p0, method="lm")

y = sigmoid(xx, *popt)

L = popt[0]
x0 = popt[1]
k = popt[2]
b = popt[3]

integral = mp * g * (b + L) * xx + mp * g * L / k * np.log(1 + np.exp(-k * (xx - x0)))

plt.plot(xx, integral)
plt.show()

# L = popt[0] * 1e6
# x0 = popt[1] * 3390e3
# k = popt[2]
# b = popt[3] * 1e6
#
# integral = mp * g * (b + L) * xx * 3390e3 + mp * g * L / k * np.log(
#     1 + np.exp(-k * (xx * 3390e3 - x0))
# )
#
# plt.plot(xx * 3390e3, integral * 1e9)
# plt.show()

plt.plot(xx, r)
plt.plot(xx, y)
plt.show()
