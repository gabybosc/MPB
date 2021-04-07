import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys


path = "../../../datos/simulacion_leonardo/"
datos = np.load(
    path + "ordenado_cut_0dot05.npy"
)  # todos los datos ordenados con y,z menores a 0.05

sys.path.append("..")
from funciones import donde


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

#
# inicio_MPB = donde(pos[:, 0], 1.2)
# fin_MPB = donde(pos[:, 0], 1.36)
# inicio_BS = donde(pos[:, 0], 1.67)
# fin_BS = donde(pos[:, 0], 1.72)
#


# def sigmoid(x, L, x0, k, b):
#     y = L / (1 + np.exp(-k * (x - x0))) + b
#     return y


i = donde(x, 1.1)
f = donde(x, 1.3)

plt.figure()
plt.plot(x[i:f], HpRho[i:f])
plt.plot(x[i:f], O2pRho[i:f])
plt.plot(x[i:f], OpRho[i:f])
plt.plot(x[i:f], CO2pRho[i:f])
plt.show()


xx = x[i:f] * 3390e3


def lineal(x, a, b):
    y = a * x + b
    return y


grav = []

for densidad in [OpRho[i:f], O2pRho[i:f], CO2pRho[i:f]]:
    logfit = np.log(densidad * 1e6)

    # coef = np.polyfit(xx, logfit, 1)
    # poly1d_fn = np.poly1d(coef)
    #
    # plt.plot(xx, logfit, "yo", xx, poly1d_fn(xx), "--k")

    guess = [-50, 63]
    popt, pcov = curve_fit(lineal, xx, logfit, guess, method="lm")
    yy = lineal(xx, *popt)

    plt.figure()
    plt.plot(xx, logfit, "o", label="datos")
    plt.plot(xx, yy, "r-", label="Fit")
    plt.legend()

    plt.figure()
    plt.plot(xx, densidad * 1e6)
    plt.plot(xx, np.exp(yy))

    plt.show()

    a = popt[0]  # tiene que tener unidades de 1/m
    b = popt[1]  # adimensional, pero e^b tiene unidades 1/m³
    int = np.exp(popt[0] * xx + popt[1]) / (popt[0])  # esto está en 1/m²
    grav.append(-int * mp * g * 1e9)  # en nPan   # es menos la integral
    print(
        "deberian ser iguales",
        sum(-int * mp * g * 1e9),
        np.trapz(densidad * 1e6 * mp * g, xx) * 1e9,
    )

P_grav = np.zeros(len(P_heavy))
P_grav[i:f] = grav[0] + grav[1] + grav[2]

P_total = P_heavy + P_B + Pe + P_ram + HP + P_grav


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
