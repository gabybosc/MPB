import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# from matplotlib.widgets import MultiCursor
from importar_datos import importar_mag, importar_swia
import sys

sys.path.append("..")

from funciones_plot import equal_axes  # , onpick1
from funciones import diezmar, donde

path = "../../../datos/simulacion_leonardo/"
datos_enteros = np.loadtxt(path + "simu_tray.sat", skiprows=8641)

"""
Mismos análisis pero para los datos sobre la trayectoria de la nave
"""

"""it year mo dy hr mn sc msc X Y (0 a 9)
Z Rho Ux Uy Uz Bx By Bz P HpRho (10 a 19)
HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP OpRho (20 a 29)
OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x (30 a 39)
b1y b1z e jx jy jz (40 a 45)
"""

# Datos de la simulación


pi = 0
pf = 600
mu0 = 4 * np.pi * 1e-7
datos = datos_enteros[pi:pf]

x = datos[:, 8]
y = datos[:, 9]
z = datos[:, 10]
B = datos[:, 15:18]

velocidad_plasma = datos[:, 12:15]
velocidad_i = datos[:, 20:23]
J = datos[:, -3:]

rho = datos[:, 11]
HpRho = datos[:, 19]
O2pRho = datos[:, 24]
OpRho = datos[:, 29]
CO2pRho = datos[:, 34]

P_plasma = datos[:, 18]
presion_H = datos[:, 23]
presion_O2 = datos[:, 28]
presion_O = datos[:, 33]
presion_CO2 = datos[:, 38]


kT_H = presion_H / HpRho
kT_O = presion_O / OpRho
kT_O2 = presion_O2 / O2pRho
kT_CO2 = presion_CO2 / CO2pRho

P_heavy = presion_O + presion_O2 + presion_CO2
P_B = np.linalg.norm(B, axis=1) ** 2 * 1e-9 / (2 * mu0)
Pe = P_plasma - P_heavy - presion_H
P_ram = 1.67e-6 * HpRho * velocidad_i[:, 0] ** 2  # nPa
P_total = P_heavy + P_B + Pe + P_ram + presion_H

beta = P_plasma / P_B

density_mean = [np.mean(HpRho[i : i + 50]) for i in range(len(HpRho) - 50)]
ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km
ion_length = np.append(ion_length, ion_length[-51:-1])
plt.plot(x, kT_H)
plt.plot(x, Pe / rho)
# plt.plot(x, kT_O)
# plt.plot(x, kT_O2)
# plt.plot(x, kT_CO2)
plt.show()


# Datos de MAVEN
mag, t, B_mag, posicion = importar_mag(2016, "03", 16, 17.7, 19)
swia, t_swia, proton_density = importar_swia(2016, "03", 16, 17.7, 19)

# Datos del análisis de MAVEN
R = [1.082, -0.064, 0.515]
normal = [0.920, -0.302, 0.251]
j_maven = 23.2  # mA/m

# hay 10 órbitas, cada una de 1620 puntos
"""Plotea nuestra MPB y la órbita de la simu que estamos usando"""
L = 0.96
x0 = 0.78
e = 0.9
theta = np.linspace(0, np.pi * 3 / 4, 100)
phi = np.linspace(0, 2 * np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(THETA))
# r1 = L / (1 + e * np.cos(THETA))

X1 = x0 + r1 * np.cos(THETA)
Y1 = r1 * np.sin(THETA) * np.cos(PHI)
Z1 = r1 * np.sin(THETA) * np.sin(PHI)

# ahora plotea
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection="3d")
ax.set_xlabel(r"$X_{MSO} (R_m)$", fontsize=14)
ax.set_ylabel(r"$Y_{MSO} (R_m)$", fontsize=14)
ax.set_zlabel(r"$Z_{MSO} (R_m)$", fontsize=14)
ax.plot_surface(
    X1,
    Y1,
    Z1,
    rstride=4,
    cstride=4,
    alpha=0.5,
    # edgecolor="gray",
    cmap=plt.get_cmap("Blues_r"),
)
ax.scatter(R[0], R[1], R[2], label="MAVEN", color="k", s=40)

u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
ax.plot_wireframe(
    np.cos(u) * np.sin(v),
    np.sin(u) * np.sin(v),
    np.cos(v),
    color="#c1440e",
    linewidth=0.5,
)
# sc = ax.scatter(x, y, z, vmin=80, vmax=150, c=ion_length, cmap="magma")
sc = ax.scatter(x, y, z, c=np.log10(beta), vmin=-1.5, vmax=1.5, cmap="PiYG")
plt.colorbar(sc)
equal_axes(ax, X1, Y1, Z1)
ax.set_title("log10(beta)")
plt.show()


# esta sección comentada es para ver el B en las distintas órbitas.

# plt.figure()
# for i in range(6):
#     minimo = int(i * 1620)
#     maximo = int(minimo + 600)
#     plt.plot(x[minimo:maximo], B[minimo:maximo, 0], label=i)
# for i in range(6):
#     m = int(10200 + 1620 * i)
#     M = int(m + 600)
#     plt.plot(x[m:M], B[m:M, 0], label=i)
#
# plt.show()

"""
Comparación de B y de la densidad de protones
"""

t1t2t3t4 = [18.2167, 18.2204, 18.235, 18.2476]
i1 = donde(t, t1t2t3t4[0])
i2 = donde(t, t1t2t3t4[1])
i3 = donde(t, t1t2t3t4[2])
i4 = donde(t, t1t2t3t4[3])

idx = diezmar(t, t_swia)

plt.figure()
plt.plot(posicion[:, 0] / 3390, np.linalg.norm(B_mag, axis=1), label="MAVEN")
plt.axvline(posicion[i1, 0] / 3390, color="C2", linestyle="--", label="t1")
plt.axvline(posicion[i2, 0] / 3390, color="C3", linestyle="--", label="t2")
plt.axvline(posicion[i3, 0] / 3390, color="C4", linestyle="--", label="t3")
plt.axvline(posicion[i4, 0] / 3390, color="C5", linestyle="--", label="t4")
plt.xlabel("x (RM)")
plt.ylabel("|B| (nT)")
plt.legend()

plt.figure()
plt.plot(
    posicion[:, 0] / 3390 + 600 / 3390, np.linalg.norm(B_mag, axis=1), label="MAVEN"
)
plt.plot(x, np.linalg.norm(B, axis=1), label="simulación")
plt.axvline(x=R[0] + 600 / 3390, color="k", linestyle="--", label="cruce MAVEN")
plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
plt.xlabel("x (RM)")
plt.ylabel("|B| (nT)")
plt.legend()

plt.figure()
plt.plot(posicion[idx, 0] / 3390, proton_density, label="swia")
plt.plot(x, HpRho, label="simulación")
plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
plt.ylim(ymax=30)
plt.xlabel("x (RM)")
plt.ylabel("Densidad de protones")
plt.legend()

plt.show()

# Descomentar esto si quiero seleccionar la MPB de la simu

# happy = False
# while not happy:
#     val = []
#     while len(val) < 4:
#         plt.clf()  # clear figure
#         fig = plt.figure(1, constrained_layout=True)
#         fig.subplots_adjust(
#             top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
#         )
#         plt.title("Spacebar when ready to click:")
#
#         ax1 = plt.subplot2grid((2, 1), (0, 0))
#         plt.plot(x, np.linalg.norm(B, axis=1))
#         ax1.set_ylabel("|B| (nT)")
#
#         ax5 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)
#         ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
#         ax5.set_xlabel("x (RM)")
#         plt.plot(x, HpRho)
#
#         fig.canvas.mpl_connect("pick_event", onpick1)
#         multi = MultiCursor(fig.canvas, (ax1, ax5), color="black", lw=1)
#
#         zoom_ok = False
#         print("\nSpacebar when ready to click:\n")
#         while not zoom_ok:
#             zoom_ok = plt.waitforbuttonpress(-1)
#         print("Click to select MPB: ")
#         val = np.asarray(plt.ginput(4))[:, 0]
#         print("Selected values: ", val)
#         outs = sorted(val)
#
#     print("Happy? Keyboard click for yes, mouse click for no.")
#     happy = plt.waitforbuttonpress()
#
# plt.show()

"""
Cálculo de J = n x (Bu-Bd)
"""
# las posiciones de entrada y salida de la MPB en el eje x en RM
x1 = 1.393042495053138
x2 = 1.3234990984508788
x3 = 1.2017981543969254
x4 = 1.1474673758014104

x23 = (x2 - x3) * 3390e3  # ancho en m
ancho_updown = 0.015 * 13000 / 3390

inicio_up = donde(x, x1 + ancho_updown)
fin_up = donde(x, x1)
inicio_down = donde(x, x4)
fin_down = donde(x, x4 - ancho_updown)
inicio_MPB = donde(x, x2)
fin_MPB = donde(x, x3)

n2 = [0.856, -0.066, 0.512]
B_upstream = np.mean(B[inicio_up:fin_up], axis=0)
B_downstream = np.mean(B[inicio_down:fin_down], axis=0)
J_s = np.cross(n2, (B_upstream - B_downstream)) / mu0  # nA/m

J_integrado = np.trapz(J[inicio_MPB:fin_MPB, 0], x[inicio_MPB:fin_MPB])
# Integramos el J de la simulación en x para obtener un Js

# plt.figure()
# plt.plot(x, np.linalg.norm(J_integrado, axis=1) * 1e-3)
# plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
# plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
# plt.axhline(y=np.linalg.norm(J_s * 1e-6), color="C1", label="J = n x (Bu-Bd)")
# plt.axhline(y=np.linalg.norm(j_maven), color="C2", label="J = n x (Bu-Bd)")
# plt.xlim(xmin=1)
# plt.xlabel("x MSO (RM)")
# plt.ylabel("J_s  (mA/m)")
# plt.legend()
# plt.show()


plt.figure()
plt.plot(x, np.linalg.norm(J, axis=1) * 1e3)
plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
plt.axhline(y=np.linalg.norm(J_s / x23), color="C1", label="J = n x (Bu-Bd) / h (simu)")
plt.ylim(ymax=100)
plt.xlabel("x (RM)")
plt.ylabel("j_v (nA/m²)")
plt.legend()
plt.show()

plt.plot(x, OpRho, label="O+")
plt.plot(x, O2pRho, label="O2+")
plt.plot(x, CO2pRho, label="CO2+")
plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
plt.legend()
plt.ylim(-1, 50)
plt.xlim(xmin=0.9)
plt.xlabel("x (RM)")
plt.ylabel("rho (mp/cc)")

"""
Análsis de las presiones
"""

plt.figure()
sc = plt.scatter(x, y, c=np.log10(beta), vmin=-4, vmax=4, s=35, cmap="coolwarm",)
plt.xlabel("X (RM)")
plt.ylabel("Y (RM)")
plt.title("log10(beta)")
plt.colorbar(sc)


plt.figure()
plt.scatter(x, P_total, label="P total")
plt.scatter(x, Pe, label="P electrónica")
plt.scatter(x, P_B, label="P magnética")
plt.scatter(x, P_ram, label="P ram")
plt.scatter(x, presion_H, label="P protones", marker=".")
plt.scatter(x, P_heavy, label="P heavies")
plt.axvspan(x2, x3, color="red", alpha=0.2, label="MPB")
plt.xlim([1, 1.5])
plt.ylim([0, 1])
plt.legend()
plt.xlabel("x (RM)")
plt.ylabel("presion (nPa)")

plt.show()
