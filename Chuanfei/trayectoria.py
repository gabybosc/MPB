import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.widgets import MultiCursor
from importar_datos import importar_mag, importar_swica
import sys

sys.path.append("..")

from funciones_plot import equal_axes, onpick1, hodograma
from funciones import diezmar, donde, Mij, angulo, ancho_mpb

# Simulación
path = "../../../datos/simulacion_chuanfei/"
datos_On = np.loadtxt(path + "sat_trajectory_HallOn_new1.sat", skiprows=2)
datos_Off = np.loadtxt(path + "sat_trajectory_HallOff.sat", skiprows=2)
datos_new = np.loadtxt(path + "sat_trajectory_HallOn_new2.sat", skiprows=2)

# Datos de MAVEN
mag, t, B_mag, posicion = importar_mag(2016, "03", 16, 17.85, 19)
swia, t_swia, proton_density = importar_swica(2016, "03", 16, 17.85, 19)

"""
Mismos análisis pero para los datos sobre la trayectoria de la nave
En este intervalo, x(t) es biyectiva así que se puede invertir para tener B(x)
"""

"""
it year mo dy hr mn sc msc X Y (0 a 9)
Z Rho Ux Uy Uz Bx By Bz Pe P (10 a 19)
HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP (20 a 29)
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP (30 a 39)
b1x b1y b1z e jx jy jz (40 a 46)
"""

# Datos de la simulación

mu0 = 4 * np.pi * 1e-7

datos = datos_On[
    1000:1250
]  # le saco los primeros mil puntos porque son el lado de atrás que no me interesa

x = datos[:, 8]  # la posicion es la misma On y Off
y = datos[:, 9]
z = datos[:, 10]
B = datos[:, 15:18]

B_off = datos_Off[1000:1250, 15:18]

velocidad_plasma = datos[:, 12:15]
velocidad_H = datos[:, 21:24]
J = datos[:, -3:]  # uA/m2

J_off = datos_Off[1000:1250, -3:]

rho = datos[:, 11]
HpRho = datos[:, 20]
O2pRho = datos[:, 25]
OpRho = datos[:, 30]
CO2pRho = datos[:, 35]

HpRho_off = datos_Off[1000:1250, 20]

P_plasma = datos[:, 19]
presion_e = datos[:, 18]
presion_H = datos[:, 24]
presion_O2 = datos[:, 29]
presion_O = datos[:, 34]
presion_CO2 = datos[:, 39]

posicion_new = datos_new[1000:1250, 8:11]
presion_new = {
    "e-": datos_new[1000:1250, 18],
    "H+": datos_new[1000:1250, 24],
    "O+": datos_new[1000:1250, 34],
    "O2+": datos_new[1000:1250, 29],
    "CO2+": datos_new[1000:1250, 39],
}
densidad_new = {
    "e-": datos_new[1000:1250, 11],
    "H+": datos_new[1000:1250, 20],
    "O+": datos_new[1000:1250, 30],
    "O2+": datos_new[1000:1250, 25],
    "CO2+": datos_new[1000:1250, 35],
}
B_new = datos_new[1000:1250, 15:18]
J_new = datos_new[1000:1250, -3:]
velocidad_H_new = datos_new[1000:1250, 21:24]
b1_new = datos_new[1000:1250, 40:43]


P_heavy = presion_O + presion_O2 + presion_CO2
P_B = np.linalg.norm(B, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram = 1.67e-6 * HpRho * velocidad_H[:, 0] ** 2  # nPa
P_total = P_heavy + P_B + presion_e + P_ram + presion_H

P_heavy_new = presion_new["O+"] + presion_new["O2+"] + presion_new["CO2+"]
P_B_new = np.linalg.norm(B_new, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram_new = 1.67e-6 * densidad_new["H+"] * velocidad_H_new[:, 0] ** 2  # nPa
P_total_new = P_heavy_new + P_B_new + presion_new["H+"] + P_ram_new + presion_new["e-"]

beta = P_plasma / P_B

density_mean = [np.mean(HpRho[i : i + 50]) for i in range(len(HpRho) - 50)]
ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km
ion_length = np.append(ion_length, ion_length[-51:-1])


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
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection="3d")
# ax.set_xlabel(r"$X_{MSO} (R_m)$", fontsize=14)
# ax.set_ylabel(r"$Y_{MSO} (R_m)$", fontsize=14)
# ax.set_zlabel(r"$Z_{MSO} (R_m)$", fontsize=14)
# ax.plot_surface(
#     X1,
#     Y1,
#     Z1,
#     rstride=4,
#     cstride=4,
#     alpha=0.5,
#     # edgecolor="gray",
#     cmap=plt.get_cmap("Blues_r"),
# )
# ax.scatter(R[0], R[1], R[2], label="MAVEN", color="k", s=40)
#
# u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
# ax.plot_wireframe(
#     np.cos(u) * np.sin(v),
#     np.sin(u) * np.sin(v),
#     np.cos(v),
#     color="#c1440e",
#     linewidth=0.5,
# )
#
# plt.plot(x, y, z, label="Órbita", color="C1")
# plt.plot(
#     posicion_new[:, 0], posicion_new[:, 1], posicion_new[:, 2], label="new", color="C2"
# )
# # sc = ax.scatter(x, y, z, vmin=80, vmax=150, c=ion_length, cmap="magma")
# # sc = ax.scatter(x, y, z, c=np.log10(beta), vmin=-1.5, vmax=1.5, cmap="PiYG")
# # plt.colorbar(sc)
# # ax.set_title("log10(beta)")
# equal_axes(ax, X1, Y1, Z1)
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

limites_simu = [
    0.8978406562499996,
    1.0030125729166666,
    1.3053818333333331,
    1.3974072604166663,
]
# las posiciones de entrada y salida de la MPB en el eje x en RM (de la simu)
x1, x2, x3, x4 = limites_simu[0], limites_simu[1], limites_simu[2], limites_simu[3]

pos_MPB_simu = np.mean(limites_simu)

# plt.figure()
# plt.plot(posicion[:, 0] / 3390, np.linalg.norm(B_mag, axis=1), label="MAVEN")
# plt.axvline(posicion[i1, 0] / 3390, color="C2", linestyle="--", label="t1")
# plt.axvline(posicion[i2, 0] / 3390, color="C3", linestyle="--", label="t2")
# plt.axvline(posicion[i3, 0] / 3390, color="C4", linestyle="--", label="t3")
# plt.axvline(posicion[i4, 0] / 3390, color="C5", linestyle="--", label="t4")
# plt.xlabel("x (RM)")
# plt.ylabel("|B| (nT)")
# plt.title("Campo medido")
# plt.legend()
#
# plt.figure()
# plt.plot(x, np.linalg.norm(B, axis=1), label="campo")
# plt.axvline(x1, color="C2", linestyle="--", label="t1")
# plt.axvline(x2, color="C3", linestyle="--", label="t2")
# plt.axvline(x3, color="C4", linestyle="--", label="t3")
# plt.axvline(x4, color="C5", linestyle="--", label="t4")
# plt.xlabel("x (RM)")
# plt.ylabel("|B| (nT)")
# plt.title("Campo de la simulación Hall")
# plt.legend()


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
#         limites_simu = sorted(val)
#
#     print("Happy? Keyboard click for yes, mouse click for no.")
#     happy = plt.waitforbuttonpress()
#
# plt.show()


plt.figure()
ax1 = plt.subplot2grid((2, 1), (0, 0))
ax2 = plt.subplot2grid((2, 1), (1, 0))

ax1.plot(posicion[:, 0] / 3390, np.linalg.norm(B_mag, axis=1), label="MAVEN")
ax1.plot(posicion_new[:, 0], np.linalg.norm(B_new, axis=1), label="Simulation\n(HR)")
ax1.plot(
    posicion_new[:, 0], np.linalg.norm(b1_new, axis=1), label="Simulation b1\n(HR)"
)
ax1.plot(x, np.linalg.norm(B, axis=1), c="C3", linewidth=2, label="Simulation\n(LR)")
ax1.set_ylabel("|B| (nT)")
plt.setp(ax1.get_xticklabels(), visible=False)

ax2.plot(posicion[idx, 0] / 3390, proton_density, label="SWIA")
ax2.plot(posicion_new[:, 0], densidad_new["H+"], label="Simulation\n(HR)")
ax2.plot(x, HpRho, c="C2", linewidth=2, label="Simulation\n(LR)")
ax2.set_xlabel("x (RM)")
ax2.set_ylabel("Proton density (mp/cc)")

ax1.set_title("Quantities over the trajectory projected onto x")

for ax in [ax1, ax2]:
    ax.set_xlim(xmin=0)
    ax.grid()
    ax.legend(loc="upper left")

plt.show()


"""
Cálculo de J = n x (Bu-Bd)
"""

# x23 = (x2 - x3) * 3390e3  # esto no es el ancho posta, debería proyectar sobre la normal
# ancho_updown = 0.015 * 13000 / 3390
#
# inicio_up = donde(x, x1 + ancho_updown)
# fin_up = donde(x, x1)
# inicio_down = donde(x, x4)
# fin_down = donde(x, x4 - ancho_updown)
# inicio_MPB = donde(x, x2)
# fin_MPB = donde(x, x3)
#
# n2 = [0.856, -0.066, 0.512]
# B_upstream = np.mean(B[inicio_up:fin_up], axis=0)
# B_downstream = np.mean(B[inicio_down:fin_down], axis=0)
# J_s = np.cross(n2, (B_upstream - B_downstream)) / mu0  # nA/m
#
# J_integrado_x = np.trapz(
#     J[inicio_MPB:fin_MPB, 0], x[inicio_MPB:fin_MPB] * 3390e3
# )  # uA/m
# J_integrado_y = np.trapz(
#     J[inicio_MPB:fin_MPB, 0], y[inicio_MPB:fin_MPB] * 3390e3
# )  # uA/m
# J_integrado_z = np.trapz(
#     J[inicio_MPB:fin_MPB, 0], z[inicio_MPB:fin_MPB] * 3390e3
# )  # uA/m
#
# J_integrado = np.array([J_integrado_x, J_integrado_y, J_integrado_z]) * 1e-3  # mA/m
# np.linalg.norm(J_integrado)
#
# # Integramos el J de la simulación en x para obtener un Js
# # pero esto no sirve! La trayectoria no está en x
#
# # plt.figure()
# # # plt.plot(x, np.linalg.norm(J_integrado, axis=1) * 1e-3)
# # plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
# # plt.axvline(x=pos_MPB_simu, color="C3", linestyle="--", label="MPB simulacion")
# # plt.axhline(y=np.linalg.norm(J_s * 1e-6), color="C1", label="J = n x (Bu-Bd)")
# # plt.axhline(y=np.linalg.norm(j_maven), color="C2", label="J = n x (Bu-Bd)")
# # plt.xlim(xmin=1)
# # plt.xlabel("x MSO (RM)")
# # plt.ylabel("J_s  (mA/m)")
# # plt.legend()
# # plt.show()
#
#
# plt.figure()
# plt.plot(x, np.linalg.norm(J, axis=1) * 1e3, c="r", label="Hall")
# plt.plot(x, np.linalg.norm(J_off, axis=1) * 1e3, c="k", label="sin Hall")
# # plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
# # plt.axvline(x=pos_MPB_simu, color="C3", linestyle="--", label="MPB simulacion")
# # plt.axhline(y=282, color="C1", label="|J_v| MAVEN")
# plt.ylim(ymax=100)
# plt.xlabel("x (RM)")
# plt.ylabel("j_v (nA/m²)")
# plt.legend()
# plt.show()
#
# plt.plot(x, OpRho, label="O+")
# plt.plot(x, O2pRho, label="O2+")
# plt.plot(x, CO2pRho, label="CO2+")
# plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
# plt.axvline(x=pos_MPB_simu, color="C3", linestyle="--", label="MPB simulacion")
# plt.legend()
# plt.ylim(-1, 50)
# plt.xlim(xmin=0.9)
# plt.xlabel("x (RM)")
# plt.ylabel("rho (mp/cc)")
#
# """
# Análsis de las presiones
# """
#
# plt.figure()
# sc = plt.scatter(x, y, c=np.log10(beta), vmin=-4, vmax=4, s=35, cmap="coolwarm",)
# plt.xlabel("X (RM)")
# plt.ylabel("Y (RM)")
# plt.title("log10(beta)")
# plt.colorbar(sc)
#
#
# plt.figure()
# plt.scatter(x, P_total, label="P total")
# plt.scatter(x, presion_e, label="P electrónica")
# plt.scatter(x, P_B, label="P magnética")
# plt.scatter(x, P_ram, label="P ram")
# plt.scatter(x, presion_H, label="P protones")
# plt.scatter(x, P_heavy, label="P heavies")
# plt.axvspan(x2, x3, color="red", alpha=0.2, label="MPB")
# plt.xlim([0.5, 1.5])
# # plt.ylim([0, 1])
# plt.legend(loc="center left")
# plt.xlabel("x (RM)")
# plt.ylabel("presion (nPa)")
# plt.title("presiones sobre la trayectoria")
# plt.show()
