import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import MultiCursor
from importar_datos import importar_mag, importar_swia
import sys

sys.path.append("..")

from funciones_plot import equal_axes, onpick1
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

datos = datos_enteros[pi:pf]

x = datos[:, 8]
y = datos[:, 9]
z = datos[:, 10]

B = datos[:, 15:18]
J = datos[:, -3:]

rho = datos[:, 11]
HpRho = datos[:, 19]
O2pRho = datos[:, 24]
OpRho = datos[:, 29]
CO2pRho = datos[:, 34]

presion_H = datos[:, 23]
presion_O2 = datos[:, 28]
presion_O = datos[:, 33]
presion_CO2 = datos[:, 38]

# Datos de MAVEN
mag, t, B_mag, posicion = importar_mag(2016, "03", 16, 18, 19)
swia, t_swia, proton_density = importar_swia(2016, "03", 16, 17.7, 19)

# Datos del análisis de MAVEN
R = [1.082, -0.064, 0.515]
normal = [0.920, -0.302, 0.251]
j_maven = 282  # nA/m²
v_x = -13000  # km/h la velocidad de MAVEN en x

x_km = x * 3390

deltat = np.abs((x_km[200] - x_km[201]) / v_x)


t1, t2, t3, t4 = 18.2167, 18.2204, 18.235, 18.2476
t_up = t1 - 0.015
t_down = t4 + 0.015

# Suponiendo que MAVEN y los datos están en el mismo lugar a las 17.7 h:
# 17.7h debería ser un punto en el viento solar

i = donde(x_km, posicion[0, 0])
x = x[i:]
x_km = x_km[i:]
B_cut = B[i:]
t_simu = np.array([18 + deltat * i for i in range(len(x_km))])

inicio_MPB = donde(t_simu, t1)
fin_MPB = donde(t_simu, t4)
inicio_up = donde(t_simu, t_up)
fin_down = donde(t_simu, t_down)


plt.figure()
ax1 = plt.subplot2grid((2, 1), (0, 0))
ax1.scatter(t_simu, x_km)
ax1.set_title('Simulación')
ax2 = plt.subplot2grid((2, 1), (1, 0))
ax2.scatter(t, posicion[:, 0])
ax2.set_title('MAVEN')
for ax in [ax1, ax2]:
    ax.axvspan(
        xmin=t_simu[inicio_MPB],
        xmax=t_simu[fin_MPB],
        facecolor="#FFC300",
        alpha=0.5,
        label='MPB'
    )
    ax.axvspan(
        xmin=t_simu[inicio_up],
        xmax=t_simu[inicio_MPB],
        facecolor="#581845",
        alpha=0.5,
        label='Upstream'
    )
    ax.axvspan(
        xmin=t_simu[fin_MPB],
        xmax=t_simu[fin_down],
        facecolor="#C70039",
        alpha=0.5,
        label='Downstream'
    )
    ax.set_xlim(18, 19)
plt.legend()
plt.show()


"""
Comparación de B y de la densidad de protones
"""

# idx = diezmar(t, t_swia)
#
#
# plt.figure()
# plt.plot(t, np.linalg.norm(B_mag, axis=1), label="MAVEN")
# plt.plot(t_simu, np.linalg.norm(B_cut, axis=1), label="simulación")
# plt.xlabel("x (RM)")
# plt.ylabel("|B| (nT)")
# plt.legend()
#
# plt.figure()
# plt.plot(posicion[:, 0] / 3390, np.linalg.norm(B_mag, axis=1), label="MAVEN")
# plt.plot(x, np.linalg.norm(B, axis=1), label="simulación")
# plt.plot(x - 0.178, np.linalg.norm(B, axis=1), label="simulación")
# plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
# plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
# plt.xlabel("x (RM)")
# plt.ylabel("|B| (nT)")
# plt.legend()
#
# plt.figure()
# plt.plot(posicion[idx, 0] / 3390, proton_density, label="swia")
# plt.plot(x, HpRho[pi:pf], label="simulación")
# plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
# plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
# plt.ylim(ymax=30)
# plt.xlabel("x (RM)")
# plt.ylabel("Densidad de protones")
# plt.legend()
#
# plt.show()

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
#         plt.plot(x, HpRho[pi:pf])
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
# x1 = 1.393042495053138
# x2 = 1.3234990984508788
# x3 = 1.2017981543969254
# x4 = 1.1474673758014104
#
# x23 = (x2 - x3) * 3390e3  # ancho en m
# ancho_updown = 0.015 * 13000 / 3390
#
# inicio_up = donde(x[200:], x1 + ancho_updown)
# fin_up = donde(x[200:], x1)
# inicio_down = donde(x[200:], x4)
# fin_down = donde(x[200:], x4 - ancho_updown)
#
# B_upstream = np.mean(B[inicio_up:fin_up], axis=0)
# B_downstream = np.mean(B[inicio_down:fin_down], axis=0)
# mu = 4 * np.pi * 1e-7
# J_v = np.cross(normal, (B_upstream - B_downstream)) / mu / x23
