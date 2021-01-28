import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import MultiCursor
from importar_datos import importar_mag, importar_swia
import sys

sys.path.append("..")

from funciones_plot import equal_axes, onpick1
from funciones import diezmar, donde, corrientes

path = "../../../datos/simulacion_leonardo/"
datos = np.loadtxt(path + "simu_tray.sat", skiprows=2)

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
mag, t, B_mag, posicion = importar_mag(2016, "03", 16, 17.7, 19)
swia, t_swia, proton_density = importar_swia(2016, "03", 16, 17.7, 19)

# Datos del análisis de MAVEN
R = [1.082, -0.064, 0.515]
normal = [0.920, -0.302, 0.251]
j_maven = 282  # nA/m²
v_x = -13000  # km/h la velocidad de MAVEN en x

# hay 10 órbitas, cada una de 1620 puntos
pi = 10400
pf = 10800


"""
Las primeras órbitas, hasta el punto 10200 aprox, no tienen una MPB. El campo
en ellas es menor a 10 nT siempre.
"""

x_cut = x[pi:pf] * 3390
B_cut = B[pi:pf, :]

t = x_cut / v_x

"""
Comparación de B y de la densidad de protones
"""

idx = diezmar(t, t_swia)

plt.figure()
plt.plot(posicion[:, 0] / 3390, np.linalg.norm(B_mag, axis=1), label="MAVEN")
plt.plot(x_cut, np.linalg.norm(B_cut, axis=1), label="simulación")
plt.axvline(x=R[0], color="k", linestyle="--", label="cruce MAVEN")
plt.axvline(x=1.26, color="C3", linestyle="--", label="MPB simulacion")
plt.xlabel("x (RM)")
plt.ylabel("|B| (nT)")
plt.legend()

plt.figure()
plt.plot(posicion[idx, 0] / 3390, proton_density, label="swia")
plt.plot(x_cut, HpRho[pi:pf], label="simulación")
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
#         plt.plot(x_cut, np.linalg.norm(B_cut, axis=1))
#         ax1.set_ylabel("|B| (nT)")
#
#         ax5 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)
#         ax5.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
#         ax5.set_xlabel("x (RM)")
#         plt.plot(x_cut, HpRho[pi:pf])
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

inicio_up = donde(x_cut[200:], x1 + ancho_updown)
fin_up = donde(x_cut[200:], x1)
inicio_down = donde(x_cut[200:], x4)
fin_down = donde(x_cut[200:], x4 - ancho_updown)

B_upstream = np.mean(B_cut[inicio_up:fin_up], axis=0)
B_downstream = np.mean(B_cut[inicio_down:fin_down], axis=0)
mu = 4 * np.pi * 1e-7
J_v = np.cross(normal, (B_upstream - B_downstream)) / mu / x23

plt.figure()
plt.plot(x_cut, np.linalg.norm(J[pi:pf], axis=1) * 1e3)
plt.axvline(x=R[0], color="k")
plt.xlabel("x (RM)")
plt.ylabel("j (nA/m²)")

plt.show()

plt.plot(x_cut, OpRho[pi:pf], label="O+")
plt.plot(x_cut, O2pRho[pi:pf], label="O2+")
plt.plot(x_cut, CO2pRho[pi:pf], label="CO2+")
plt.axvline(x=R[0], color="k")
plt.legend()
plt.ylim(-1, 50)
plt.xlim(xmin=0.9)
plt.xlabel("x (RM)")
plt.ylabel("rho (mp/cc)")

plt.figure()
plt.plot(x_cut, presion_H[pi:pf], label="H+")
plt.plot(x_cut, presion_O[pi:pf], label="O+")
plt.plot(x_cut, presion_O2[pi:pf], label="O2+")
plt.plot(x_cut, presion_CO2[pi:pf], label="CO2+")
plt.axvline(x=R[0], color="k")
plt.legend()
plt.xlabel("x (RM)")
plt.ylabel("presion (nPa)")

plt.show()
