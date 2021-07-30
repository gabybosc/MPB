import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from sklearn.linear_model import LinearRegression
import sys

sys.path.append("..")

from funciones import donde


path = "../../../datos/simulacion_chuanfei/"
datos_hr = np.loadtxt(path + "ejex_new2_+.gz")  # high resolution
datos_lr = np.loadtxt(path + "ejex_new1_+.gz")  # low resolution

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
e_SI = 1.6e-19  # C
g = 3.7  # Mars surface gravity, m/s²

"""
r1 r2 rho Hrho Orho O2rho CO2rho H_vx H_vy H_vz
Bx By Bz b1x b1y b1z Pe P HP OP
O2P CO2P jx jy jz gradx grady gradz O_vx O_vy
O_vz O2_vx O2_vy O2_vz CO2_vx CO2_vy CO2_vz
"""

x = datos_hr[:, 0]  # RM
z = datos_hr[:, 1]  # RM
B_hr = datos_hr[:, 10:13]  # nT
b1_hr = datos_hr[:, 13:16]  # nT
J_hr = datos_hr[:, 22:25]  # ua/m2
grad_p_hr = datos_hr[:, 25:28]  # nPa/m
presion = {
    "e-": datos_hr[:, 16],
    "H+": datos_hr[:, 18],
    "O+": datos_hr[:, 19],
    "O2+": datos_hr[:, 20],
    "CO2+": datos_hr[:, 21],
}  # nPa
densities = {
    "e-": datos_hr[:, 2],
    "H+": datos_hr[:, 3],
    "O+": datos_hr[:, 4],
    "O2+": datos_hr[:, 5],
    "CO2+": datos_hr[:, 6],
}  # Mp/cc
velocities = {
    "H+": datos_hr[:, 7:10],
    "O+": datos_hr[:, 28:31],
    "O2+": datos_hr[:, 31:34],
    "CO2+": datos_hr[:, 34:],
}  # km/s


x_lr = datos_lr[:, 0]  # RM
z_lr = datos_lr[:, 1]  # RM
B_lr = datos_lr[:, 10:13]  # nT
b1_lr = datos_lr[:, 13:16]  # nT
J_lr = datos_lr[:, 22:25]  # ua/m2
grad_p_lr = datos_lr[:, 25:28]  # nPa/m

presion_lr = {
    "e-": datos_lr[:, 16],
    "H+": datos_lr[:, 18],
    "O+": datos_lr[:, 19],
    "O2+": datos_lr[:, 20],
    "CO2+": datos_lr[:, 21],
}  # nPa
densities_lr = {
    "e-": datos_lr[:, 2],
    "H+": datos_lr[:, 3],
    "O+": datos_lr[:, 4],
    "O2+": datos_lr[:, 5],
    "CO2+": datos_lr[:, 6],
}  # Mp/cc
velocities_lr = {
    "H+": datos_lr[:, 7:10],
    "O+": datos_lr[:, 28:31],
    "O2+": datos_lr[:, 31:34],
    "CO2+": datos_lr[:, 34:],
}  # km/s

# High res
B_norm = np.linalg.norm(B_hr, axis=1)

v_plus = np.zeros((len(x), 3))
for ion in ["H+", "O+", "O2+", "CO2+"]:
    for i in range(3):
        v_plus[:, i] += densities[ion] * velocities[ion][:, i] / densities["e-"]

v_SI_hr = v_plus * 1e3  # m/s
B_SI_hr = B_hr * 1e-9  # T
n_SI_hr = densities["H+"] * 1e6  # 1/m3
J_SI_hr = J_hr * 1e-6  # A/m2
grad_p_SI_hr = grad_p_hr * 1e-9

Ecv_hr = np.array([-np.cross(v_SI_hr[i], B_SI_hr[i]) for i in range(len(B_hr))])  # V/m
Ehall_hr = np.array(
    [
        1 / (e_SI * n_SI_hr[i]) * np.cross(J_SI_hr[i], B_SI_hr[i])
        for i in range(len(B_hr))
    ]
)
Ep_hr = np.array(
    [-1 / (e_SI * n_SI_hr[i]) * grad_p_SI_hr[i, :] for i in range(len(grad_p_hr))]
)


P_heavy_hr = presion["O+"] + presion["O2+"] + presion["CO2+"]
P_B_hr = np.linalg.norm(B_hr, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram_hr = 1.67e-6 * densities["H+"] * velocities["H+"][:, 0] ** 2  # nPa
P_total_hr = P_heavy_hr + P_B_hr + presion["e-"] + P_ram_hr + presion["H+"]

rho_heavies_hr = np.zeros(len(x))
for ion in ["O+", "O2+", "CO2+"]:
    rho_heavies_hr += densities[ion]

# new
B_norm_lr = np.linalg.norm(B_lr, axis=1)

v_plus_lr = np.zeros((len(x_lr), 3))
for ion in ["H+", "O+", "O2+", "CO2+"]:
    for i in range(3):
        v_plus_lr[:, i] += (
            densities_lr[ion] * velocities_lr[ion][:, i] / densities_lr["e-"]
        )

v_SI_lr = v_plus_lr * 1e3  # m/s
B_SI_lr = B_lr * 1e-9  # T
n_SI_lr = densities_lr["H+"] * 1e6  # 1/m3
J_SI_lr = J_lr * 1e-6  # A/m2
grad_p_SI_lr = grad_p_lr * 1e-9

Ecv_lr = np.array([-np.cross(v_SI_lr[i], B_SI_lr[i]) for i in range(len(B_lr))])
Ehall_lr = np.array(
    [
        1 / (e_SI * n_SI_lr[i]) * np.cross(J_SI_lr[i], B_SI_lr[i])
        for i in range(len(B_lr))
    ]
)
Ep_lr = np.array(
    [-1 / (e_SI * n_SI_lr[i]) * grad_p_SI_lr[i, :] for i in range(len(grad_p_lr))]
)  # V/m


P_heavy_lr = presion_lr["O+"] + presion_lr["O2+"] + presion_lr["CO2+"]
P_B_lr = np.linalg.norm(B_lr, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram_lr = 1.67e-6 * densities_lr["H+"] * velocities_lr["H+"][:, 0] ** 2  # nPa
P_total_lr = P_heavy_lr + P_B_lr + presion_lr["e-"] + P_ram_lr + presion_lr["H+"]

rho_heavies_lr = np.zeros(len(x_lr))
for ion in ["O+", "O2+", "CO2+"]:
    rho_heavies_lr += densities_lr[ion]


""" derivada de B en el eje x
Busco el punto máximo de la derivada y luego ajusto con una lineal el campo B
pero que pase por el punto máximo
"""
dBx_dx = np.gradient(b1_hr[:, 0], np.abs(x[0] - x[2]) * 3390e3)
dBy_dx = np.gradient(b1_hr[:, 1], np.abs(x[0] - x[2]) * 3390e3)
dBz_dx = np.gradient(b1_hr[:, 2], np.abs(x[0] - x[2]) * 3390e3)

dBdx = np.transpose(np.vstack((dBx_dx, dBy_dx, dBz_dx)))

plt.figure()
plt.plot(x, dBdx, ".")
plt.plot(x, np.linalg.norm(dBdx, axis=1), ".")
plt.legend(["x", "y", "z", "norm"])
plt.title("Derivada de B")
plt.show()

fig = plt.figure()
fig.subplots()
plt.title("Campo magnético y su derivada en x")

ax1 = plt.subplot2grid((2, 1), (0, 0))
plt.plot(x, B_norm, ".")
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("|B|")
ax1.grid()

ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)
ax2.plot(x, np.linalg.norm(dBdx, axis=1), ".")
ax2.set_ylabel("dB/dx")
ax2.set_xlabel("x (RM)")
ax2.grid()

multi = MultiCursor(fig.canvas, (ax1, ax2), color="black", lw=1)
plt.show()

"""
el ancho de la MPB es fácil porque hay un punto de inflexión bastante marcado
entonces tomo como límites los valores donde veo que la derivada forma la campana
"""

x_cut = x[donde(x, 1.174) : donde(x, 1.22)]
B_cut = B_norm[donde(x, 1.174) : donde(x, 1.22)]

coef = np.polynomial.polynomial.polyfit(x_cut, B_cut, deg=1)

MPR = np.mean(B_norm[donde(x, 1.16) : donde(x, 1.174)])

MPB_inicio = x[donde(coef[0] + x * coef[1], MPR)]
MPB_fin = x[donde(x, 1.22)]

plt.figure()
plt.title("Campo magnético y ajuste")
plt.plot(x, B_norm, ".")
plt.plot(x, coef[0] + x * coef[1])
plt.axhline(y=MPR, c="k")
plt.axhline(y=B_norm[donde(x, 1.22)], c="k")
plt.ylabel("|B|")
plt.xlabel("x (RM)")
plt.ylim(ymin=0)
plt.grid()
plt.show()

ancho_mpb = (MPB_fin - MPB_inicio) * 3390  # km
print(f"ancho de la mpb = {ancho_mpb}")

# longitud inercial de protones
paso = 20

density_mean_hr = [
    np.mean(densities["H+"][i : i + paso]) for i in range(len(densities["H+"]) - paso)
]

ion_length = 2.28e07 / np.sqrt(density_mean_hr) * 1e-5  # km

density_mean_lr = [
    np.mean(densities_lr["H+"][i : i + paso])
    for i in range(len(densities_lr["H+"]) - paso)
]

ion_length_lr = 2.28e07 / np.sqrt(density_mean_lr) * 1e-5  # km

fig, ax = plt.subplots()
ax2 = ax.twinx()

ax.plot(x[:-paso], ion_length, ".")
ax.set_ylabel("proton inertial length (km)")
ax2.set_ylabel("proton inertial length (RM)")

ax.set_ylim([50, 250])
ax2.set_ylim([50 / 3390, 250 / 3390])
# set an invisible artist to twin axes
# to prevent falling back to initial values on rescale events
ax2.plot([], [])
ax.set_xlabel("X (RM)")
ax.set_xlim([1.1, 1.7])

ax.grid()
plt.show()


"""
Giroradio
Habría que probar usando la v térmica en lugar de la v total que me da la simu
Pero el problema es que v térmica no es un vector
Puedo proyectar v_th en v_total y después hacer el v_perp pero vuelvo a tener el
problema de que la v total está mayoritariamente en la dirección de B
"""

# v_th_norm = np.sqrt(presion["H+"] * 1e-21 / (densities["H+"] * mp))  # km/s
#
# v_normalizado = np.array(
#     [
#         velocities["H+"][i, :] / np.linalg.norm(velocities["H+"][i, :])
#         for i in range(len(velocities["H+"]))
#     ]
# )
#
# v_th = np.array([v_th_norm[i] * v_normalizado[i, :] for i in range(len(v_th_norm))])
#
# B_medio = np.mean(B_hr, axis=0)
# B_medio_normalizado = B_medio / np.linalg.norm(B_medio)
#
# dot = np.dot(v_th, B_medio_normalizado)
# N = np.zeros((len(dot), len(B_medio)))
#
# for i in range(len(N)):
#     N[i, :] = dot[i] * B_medio_normalizado
#
# v_perp = np.mean(v_th - N, axis=0)  # v perp B
#
# # el giroradio entonces:
# rg = mp * np.linalg.norm(v_th_norm) / (e_SI * np.linalg.norm(B_medio)) * 1e9  # km
#
# print(f"El radio de Larmor calculado usando v_th es {rg:1.3g} km")
#
#
# # punto a punto
# N = np.zeros((len(dot), len(B_medio)))
# rg_pp = np.zeros(len(B_hr))
# for i in range(len(B_hr) - paso):
#     B_norm = B_hr[i, :] / np.linalg.norm(B_hr[i, :])
#     dot_rg = np.dot(velocities["H+"][i, :], B_norm)
#     v_perp = velocities["H+"][i, :] - dot_rg  # v perp B
#     rg_pp[i] = mp * np.linalg.norm(v_perp) / (e_SI * np.linalg.norm(B_hr)) * 1e9  # km
#
#
# fig, ax = plt.subplots()
# ax2 = ax.twinx()
#
# ax.plot(x, rg_pp, ".")
# ax.set_ylabel("proton gyroradius (km)")
# ax2.set_ylabel("proton gyroradius (RM)")
#
# # ax.set_ylim([50, 250])
# # ax2.set_ylim([50 / 3390, 250 / 3390])
# # set an invisible artist to twin axes
# # to prevent falling back to initial values on rescale events
# ax2.plot([], [])
# ax.set_xlabel("X (RM)")
# ax.set_xlim([1.1, 1.7])
#
# ax.grid()
# plt.show()

# plots

# Presion
plt.figure()
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))

ax1.scatter(x_lr, presion_lr["e-"])
ax1.scatter(x_lr, P_heavy_lr)
ax1.scatter(x_lr, P_B_lr)
ax1.scatter(x_lr, P_ram_lr)
ax1.scatter(x_lr, P_total_lr)
ax1.scatter(x_lr, presion_lr["H+"])

ax2.scatter(x, presion["e-"], label="e-")
ax2.scatter(x, P_heavy_hr, label="P heavies")
ax2.scatter(x, P_B_hr, label="mag")
ax2.scatter(x, P_ram_hr, label="ram")
ax2.scatter(x, P_total_hr, label="total")
ax2.scatter(x, presion["H+"], label="H+", marker=".")

plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_title("Hall LR")
ax2.set_title("Hall HR")
ax1.set_ylabel("P (nPa)")
ax1.set_xlabel("x (RM)")

for ax in [ax1, ax2]:
    ax.set_xlabel("x (RM)")
    ax.set_ylim([-0.1, 1.5])
    ax.set_xlim([1.1, 1.7])
    ax.grid()

ax2.legend(loc="upper right")
plt.show()


# Densidades

plt.figure()
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))

ax1.scatter(x_lr, densities_lr["H+"])
ax1.scatter(x_lr, rho_heavies_lr)
ax1.scatter(x_lr, densities_lr["e-"])

ax2.scatter(x, densities["H+"], label="H+")
ax2.scatter(x, rho_heavies_hr, label="heavies")
ax2.scatter(x, densities["e-"], label="e-")

plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_title("Hall LR")
ax2.set_title("Hall HR")
ax1.set_ylabel("Particle density (mp/cc)")

ax2.legend()

for ax in [ax1, ax2]:
    ax.set_xlabel("x (RM)")
    ax.set_ylim([-0.1, 10])
    ax.set_xlim([1.1, 1.7])
    ax.grid()

plt.show()

# Corrientes y B

plt.figure()
ax2 = plt.subplot2grid((2, 1), (0, 0))
ax3 = plt.subplot2grid((2, 1), (1, 0))

ax2.plot(x, B_hr, ".")
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_ylim([-10, 50])
ax2.set_title("Magnetic field and volume current density (high resolution)")

ax3.plot(x, J_hr, ".")
ax3.set_ylabel("J (uA/m²)")
ax3.set_ylim([-0.1, 0.2])
ax3.set_xlabel("x (RM)")

for ax in [ax2, ax3]:
    ax.grid()
    ax.set_xlim([1.1, 1.7])
    ax.axvspan(xmin=MPB_inicio, xmax=MPB_fin, facecolor="b", alpha=0.2)

ax3.legend(["x", "y", "z", "MPB"])
plt.show()

# B low vs B hi
plt.figure()
plt.plot(x, np.linalg.norm(B_hr, axis=1), ".", label="HR")
plt.plot(x_lr, np.linalg.norm(B_lr, axis=1), ".", label="LR")
plt.legend()
plt.show()

plt.figure()
ax2 = plt.subplot2grid((2, 1), (0, 0))
ax3 = plt.subplot2grid((2, 1), (1, 0))

ax2.plot(x, B_hr, ".")
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (HR)")
ax2.set_ylim([-10, 50])
ax2.set_title("Magnetic field HR vs LR")

ax3.plot(x_lr, B_lr, ".")
ax3.set_ylabel("B (LR))")
ax3.set_ylim([-10, 50])
ax3.set_xlabel("x (RM)")

for ax in [ax2, ax3]:
    ax.grid()
    # ax.set_xlim([1.1, 1.7])
    ax.axvspan(xmin=MPB_inicio, xmax=MPB_fin, facecolor="b", alpha=0.2)

ax3.legend(["x", "y", "z", "MPB"])
plt.show()

# E norm
plt.figure()
plt.plot(x, np.linalg.norm(Ecv_hr, axis=1) * 1e3, ".", label="Ecv")
plt.plot(x, np.linalg.norm(Ehall_hr, axis=1) * 1e3, ".", label="Ehall")
plt.plot(x, np.linalg.norm(Ep_hr, axis=1) * 1e3, ".", label="Ep")
plt.axvspan(xmin=MPB_inicio, xmax=MPB_fin, facecolor="b", alpha=0.2, label="MPB")

plt.xlabel("x (RM)")
plt.ylabel("E (mV/m)")
plt.title("Ohm's Law terms (HR)")
plt.ylim([-0.1, 6])
plt.xlim([1.1, 1.7])
plt.grid()
plt.legend()

plt.show()


# E componentes
plt.figure()
ax4 = plt.subplot2grid((3, 1), (0, 0))
ax5 = plt.subplot2grid((3, 1), (1, 0))
ax6 = plt.subplot2grid((3, 1), (2, 0))

ax4.set_title("Electric fields")
ax4.plot(x, Ecv_hr * 1e3, ".")
ax4.set_ylabel("Ecv (mV/m)")
plt.setp(ax4.get_xticklabels(), visible=False)

ax5.plot(x, Ehall_hr * 1e3, ".")
ax5.set_ylabel("Ehall (mV/m)")
plt.setp(ax5.get_xticklabels(), visible=False)

ax6.plot(x, Ep_hr * 1e3, ".")
ax6.set_ylabel("Ep (mV/m)")

for ax in [ax6, ax4, ax5]:
    ax.set_ylim([-5, 5])
    ax.set_xlim([1.1, 1.7])
    ax.grid()
    ax.axvspan(xmin=MPB_inicio, xmax=MPB_fin, facecolor="b", alpha=0.2)

ax6.set_xlabel("x (RM)")
ax6.legend(["x", "y", "z", "MPB"])


plt.show()
