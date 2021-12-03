import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from cycler import cycler
import sys

sys.path.append("..")

from funciones import donde

plt.rcParams.update({"font.size": 12})
plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

path = "../../../datos/simulacion_chuanfei/"
datos_hr = np.loadtxt(path + "ejex_new3_+.gz")  # high resolution
datos_lr = np.loadtxt(path + "ejex_new1_+.gz")  # low resolution

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
e_SI = 1.6e-19  # C
g = 3.7  # Mars surface gravity, m/s²

"""
r1 r2 rho Hrho Orho O2rho CO2rho H_vx H_vy H_vz
Bx By Bz b1x b1y b1z Pe P HP OP
O2P CO2P jx jy jz gradx grady gradz O_vx O_vy
O_vz O2_vx O2_vy O2_vz CO2_vx CO2_vy CO2_vz e_vx e_vy e_vz
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
    "CO2+": datos_hr[:, 34:37],
    "e-": datos_hr[:, 37:],
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
P_B_hr = np.linalg.norm(b1_hr, axis=1) ** 2 * 1e-9 / (2 * mu0)  # sin corticales
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


ji_z = e_SI * densities["H+"] * velocities["H+"][:, 2] * 1e18  # nA/m²
je_z = e_SI * densities["e-"] * velocities["e-"][:, 2] * 1e18  # nA/m²

ji_y = e_SI * densities["H+"] * velocities["H+"][:, 1] * 1e18  # nA/m²
je_y = e_SI * densities["e-"] * velocities["e-"][:, 1] * 1e18  # nA/m²

je_y = e_SI * densities["e-"] * velocities["e-"][:, 1] * 1e18  # nA/m²

plt.plot(x, J_hr[:, 2] * 1e3, ".", label="corriente")
plt.plot(x, je_z, label="electrones")
plt.plot(x, ji_z, label="protones")
plt.plot(x, ji_z - je_z, label="resta")
plt.legend()
plt.xlabel("x")
plt.ylabel("Jz")
plt.xlim([1.15, 1.6])
plt.ylim([-100, 100])

plt.figure()
plt.plot(x, J_hr[:, 1] * 1e3, ".", label="corriente")
plt.plot(x, je_y, label="electrones")
plt.plot(x, ji_y, label="protones")
plt.plot(x, ji_y - je_y, label="resta")
plt.legend()
plt.xlabel("x")
plt.ylabel("Jy")
plt.xlim([1.15, 1.6])
plt.ylim([0, 100])
plt.show()


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

MPB_inicio = donde(coef[0] + x * coef[1], MPR)
MPB_fin = donde(x, 1.22)

plt.figure()


plt.title("MPB thickness")
plt.plot(x, B_norm, ".")
plt.plot(x, coef[0] + x * coef[1], label="fit")
plt.axhline(y=MPR, c="k")
plt.axhline(y=B_norm[donde(x, 1.22)], c="k")
plt.ylabel("|B|")
plt.xlabel("x (RM)")
plt.ylim(ymin=0)
plt.grid()

plt.show()

ancho_mpb = (x[MPB_fin] - x[MPB_inicio]) * 3390  # km
print(f"ancho de la mpb = {ancho_mpb:.3g} km")

j_media = np.mean(J_hr[MPB_inicio:MPB_fin]) * 1e3
print(f"la corriente media en la MPB de la simu es {np.abs(j_media):.3g} nA/m²")

# normal = [0.920, -0.302, 0.251]
normal = [1, 0, 0]
Bup = np.mean(B_hr[donde(x, 1.16) : donde(x, 1.174)], axis=0) * 1e-9
Bdown = np.mean(B_hr[donde(x, 1.22) : donde(x, 1.234)], axis=0) * 1e-9
J_salto = 1 / (mu0 * ancho_mpb * 1e3) * np.cross(normal, Bup - Bdown) * 1e9

print(
    f"la corriente con la condición de salto en la MPB de la simu es {np.linalg.norm(J_salto):.3g} nA/m²"
)
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

f, (ax1, ax) = plt.subplots(2, 1, sharex=True)
ax1.plot(x, B_norm, ".")
ax1.plot(x, coef[0] + x * coef[1], label="fit")
ax1.axhline(y=MPR, c="k")
ax1.axhline(y=B_norm[donde(x, 1.22)], c="k")
ax1.set_title("MPB thickness")
ax1.set_ylabel("|B|")
ax1.set_ylim([0, 50])
ax1.grid()

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
ax.set_xlim([1.15, 1.6])

ax.grid()
plt.show()

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
ax.set_xlim([1.15, 1.6])

ax.grid()
plt.show()


"""
Giroradio térmico en la MPB
"""

v_th = np.sqrt(
    np.mean(presion["H+"][MPB_inicio:MPB_fin], axis=0)
    * 1e-21
    / (np.mean(densities["H+"][MPB_inicio:MPB_fin], axis=0) * mp)
)  # km/s

B_medio = np.mean(B_hr[MPB_inicio:MPB_fin], axis=0)

rg = mp * v_th / (e_SI * np.linalg.norm(B_medio)) * 1e9  # km

print(f"El giroradio térmico medio en la MPB es {rg:1.3g} km")

# punto a punto
v_th = np.sqrt(presion["H+"] * 1e-21 / (densities["H+"] * mp))  # km/s

rg = [
    mp * v_th[i] / (e_SI * np.linalg.norm(B_hr[i, :])) * 1e9 for i in range(len(B_hr))
]  # km

fig, ax = plt.subplots()
ax2 = ax.twinx()

ax.plot(x, rg, ".")
ax.set_ylabel("proton gyroradius (km)")
ax2.set_ylabel("proton gyroradius (RM)")

ax.set_ylim([0, 250])
ax2.set_ylim([0 / 3390, 250 / 3390])
# set an invisible artist to twin axes
# to prevent falling back to initial values on rescale events
ax2.plot([], [])
ax.set_xlabel("X (RM)")
ax.set_xlim([1.15, 1.6])

ax.grid()
plt.show()

# plots

# Presion
plt.figure()

# plt.scatter(x, presion["H+"], label="H+")
# plt.scatter(x, presion["e-"], label="e-")
# plt.scatter(x, P_heavy_hr, label="heavies")
plt.scatter(x, presion["H+"] + presion["e-"] + P_heavy_hr, label="thermal")
plt.scatter(x, P_B_hr, label="magnetic")
plt.scatter(x, P_ram_hr, label="dynamic")
plt.scatter(x, P_total_hr, label="total")

plt.title("Simulation pressure results")
plt.ylabel("P (nPa)")
plt.xlabel("x (RM)")
plt.ylim([-0.1, 1])
plt.xlim([1.15, 1.6])
plt.grid()
plt.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2, label="MPB")

plt.legend(loc="best")
plt.show()


# Densidades

plt.figure()

plt.scatter(x, densities["H+"], label="H+")
plt.scatter(x, rho_heavies_hr, label="heavies")
plt.scatter(x, densities["e-"], label="e-")
plt.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2, label="MPB")

plt.title("Simulation density results")
plt.ylabel("Particle density (mp/cc)")
plt.xlabel("x (RM)")
plt.legend()
plt.ylim([-0.1, 20])
plt.xlim([1.15, 1.6])
plt.grid()

plt.show()

# Corrientes y B

plt.figure()
ax2 = plt.subplot2grid((2, 1), (0, 0))
ax3 = plt.subplot2grid((2, 1), (1, 0))

ax2.plot(x, B_hr, ".")
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_ylim([-10, 50])
ax2.set_title("Magnetic field and volume current density")

ax3.plot(x, J_hr * 1e3, ".")
ax3.set_ylabel("J (nA/m²)")
ax3.set_ylim([-100, 100])
ax3.set_xlabel("x (RM)")

for ax in [ax2, ax3]:
    ax.grid()
    ax.set_xlim([1.15, 1.6])
    ax.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2)

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
    # ax.set_xlim([1.15, 1.6])
    ax.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2)

ax3.legend(["x", "y", "z", "MPB"])
plt.show()

# E norm
plt.figure()
plt.plot(x, np.linalg.norm(Ecv_hr, axis=1) * 1e3, ".", label="|E| cv")
plt.plot(x, np.linalg.norm(Ehall_hr, axis=1) * 1e3, ".", label="|E| Hall")
plt.plot(x, np.linalg.norm(Ep_hr, axis=1) * 1e3, ".", label=r"|E| p$_e$")
plt.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2, label="MPB")

plt.xlabel("x (RM)")
plt.ylabel("|E| (mV/m)")
plt.title("Ohm's Law terms")
plt.ylim([-0.1, 4])
plt.xlim([1.15, 1.6])
plt.grid()
plt.legend()

plt.show()


# E componentes
plt.figure()
ax4 = plt.subplot2grid((3, 1), (0, 0))
ax5 = plt.subplot2grid((3, 1), (1, 0))
ax6 = plt.subplot2grid((3, 1), (2, 0))

ax4.set_title("Electric field components")
ax4.plot(x, Ecv_hr * 1e3, ".")
ax4.set_ylabel("Ecv (mV/m)")
plt.setp(ax4.get_xticklabels(), visible=False)

ax5.plot(x, Ehall_hr * 1e3, ".")
ax5.set_ylabel("Ehall (mV/m)")
plt.setp(ax5.get_xticklabels(), visible=False)

ax6.plot(x, Ep_hr * 1e3, ".")
ax6.set_ylabel("Ep (mV/m)")

for ax in [ax6, ax4, ax5]:
    ax.set_ylim([-4, 4])
    ax.set_xlim([1.15, 1.6])
    ax.grid()
    ax.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2)

ax6.set_xlabel("x (RM)")
ax6.legend(["x", "y", "z", "MPB"])


plt.show()

plt.figure()
plt.plot(x, Ehall_hr * 1e3, ".")
plt.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2)
plt.ylabel("E Hall (mV/m)")
plt.xlabel("x (RM)")
plt.title("E Hall components")
plt.legend(["x", "y", "z", "MPB"])
plt.grid()
plt.ylim([-1, 4])
plt.xlim([1.15, 1.6])
plt.show()
plt.xlim([1.15, 1.6])
plt.show()
