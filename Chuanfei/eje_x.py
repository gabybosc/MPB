import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
import sys

sys.path.append("..")

from funciones import donde


path = "../../../datos/simulacion_chuanfei/"
datos_on = np.loadtxt(path + "ejex_HallOn_+.gz")
datos_off = np.loadtxt(path + "ejex_HallOff_+.gz")

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

x = datos_on[:, 0]  # RM
z = datos_on[:, 1]  # RM
B_on = datos_on[:, 10:13]  # nT
b1_on = datos_on[:, 13:16]  # nT
J_on = datos_on[:, 22:25]  # ua/m2
grad_p_on = datos_on[:, 25:28]  # nPa/m
presion = {
    "e-": datos_on[:, 16],
    "H+": datos_on[:, 18],
    "O+": datos_on[:, 19],
    "O2+": datos_on[:, 20],
    "CO2+": datos_on[:, 21],
}  # nPa
densities = {
    "e-": datos_on[:, 2],
    "H+": datos_on[:, 3],
    "O+": datos_on[:, 4],
    "O2+": datos_on[:, 5],
    "CO2+": datos_on[:, 6],
}  # Mp/cc
velocities = {
    "H+": datos_on[:, 7:10],
    "O+": datos_on[:, 28:31],
    "O2+": datos_on[:, 31:34],
    "CO2+": datos_on[:, 34:],
}  # km/s


Rho_off = datos_off[:, 2]  # mp/cc
HpRho_off = datos_off[:, 3]  # Mp/cc
OpRho_off = datos_off[:, 4]  # Mp/cc
O2pRho_off = datos_off[:, 5]  # Mp/cc
CO2pRho_off = datos_off[:, 6]  # Mp/cc
B_off = datos_off[:, 10:13]  # nT
b1_off = datos_off[:, 13:16]  # nT
Pe_off = datos_off[:, 16]  # nPa
HP_off = datos_off[:, 18]  # nPa
OP_off = datos_off[:, 19]  # nPa
O2P_off = datos_off[:, 20]  # nPa
CO2P_off = datos_off[:, 21]  # nPa
J_off = datos_off[:, 22:25]  # ua/m2
grad_p_off = datos_off[:, 25:28]  # nPa/m
v_i_off = datos_off[:, 7:10]  # km/s
v_O_off = datos_off[:, 28:31]  # km/s
v_O2_off = datos_off[:, 31:34]  # km/s
v_CO2_off = datos_off[:, 34:]  # km/s

# On
B_norm = np.linalg.norm(B_on, axis=1)

v_plus = np.zeros((len(x), 3))
for ion in ["H+", "O+", "O2+", "CO2+"]:
    for i in range(3):
        v_plus[:, i] += densities[ion] * velocities[ion][:, i] / densities["e-"]

v_SI_on = v_plus * 1e3  # m/s
B_SI_on = B_on * 1e-9  # T
n_SI_on = densities["H+"] * 1e6  # 1/m3
J_SI_on = J_on * 1e-6  # A/m2
grad_p_SI_on = grad_p_on * 1e-9

Ecv_on = np.array([-np.cross(v_SI_on[i], B_SI_on[i]) for i in range(len(B_on))])  # V/m
Ehall_on = np.array(
    [
        1 / (e_SI * n_SI_on[i]) * np.cross(J_SI_on[i], B_SI_on[i])
        for i in range(len(B_on))
    ]
)
Ep_on = np.array(
    [-1 / (e_SI * n_SI_on[i]) * grad_p_SI_on[i, :] for i in range(len(grad_p_on))]
)


P_heavy_on = presion["O+"] + presion["O2+"] + presion["CO2+"]
P_B_on = np.linalg.norm(B_on, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram_on = 1.67e-6 * densities["H+"] * velocities["H+"][:, 0] ** 2  # nPa
P_total_on = P_heavy_on + P_B_on + presion["e-"] + P_ram_on + presion["H+"]

rho_heavies_on = np.zeros(len(x))
for ion in ["O+", "O2+", "CO2+"]:
    rho_heavies_on += densities[ion]

# Off
v_SI_off = v_i_off * 1e3  # m/s
B_SI_off = B_off * 1e-9  # T
n_SI_off = HpRho_off * 1e6  # 1/m3
J_SI_off = J_off * 1e-6  # A/m2
grad_p_SI_off = grad_p_off * 1e-9

Ecv_off = np.array(
    [-np.cross(v_SI_off[i], B_SI_off[i]) for i in range(len(B_off))]
)  # V/m

Ehall_off = np.array(
    [
        1 / (e_SI * n_SI_off[i]) * np.cross(J_SI_off[i], B_SI_off[i])
        for i in range(len(B_off))
    ]
)
Ep_off = np.array(
    [-1 / (e_SI * n_SI_off[i]) * grad_p_SI_off[i, :] for i in range(len(grad_p_off))]
)

P_heavy_off = OP_off + O2P_off + CO2P_off
P_B_off = np.linalg.norm(B_off, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram_off = 1.67e-6 * HpRho_off * v_i_off[:, 0] ** 2  # nPa
P_total_off = P_heavy_off + P_B_off + Pe_off + P_ram_off + HP_off


v_total_off = (v_i_off + v_O_off + v_O2_off + v_CO2_off) * 1e3  # m/s

Ecv_off_total = np.array(
    [-np.cross(v_total_off[i], B_SI_off[i]) for i in range(len(B_off))]
)  # V/m

rho_heavies_off = OpRho_off + O2pRho_off + CO2pRho_off

""" derivada de B en el eje x
Busco el punto máximo de la derivada y luego ajusto con una lineal el campo B
pero que pase por el punto máximo
"""
dBx_dx = np.gradient(b1_on[:, 0], np.abs(x[0] - x[2]) * 3390e3)
dBy_dx = np.gradient(b1_on[:, 1], np.abs(x[0] - x[2]) * 3390e3)
dBz_dx = np.gradient(b1_on[:, 2], np.abs(x[0] - x[2]) * 3390e3)

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

x_cut = x[donde(x, 1.20464) : donde(x, 1.26368)]
B_cut = B_norm[donde(x, 1.20464) : donde(x, 1.26368)]

coef = np.polynomial.polynomial.polyfit(x_cut, B_cut, deg=1)

MPR = np.mean(B_norm[donde(x, 1.17) : donde(x, 1.20)])
MS = np.mean(B_norm[donde(x, 1.26) : donde(x, 1.4)])

plt.figure()
plt.title("Campo magnético y ajuste")
plt.plot(x, B_norm, ".")
plt.plot(x, coef[0] + x * coef[1])
plt.axhline(y=MPR, c="k")
plt.axhline(y=MS, c="k")
plt.ylabel("|B|")
plt.xlabel("x (RM)")
plt.ylim(ymin=0)
plt.grid()
plt.show()

# longitud inercial de protones
paso = 20

density_mean_on = [
    np.mean(densities["H+"][i : i + paso]) for i in range(len(densities["H+"]) - paso)
]

ion_length = 2.28e07 / np.sqrt(density_mean_on) * 1e-5  # km

density_mean_off = [
    np.mean(HpRho_off[i : i + paso]) for i in range(len(HpRho_off) - paso)
]

ion_length_off = 2.28e07 / np.sqrt(density_mean_off) * 1e-5  # km

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
Giroradio térmico
"""
v_th = np.sqrt(
    np.mean(presion["H+"][MPB_inicio:MPB_fin], axis=0)
    * 1e-21
    / (np.mean(densities["H+"][MPB_inicio:MPB_fin], axis=0) * mp)
)  # km/s

B_medio = np.mean(B_hr[MPB_inicio:MPB_fin], axis=0)

rg = mp * np.linalg.norm(v_th) / (e_SI * np.linalg.norm(B_medio)) * 1e9  # km

print(f"El giroradio térmico es {rg:1.3g} km")

# punto a punto
v_th = np.sqrt(presion["H+"] * 1e-21 / (densities["H+"] * mp))  # km/s

rg = [
    mp * v_th[i] / (e_SI * np.linalg.norm(B_on[i, :])) * 1e9 for i in range(len(B_on))
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
ax.set_xlim([1.1, 1.7])

ax.grid()
plt.show()


# plots

# Presion
plt.figure()
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))

ax1.scatter(x, Pe_off)
ax1.scatter(x, P_heavy_off)
ax1.scatter(x, P_B_off)
ax1.scatter(x, P_ram_off)
ax1.scatter(x, P_total_off)
ax1.scatter(x, HP_off)

ax2.scatter(x, presion["e-"], label="e-")
ax2.scatter(x, P_heavy_on, label="P heavies")
ax2.scatter(x, P_B_on, label="mag")
ax2.scatter(x, P_ram_on, label="ram")
ax2.scatter(x, P_total_on, label="total")
ax2.scatter(x, presion["H+"], label="H+", marker=".")

plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_title("Hall Off")
ax2.set_title("Hall On")
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

ax1.scatter(x, HpRho_off)
ax1.scatter(x, rho_heavies_off)
ax1.scatter(x, Rho_off)

ax2.scatter(x, densities["H+"], label="H+")
ax2.scatter(x, rho_heavies_on, label="heavies")
ax2.scatter(x, densities["e-"], label="e-")

plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_title("Hall Off")
ax2.set_title("Hall On")
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
plt.title("Magnetic field and volume current density")
ax2 = plt.subplot2grid((2, 1), (0, 0))
ax3 = plt.subplot2grid((2, 1), (1, 0))

ax2.plot(x, B_on, ".")
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_ylim([-10, 50])

ax3.plot(x, J_on, ".")
ax3.legend(["x", "y", "z"])
ax3.set_ylabel("J (uA/m²)")
ax3.set_ylim([-0.1, 0.2])
ax3.set_xlabel("x (RM)")

for ax in [ax2, ax3]:
    ax.grid()
    ax.set_xlim([1.1, 1.7])

plt.show()


# E norm
plt.figure()
plt.plot(x, np.linalg.norm(Ecv_on, axis=1) * 1e3, ".", label="Ecv")
plt.plot(x, np.linalg.norm(Ehall_on, axis=1) * 1e3, ".", label="Ehall")
plt.plot(x, np.linalg.norm(Ep_on, axis=1) * 1e3, ".", label="Ep")
plt.legend()

plt.xlabel("x (RM)")
plt.ylabel("E (mV/m)")
plt.ylim([-0.1, 10])
plt.xlim([1.1, 1.7])
plt.grid()

plt.show()


# E componentes
plt.figure()
ax4 = plt.subplot2grid((3, 1), (0, 0))
ax5 = plt.subplot2grid((3, 1), (1, 0))
ax6 = plt.subplot2grid((3, 1), (2, 0))

ax4.set_title("Electric fields")
ax4.plot(x, Ecv_on * 1e3, ".")
ax4.set_ylabel("Ecv (mV/m)")
plt.setp(ax4.get_xticklabels(), visible=False)

ax5.plot(x, Ehall_on * 1e3, ".")
ax5.set_ylabel("Ehall (mV/m)")
plt.setp(ax5.get_xticklabels(), visible=False)

ax6.plot(x, Ep_on * 1e3, ".")
ax6.set_ylabel("Ep (mV/m)")
ax6.legend(["x", "y", "z"])

for ax in [ax6, ax4, ax5]:
    ax.set_ylim([-5, 5])
    ax.set_xlim([1.1, 1.7])
    ax.grid()

ax6.set_xlabel("x (RM)")


plt.show()
