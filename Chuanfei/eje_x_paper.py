import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
import sys
from cycler import cycler

sys.path.append("..")

from funciones import donde

plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

path = "../../../datos/simulacion_chuanfei/"
datos = np.loadtxt(path + "nueva_simu/ejex_nuevasimu+.gz")  # high resolution

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

x = datos[:, 0]  # RM
z = datos[:, 1]  # RM
B = datos[:, 10:13]  # nT
b1 = datos[:, 13:16]  # nT
J = datos[:, 22:25]  # ua/m2
grad_p = datos[:, 25:28]  # nPa/m
presion = {
    "e-": datos[:, 16],
    "H+": datos[:, 18],
    "O+": datos[:, 19],
    "O2+": datos[:, 20],
    "CO2+": datos[:, 21],
}  # nPa
densities = {
    "e-": datos[:, 2],
    "H+": datos[:, 3],
    "O+": datos[:, 4],
    "O2+": datos[:, 5],
    "CO2+": datos[:, 6],
}  # Mp/cc
velocities = {
    "H+": datos[:, 7:10],
    "O+": datos[:, 28:31],
    "O2+": datos[:, 31:34],
    "CO2+": datos[:, 34:],
}  # km/s

B_norm = np.linalg.norm(B, axis=1)

v_plus = np.zeros((len(x), 3))
for ion in ["H+", "O+", "O2+", "CO2+"]:
    for i in range(3):
        v_plus[:, i] += densities[ion] * velocities[ion][:, i] / densities["e-"]

v_SI = v_plus * 1e3  # m/s
B_SI = B * 1e-9  # T
n_SI = densities["H+"] * 1e6  # 1/m3
J_SI = J * 1e-6  # A/m2
grad_p_SI = grad_p * 1e-9

Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)
Ep = np.array([-1 / (e_SI * n_SI[i]) * grad_p_SI[i, :] for i in range(len(grad_p))])

P_heavy = presion["O+"] + presion["O2+"] + presion["CO2+"]
P_B = np.linalg.norm(b1, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram = 1.67e-6 * densities["H+"] * velocities["H+"][:, 0] ** 2  # nPa
P_total = P_heavy + P_B + presion["e-"] + P_ram + presion["H+"]

rho_heavies = np.zeros(len(x))
for ion in ["O+", "O2+", "CO2+"]:
    rho_heavies += densities[ion]


"""Ancho MPB: hay una función al fondo que uso para elegir los límites que
tomé acá"""
lim_sup = 1.20
lim_inf = 1.175

x_cut = x[donde(x, lim_inf) : donde(x, lim_sup)]
B_cut = B_norm[donde(x, lim_inf) : donde(x, lim_sup)]

coef = np.polynomial.polynomial.polyfit(x_cut, B_cut, deg=1)

MPR = np.mean(B_norm[donde(x, lim_inf - 0.015) : donde(x, lim_inf)])

MPB_inicio = donde(coef[0] + x * coef[1], MPR)
MPB_fin = donde(x, lim_sup)

plt.figure()
plt.title("Campo magnético y ajuste")
plt.plot(x, B_norm, ".")
plt.plot(x, coef[0] + x * coef[1])
plt.axhline(y=MPR, c="k")
plt.axhline(y=B_norm[donde(x, lim_sup)], c="k")
plt.ylabel("|B|")
plt.xlabel("x (RM)")
plt.ylim(ymin=0)
plt.grid()
plt.show()

ancho_mpb = (x[MPB_fin] - x[MPB_inicio]) * 3390  # km
print(f"ancho de la mpb = {ancho_mpb:.3g} km")

# valores medios de los campos en la MPB
EH_medio = np.mean(np.linalg.norm(Ehall[MPB_inicio:MPB_fin], axis=1)) * 1e3
Ecv_medio = np.mean(np.linalg.norm(Ecv[MPB_inicio:MPB_fin], axis=1)) * 1e3
Ep_medio = np.mean(np.linalg.norm(Ep[MPB_inicio:MPB_fin], axis=1)) * 1e3
J_medio = np.mean(np.linalg.norm(J[MPB_inicio:MPB_fin], axis=1)) * 1e3
normal = [1, 0, 0]
Bdown = np.mean(B[donde(x, lim_inf - 0.015) : donde(x, lim_inf)], axis=0) * 1e-9
Bup = np.mean(B[donde(x, lim_sup) : donde(x, lim_sup + 0.015)], axis=0) * 1e-9
J_salto = 1 / (mu0 * ancho_mpb * 1e3) * np.cross(normal, Bup - Bdown) * 1e9

print(
    f"EHall medio = {EH_medio:.3g} mV/m\nEcv medio = {Ecv_medio:.3g} mV/m\nEp medio = {Ep_medio:.3g} mV/m"
)
print(f"J medio = {J_medio:.3g} nA/m², J jump = {np.linalg.norm(J_salto):.3g} nA/m²")

# longitud inercial de protones
paso = 20

density_mean = [
    np.mean(densities["H+"][i : i + paso]) for i in range(len(densities["H+"]) - paso)
]

ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km

print(
    f"la longitud inercial de iones en la zona Upstream es {np.mean(ion_length[donde(x, lim_sup) : donde(x, lim_sup + 0.015)]):.3g} km"
)
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
Giroradio térmico en la MPB
"""

v_th = np.sqrt(
    np.mean(presion["H+"][MPB_inicio:MPB_fin], axis=0)
    * 1e-21
    / (np.mean(densities["H+"][MPB_inicio:MPB_fin], axis=0) * mp)
)  # km/s

B_medio = np.mean(B[MPB_inicio:MPB_fin], axis=0)

rg = mp * v_th / (e_SI * np.linalg.norm(B_medio)) * 1e9  # km

print(f"El giroradio térmico medio en la MPB es {rg:1.3g} km")

# punto a punto
v_th = np.sqrt(presion["H+"] * 1e-21 / (densities["H+"] * mp))  # km/s

rg = [
    mp * v_th[i] / (e_SI * np.linalg.norm(B[i, :])) * 1e9 for i in range(len(B))
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
plt.rcParams.update({"font.size": 6})
plt.figure()
ax1 = plt.subplot2grid((6, 1), (0, 0))
ax2 = plt.subplot2grid((6, 1), (3, 0))
ax3 = plt.subplot2grid((6, 1), (2, 0))
ax4 = plt.subplot2grid((6, 1), (1, 0))
ax5 = plt.subplot2grid((6, 1), (4, 0))
ax6 = plt.subplot2grid((6, 1), (5, 0))

ax1.plot(x, B[:, 0], label="Bx")
ax1.plot(x, B[:, 1], label="By")
ax1.plot(x, B[:, 2], label="Bz")
ax1.plot(x, np.linalg.norm(B, axis=1), label="|B|")
ax1.set_ylabel("B (nT)")

ax2.plot(x, presion["H+"], label="H+")
ax2.plot(x, presion["e-"], label="e-")
ax2.plot(x, P_heavy, label="heavies")
ax2.plot(x, P_B, label="mag")
ax2.plot(x, P_ram, label="ram")
ax2.plot(x, P_total, label="total")
ax2.set_ylim([0, 1])
ax2.set_ylabel("Pressure (nPa)")
ax2.legend(loc="upper right")

ax3.plot(x, densities["H+"], label="H+")
ax3.plot(x, rho_heavies, label="heavies")
ax3.plot(x, densities["e-"], label="e-")
ax3.set_ylabel("Density (mp/cc)")
ax3.set_ylim([0, 50])

ax4.plot(x, J[:, 0], label="Jx")
ax4.plot(x, J[:, 1], label="Jy")
ax4.plot(x, J[:, 2], label="Jz")
ax4.plot(x, np.linalg.norm(J, axis=1), label="|J|")
ax4.set_ylabel("J (uA/m²)")
ax4.set_ylim([-0.1, 0.2])

ax5.plot(x, np.linalg.norm(Ecv, axis=1) * 1e3, label="Ecv")
ax5.plot(x, np.linalg.norm(Ehall, axis=1) * 1e3, label="Ehall")
ax5.plot(x, np.linalg.norm(Ep, axis=1) * 1e3, label="Ep")
ax5.set_ylabel("E (mV/m)")
ax5.set_ylim([0, 5])

# ax6.scatter(x, Ehall[:, 0] * 1e3, s=1, label="x")
# ax6.scatter(x, Ehall[:, 1] * 1e3, s=1, label="y")
# ax6.scatter(x, Ehall[:, 2] * 1e3, s=1, label="z")
ax6.plot(x, Ehall * 1e3)
ax6.set_ylabel("Ehall (mV/m)")
ax6.set_xlabel("x (RM)")
ax6.legend(["x", "y", "z"], loc="upper right")
ax6.set_xlim([1.15, 1.3])
ax6.set_ylim([-1, 4])
ax6.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2, label="MPB")
ax6.grid()

for ax in [ax1, ax2, ax3, ax4, ax5]:
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_xlim([1.15, 1.3])
    ax.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2)
    ax.legend(loc="upper right")
    ax.grid()

plt.show()


# E componentes
plt.figure()
ax4 = plt.subplot2grid((3, 1), (0, 0))
ax5 = plt.subplot2grid((3, 1), (1, 0))
ax6 = plt.subplot2grid((3, 1), (2, 0))

ax4.set_title("Electric fields")
ax4.plot(x, Ecv * 1e3, ".")
ax4.set_ylabel("Ecv (mV/m)")
plt.setp(ax4.get_xticklabels(), visible=False)

ax5.plot(x, Ehall * 1e3, ".")
ax5.set_ylabel("Ehall (mV/m)")
plt.setp(ax5.get_xticklabels(), visible=False)

ax6.plot(x, Ep * 1e3, ".")
ax6.set_ylabel("Ep (mV/m)")

for ax in [ax6, ax4, ax5]:
    ax.set_ylim([-5, 5])
    ax.set_xlim([1.1, 1.7])
    ax.grid()
    ax.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="k", alpha=0.2)

ax6.set_xlabel("x (RM)")
ax6.legend(["x", "y", "z", "MPB"])


plt.show()


def ancho(x, B):
    """
    Busco el punto máximo de la derivada y luego ajusto con una lineal el campo B
    pero que pase por el punto máximo
    """
    dBx_dx = np.gradient(B[:, 0], np.abs(x[0] - x[2]) * 3390e3)
    dBy_dx = np.gradient(B[:, 1], np.abs(x[0] - x[2]) * 3390e3)
    dBz_dx = np.gradient(B[:, 2], np.abs(x[0] - x[2]) * 3390e3)
    dBdx = np.transpose(np.vstack((dBx_dx, dBy_dx, dBz_dx)))
    B_norm = np.linalg.norm(B, axis=1)

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
