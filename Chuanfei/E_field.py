import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

from funciones import donde


"""
Vamos a escribir cada término de la ley de Ohm generalizada.
Necesito los valores de v, B, j y p en una región más amplia que el eje x nomás
"""
path = "../../../datos/simulacion_chuanfei/"
datos_on = np.loadtxt(path + "ejex_HallOn.gz")
datos_off = np.loadtxt(path + "ejex_HallOff.gz")

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
e_SI = 1.6e-19  # C
g = 3.7  # Mars surface gravity, m/s²

"""
r1 r2 rho Hrho Orho O2rho CO2rho H_vx H_vy H_vz
Bx By Bz b1x b1y b1z Pe P HP OP
O2P CO2P jx jy jz gradx grady gradz
"""

x = datos_on[:, 0]  # RM
B_on = datos_on[:, 10:13]  # nT
b1_on = datos_on[:, 13:16]  # nT
J_on = datos_on[:, 22:25]  # ua/m2
grad_p_on = datos_on[:, 25:28]  # nPa/m

HpRho_off = datos_off[:, 3]  # Mp/cc
B_off = datos_off[:, 10:13]  # nT
b1_off = datos_off[:, 13:16]  # nT
J_off = datos_off[:, 22:25]  # ua/m2
grad_p_off = datos_off[:, 25:28]  # nPa/m
v_i_off = datos_off[:, 7:10]  # km/s
v_O_off = datos_off[:, 28:31]  # km/s
v_O2_off = datos_off[:, 31:34]  # km/s
v_CO2_off = datos_off[:, 34:]  # km/s

# On

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


v_plus = np.zeros((568, 3))
for ion in ["H+", "O+", "O2+", "CO2+"]:
    for i in range(3):
        v_plus[:, i] += densities[ion] * velocities[ion][:, i] / densities["e-"]

v_SI_on = v_plus * 1e3  # m/s
B_SI_on = B_on * 1e-9  # T
n_SI_on = densities["H+"] * 1e6  # 1/m3
J_SI_on = J_on * 1e-6  # A/m2
grad_p_SI_on = grad_p_on * 1e-9

Ecv_on = np.array(
    [-np.cross(velocities["H+"][i] * 1e3, B_SI_on[i]) for i in range(len(B_on))]
)  # V/m
Ecv_on_plus = np.array(
    [-np.cross(v_SI_on[i], B_SI_on[i]) for i in range(len(B_on))]
)  # V/m
Ehall_on = np.array(
    [
        1 / (e_SI * n_SI_on[i]) * np.cross(J_SI_on[i], B_SI_on[i])
        for i in range(len(B_on))
    ]
)
Ep_on = np.array(
    [-1 / (e_SI * n_SI_on[i]) * grad_p_SI_on[i, :] for i in range(len(grad_p_on))]
)


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

v_total_off = (v_i_off + v_O_off + v_O2_off + v_CO2_off) * 1e3  # m/s

Ecv_off_total = np.array(
    [-np.cross(v_total_off[i], B_SI_off[i]) for i in range(len(B_off))]
)  # V/m


# J_mean = np.mean(J[inicio_MPB:fin_MPB], axis=0)
# B_mean = np.mean(B[inicio_MPB:fin_MPB], axis=0)
#
# print(
#     f"corriente simu j = {J_mean*1e3} nA/m², corriente MAVEN j = (-34.773,-233.373,-153.737) nA/m²"
# )
# print(f"campo simu B = {B_mean} nT, campo MAVEN B = (11.71,26.29,-19.55) nT")

plt.figure()
ax1 = plt.subplot2grid((3, 2), (0, 0))
ax2 = plt.subplot2grid((3, 2), (1, 0))
ax3 = plt.subplot2grid((3, 2), (2, 0))
ax4 = plt.subplot2grid((3, 2), (0, 1))
ax5 = plt.subplot2grid((3, 2), (1, 1))
ax6 = plt.subplot2grid((3, 2), (2, 1))

ax1.plot(x, Ecv_on * 1e3)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylim([-5, 5])
ax1.set_ylabel("E_cv (mV/m)")
ax1.set_title("Protones")

ax2.plot(x, B_on)
ax2.set_ylim([-50, 50])
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")

ax3.plot(x, velocities["H+"])
ax3.set_xlabel("x (RM)")
ax3.set_ylim([-500, 500])
ax3.set_ylabel("v_i (km/s)")

ax4.plot(x, Ecv_on_plus * 1e3)
ax4.set_ylim([-5, 5])
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.set_title("U+ (Dong et al., 2014)")

ax5.plot(x, B_on)
ax5.set_ylim([-50, 50])
plt.setp(ax5.get_xticklabels(), visible=False)

ax6.plot(x, v_plus)
ax6.set_ylim([-500, 500])
ax6.set_xlabel("x (RM)")
ax6.legend(["x", "y", "z"])

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.grid()
    ax.set_xlim([1.1, 1.9])

plt.show()

plt.figure()
ax1 = plt.subplot2grid((3, 1), (0, 0))
ax2 = plt.subplot2grid((3, 1), (1, 0))
ax3 = plt.subplot2grid((3, 1), (2, 0))
# plt.plot(x[inicio_MPB:fin_BS], np.linalg.norm(Ecv[inicio_MPB:fin_BS, :], axis=1) * 1e3)

ax1.plot(x, Ehall_on * 1e3)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("E_hall (mV/m)")
ax1.set_ylim([-5, 5])
ax1.set_title("Hall On")

ax2.plot(x, B_on)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("B (nT)")
ax2.set_ylim(ymax=51)

ax3.plot(x, J_on)
ax3.legend(["x", "y", "z"])
ax3.set_ylabel("J (uA/m²)")
ax3.set_ylim([-0.1, 0.2])
ax3.set_xlabel("x (RM)")

for ax in [ax1, ax2, ax3]:
    ax.grid()
    ax.set_xlim([1.1, 1.9])

plt.show()

plt.figure()
plt.plot(x, Ep_on * 1e3)
plt.grid()
plt.ylim([-1.5, 2.5])
plt.xlim([1.1, 1.9])
plt.xlabel("x (R_M)")
plt.ylabel("E_p (mV/m)")
plt.title("Hall On")

plt.figure()
plt.plot(x, np.linalg.norm(Ecv_on_plus, axis=1) * 1e3, ".")
plt.plot(x, np.linalg.norm(Ehall_on, axis=1) * 1e3, ".")
plt.plot(x, np.linalg.norm(Ep_on, axis=1) * 1e3, ".")
plt.legend(["|E cv|", "|E hall|", "|Ep|"])
plt.xlim([1.1, 1.9])
plt.ylim([0, 5])
plt.title("Hall On")
plt.xlabel("x (R_M)")
plt.ylabel("E (mV/m)")
plt.grid()
plt.show()
