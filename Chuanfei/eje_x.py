import numpy as np
import matplotlib.pyplot as plt

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

x_on = datos_on[:, 0]  # RM
HpRho_on = datos_on[:, 3]  # Mp/cc
v_i_on = datos_on[:, 7:10]  # km/s
B_on = datos_on[:, 10:13]  # nT
b1_on = datos_on[:, 13:16]  # nT
Pe_on = datos_on[:, 16]  # nPa
HP_on = datos_on[:, 18]  # nPa
OP_on = datos_on[:, 19]  # nPa
O2P_on = datos_on[:, 20]  # nPa
CO2P_on = datos_on[:, 21]  # nPa
J_on = datos_on[:, 22:25]  # ua/m2
grad_p_on = datos_on[:, 25:]  # nPa/m

x_off = datos_off[:, 0]  # RM
HpRho_off = datos_off[:, 3]  # Mp/cc
v_i_off = datos_off[:, 7:10]  # km/s
B_off = datos_off[:, 10:13]  # nT
b1_off = datos_off[:, 13:16]  # nT
Pe_off = datos_off[:, 16]  # nPa
HP_off = datos_off[:, 18]  # nPa
OP_off = datos_off[:, 19]  # nPa
O2P_off = datos_off[:, 20]  # nPa
CO2P_off = datos_off[:, 21]  # nPa
J_off = datos_off[:, 22:25]  # ua/m2
grad_p_off = datos_off[:, 25:]  # nPa/m

# On
v_SI_on = v_i_on * 1e3  # m/s
B_SI_on = B_on * 1e-9  # T
n_SI_on = HpRho_on * 1e6  # 1/m3
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

P_heavy_on = OP_on + O2P_on + CO2P_on
P_B_on = np.linalg.norm(B_on, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram_on = 1.67e-6 * HpRho_on * v_i_on[:, 0] ** 2  # nPa
P_total_on = P_heavy_on + P_B_on + Pe_on + P_ram_on + HP_on

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


# plots
plt.figure()
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))

ax1.scatter(x_off, Pe_off)
ax1.scatter(x_off, P_heavy_off)
ax1.scatter(x_off, P_B_off)
ax1.scatter(x_off, P_ram_off)
ax1.scatter(x_off, P_total_off)
ax1.scatter(x_off, HP_off)

ax2.scatter(x_on, Pe_on, label="Presion electronica")
ax2.scatter(x_on, P_heavy_on, label="Presion heavies")
ax2.scatter(x_on, P_B_on, label="Presion B")
ax2.scatter(x_on, P_ram_on, label="Presion ram")
ax2.scatter(x_on, P_total_on, label="presión total")
ax2.scatter(x_on, HP_on, label="Presion H+", marker=".")

plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_title("Hall Off")
ax2.set_title("Hall On")
ax1.set_ylabel("P (nPa)")
ax1.set_xlabel("x (RM)")
ax2.set_xlabel("x (RM)")

ax1.set_ylim([-0.1, 1.5])
ax2.set_ylim([-0.1, 1.5])

ax2.legend()
ax1.grid()
ax2.grid()
plt.show()


plt.figure()
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))

ax1.plot(x_off, J_off, ".")
ax1.scatter(x_off, np.linalg.norm(J_off, axis=1), c="k")

ax2.plot(x_on, J_on, ".")
ax2.scatter(x_on, np.linalg.norm(J_on, axis=1), c="k")

plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_title("Hall Off")
ax2.set_title("Hall On")
ax1.set_ylabel("corrientes (uA/m2)")
ax1.set_xlabel("x (RM)")
ax2.set_xlabel("x (RM)")
ax1.set_ylim([-0.2, 0.2])
ax2.set_ylim([-0.2, 0.2])
ax1.set_xlim([1, 1.8])
ax2.set_xlim([1, 1.8])

ax2.legend(["Jx", "Jy", "Jz", "|J|"])
ax1.grid()
ax2.grid()
plt.show()


plt.figure()
ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

ax1.plot(x_off, B_off, ".")
ax1.scatter(x_off, np.linalg.norm(B_off, axis=1), c="k")

ax2.plot(x_on, B_on, ".")
ax2.scatter(x_on, np.linalg.norm(B_on, axis=1), c="k")

ax3.plot(x_off, b1_off, ".")
ax3.scatter(x_off, np.linalg.norm(B_off, axis=1), c="k")

ax4.plot(x_on, b1_on, ".")
ax4.scatter(x_on, np.linalg.norm(B_on, axis=1), c="k")


plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_title("Hall Off")
ax2.set_title("Hall On")
ax1.set_ylabel("campo magnético (nT)")
ax3.set_ylabel("campo magnético \n sin corticales (nT)")
ax1.set_xlabel("x (RM)")
ax2.set_xlabel("x (RM)")
ax2.legend(["Bx", "By", "Bz", "|B|"])

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xlim([1, 1.8])
    ax.grid()

plt.show()


plt.figure()
ax1 = plt.subplot2grid((3, 2), (0, 0))
ax2 = plt.subplot2grid((3, 2), (1, 0))
ax3 = plt.subplot2grid((3, 2), (2, 0))
ax4 = plt.subplot2grid((3, 2), (0, 1))
ax5 = plt.subplot2grid((3, 2), (1, 1))
ax6 = plt.subplot2grid((3, 2), (2, 1))

ax1.plot(x_off, np.linalg.norm(Ecv_off, axis=1) * 1e3, ".")
ax1.set_ylabel("E_cv (mV/m)")
ax1.set_xlim([1.1, 1.8])
ax1.set_title("Hall Off")

ax2.plot(x_off, np.linalg.norm(Ehall_off, axis=1) * 1e3, ".")
ax2.set_ylabel("Ehall (mV/m)")
ax2.set_xlim([1.1, 1.8])
# ax2.set_ylim([-1, 6])

ax3.plot(x_off, np.linalg.norm(Ep_off, axis=1) * 1e3, ".")
ax3.set_xlim([1.1, 1.8])
ax3.set_ylim([0, 2])
ax3.legend(["x ", "y", "z"])
ax3.set_ylabel("Ep (mV/m)")

ax4.plot(x_on, np.linalg.norm(Ecv_on, axis=1) * 1e3, ".")
ax4.set_xlim([1.1, 1.8])
ax4.set_title("Hall On")

ax5.plot(x_on, np.linalg.norm(Ehall_on, axis=1) * 1e3, ".")
ax5.set_xlim([1.1, 1.8])
# ax5.set_ylim([-1, 6])

ax6.plot(x_on, np.linalg.norm(Ep_on, axis=1) * 1e3, ".")
ax6.set_xlim([1.1, 1.8])
ax6.set_ylim([0, 2])

for ax in [ax1, ax2, ax4, ax5]:
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_ylim([0, 2])
    ax.grid()

ax3.grid()
ax6.grid()
ax3.set_xlabel("x (RM)")
ax6.set_xlabel("x (RM)")


plt.show()

plt.figure()
ax4 = plt.subplot2grid((3, 1), (0, 0))
ax5 = plt.subplot2grid((3, 1), (1, 0))
ax6 = plt.subplot2grid((3, 1), (2, 0))

ax4.plot(x_on, Ecv_on[:, 0] * 1e3, ".")
ax4.set_xlim([1.1, 1.3])
ax4.set_title("Hall On")
ax4.set_ylim([-6, 1])

ax5.plot(x_on, Ehall_on[:, 0] * 1e3, ".")
ax5.set_xlim([1.1, 1.3])
ax5.set_ylim([-1, 6])

ax6.plot(x_on, Ep_on[:, 0] * 1e3, ".")
ax6.set_xlim([1.1, 1.3])
ax6.set_ylim([-1, 6])

plt.show()
