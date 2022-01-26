import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.interpolate
from scipy.stats import multivariate_normal
import seaborn as sns
import numpy.ma as ma
import matplotlib.cm as cm
from matplotlib.colors import Normalize

sns.set()

path = "../../../datos/simulacion_chuanfei/"
datos_y = np.loadtxt(path + "y=0_HallOn_new2.gz")
datos_z = np.loadtxt(path + "z=0_HallOn_new2.gz")

sys.path.append("..")
from funciones import donde, angulo
from funciones_plot import plot_2d

"""
r1 r2 rho Hrho Orho O2rho CO2rho H_vx H_vy H_vz
Bx By Bz b1x b1y b1z Pe P HP OP
O2P CO2P jx jy jz gradx grady gradz O_vx O_vy
O_vz O2_vx O2_vy O2_vz CO2_vx CO2_vy CO2_vz e_vx e_vy e_vz
"""

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
g = 3.7  # Mars surface gravity, m/s²
e_SI = 1.6e-19  # C

# vectores
x = datos_y[:, 0]  # RM
y = datos_z[:, 1]  # RM
z = datos_y[:, 1]  # RM

B_y = datos_y[:, 10:13]  # nT
b1_y = datos_y[:, 13:16]  # nT
J_y = datos_y[:, 22:25]  # ua/m2

B_z = datos_z[:, 10:13]  # nT
b1_z = datos_z[:, 13:16]  # nT
J_z = datos_z[:, 22:25]  # ua/m2

presion_y = {
    "e": datos_y[:, 16],
    "H": datos_y[:, 18],
    "O": datos_y[:, 19],
    "O2": datos_y[:, 20],
    "CO2": datos_y[:, 21],
}  # nPa
densidad_y = {
    "e": datos_y[:, 2],
    "H": datos_y[:, 3],
    "O": datos_y[:, 4],
    "O2": datos_y[:, 5],
    "CO2": datos_y[:, 6],
}  # Mp/cc
velocidad_y = {
    "H": datos_y[:, 7:10],
    "O": datos_y[:, 28:31],
    "O2": datos_y[:, 31:34],
    "CO2": datos_y[:, 34:37],
    "e": datos_y[:, 37:],
}  # km/s

presion_z = {
    "e": datos_z[:, 16],
    "H": datos_z[:, 18],
    "O": datos_z[:, 19],
    "O2": datos_z[:, 20],
    "CO2": datos_z[:, 21],
}  # nPa
densidad_z = {
    "e": datos_z[:, 2],
    "H": datos_z[:, 3],
    "O": datos_z[:, 4],
    "O2": datos_z[:, 5],
    "CO2": datos_z[:, 6],
}  # Mp/cc
velocidad_z = {
    "H": datos_z[:, 7:10],
    "O": datos_z[:, 28:31],
    "O2": datos_z[:, 31:34],
    "CO2": datos_z[:, 34:37],
    "e": datos_z[:, 37:],
}  # km/s

P_heavy_z = presion_z["O"] + presion_z["O2"] + presion_z["CO2"]
P_mag_z = np.linalg.norm(b1_z, axis=1) ** 2 * 1e-9 / (2 * mu0)  # sin corticales
P_mag_cort_z = np.linalg.norm(B_z, axis=1) ** 2 * 1e-9 / (2 * mu0)  # con corticales
P_dyn_z = 1.67e-6 * densidad_z["H"] * velocidad_z["H"][:, 0] ** 2  # nPa
P_total_z = P_heavy_z + P_mag_z + presion_z["e"] + P_dyn_z + presion_z["H"]

beta_z = (P_heavy_z + presion_z["H"]) / P_mag_z
beta_str_z = (P_heavy_z + presion_z["H"] + P_dyn_z) / P_mag_z

P_heavy_y = presion_y["O"] + presion_y["O2"] + presion_y["CO2"]
P_mag_y = np.linalg.norm(b1_y, axis=1) ** 2 * 1e-9 / (2 * mu0)  # sin corticales
P_mag_cort_y = np.linalg.norm(B_y, axis=1) ** 2 * 1e-9 / (2 * mu0)  # con corticales
P_dyn_y = 1.67e-6 * densidad_y["H"] * velocidad_y["H"][:, 0] ** 2  # nPa
P_total_y = P_heavy_y + P_mag_y + presion_y["e"] + P_dyn_y + presion_y["H"]

beta_y = (P_heavy_y + presion_y["H"]) / P_mag_y
beta_str_y = (P_heavy_y + presion_y["H"] + P_dyn_y) / P_mag_y

# beta_cort = (P_heavy + presion_z["H"]) / P_mag_cort
# beta_str_cort = (P_heavy + presion_z["H"] + P_dyn) / P_mag_cort

densidad_z_heavies = densidad_z["O"] + densidad_z["O2"] + densidad_z["CO2"]
density_ratio = densidad_z["H"] / densidad_z_heavies
mass_ratio = densidad_z["H"] / (
    densidad_z["O"] * 16 + densidad_z["O2"] * 32 + densidad_z["CO2"] * 44
)

v_plus = np.zeros((len(x), 3))
for ion in ["H", "O", "O2", "CO2"]:
    for i in range(3):
        v_plus[:, i] += densidad_z[ion] * velocidad_z[ion][:, i] / densidad_z["e"]

# v_SI = v_plus * 1e3  # m/s
# B_SI = B_z * 1e-9  # T
# n_SI = densidad_z["H"] * 1e6  # 1/m3
# J_SI = J_z * 1e-6  # A/m2
# grad_p_SI = grad_p * 1e-9
# Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
# Ehall = np.array(
#     [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
# )
# Ep = np.array([-1 / (e_SI * n_SI[i]) * grad_p_SI[i, :] for i in range(len(grad_p))])

# MPB
x0 = 0.5
e = 0.9

# R = [1.082, -0.064, 0.515]
R = [1.2, 0, 0]  # fit de la MPB simulada
theta = np.linspace(0, np.pi * 2, 100)

r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(theta))

X1 = x0 + r1 * np.cos(theta)
Y1 = r1 * np.sin(theta)

plt.figure()
plot_2d(x, y, beta_z, 0, 2, "coolwarm")
plt.plot(X1, Y1, c="k")
plt.title("beta in z=0")
plt.xlabel("x (RM)")
plt.ylabel("y (RM)")
plt.xlim([1, 1.6])
plt.ylim([-1, 1])
plt.show()


def subplot_2d(
    x, y, z, zmin, zmax, ax, i, j, titulo, colormap="inferno", method="linear"
):
    xy = np.column_stack([x.flat, y.flat])  # Create a (N, 2) array of (x, y) pairs.
    # Interpolate and generate heatmap:
    grid_x, grid_y = np.mgrid[x.min() : x.max() : 1000j, y.min() : y.max() : 1000j]

    # interpolation method can be linear or cubic
    grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method=method)

    ax[i, j].pcolormesh(
        grid_x, grid_y, ma.masked_invalid(grid_z), cmap=colormap, vmin=zmin, vmax=zmax
    )
    ax[i, j].plot(X1, Y1, c="k", linestyle="--")
    ax[i, j].set_xlim([1, 1.6])
    ax[i, j].set_ylim([-1, 1])
    ax[i, j].set_title(titulo)


# betas
fig, ax = plt.subplots(2, 2)
subplot_2d(x, y, beta_z, 0, 2, ax, 0, 0, "beta z=0", "coolwarm")
subplot_2d(x, y, beta_str_z, 0, 2, ax, 0, 1, "beta* z=0", "coolwarm")
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1]:
        ax[0, k].quiver(x[i], y[i], B_z[i, 0], B_z[i, 1], color="k", alpha=0.5)
subplot_2d(x, z, beta_y, 0, 2, ax, 1, 0, "beta y=0", "coolwarm")
subplot_2d(x, z, beta_str_y, 0, 2, ax, 1, 1, "beta* y=0", "coolwarm")
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1]:
        ax[1, k].quiver(x[i], z[i], B_y[i, 0], B_y[i, 2], color="k", alpha=0.5)
for i in [0, 1]:
    plt.setp(ax[0, i].get_xticklabels(), visible=False)
for i in [0, 1]:
    plt.setp(ax[i, 1].get_yticklabels(), visible=False)
cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
fig.colorbar(cm.ScalarMappable(norm=Normalize(0, 2), cmap="coolwarm"), cax=cbar_ax)
ax[0, 0].set_ylabel("y (RM)")
ax[1, 0].set_ylabel("z (RM)")
ax[1, 0].set_xlabel("x (RM)")
ax[1, 1].set_xlabel("x (RM)")
plt.show()


# Corrientes
fig, ax = plt.subplots(2, 3)
subplot_2d(x, z, J_y[:, 0] * 1e3, -100, 100, ax, 0, 0, "Jx y=0", "coolwarm")
subplot_2d(x, z, J_y[:, 1] * 1e3, -100, 100, ax, 0, 1, "Jy y=0", "coolwarm")
subplot_2d(x, z, J_y[:, 2] * 1e3, -100, 100, ax, 0, 2, "Jz y=0", "coolwarm")
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1, 2]:
        ax[0, k].quiver(x[i], z[i], J_y[i, 0], J_y[i, 2], color="k", alpha=0.5)
subplot_2d(x, y, J_z[:, 0] * 1e3, -100, 100, ax, 1, 0, "Jx z=0", "coolwarm")
subplot_2d(x, y, J_z[:, 1] * 1e3, -100, 100, ax, 1, 1, "Jy z=0", "coolwarm")
subplot_2d(x, y, J_z[:, 2] * 1e3, -100, 100, ax, 1, 2, "Jz z=0", "coolwarm")
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1, 2]:
        ax[1, k].quiver(x[i], y[i], J_z[i, 0], J_z[i, 1], color="k", alpha=0.5)
for i in [0, 1, 2]:
    plt.setp(ax[0, i].get_xticklabels(), visible=False)
    ax[1, i].set_xlabel("x (RM)")
for i in [0, 1]:
    plt.setp(ax[i, 1].get_yticklabels(), visible=False)
    plt.setp(ax[i, 2].get_yticklabels(), visible=False)
ax[0, 0].set_ylabel("z (RM)")
ax[1, 0].set_ylabel("y (RM)")

cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.85])  # [left, bottom, width, height]
fig.colorbar(cm.ScalarMappable(norm=Normalize(-100, 100), cmap="coolwarm"), cax=cbar_ax)
plt.show()


# densidad
fig, ax = plt.subplots(2, 2)
subplot_2d(x, y, densidad_z["H"], 0, 20, ax, 0, 0, "H z=0", "inferno")
subplot_2d(x, y, densidad_z["e"], 0, 20, ax, 0, 1, "e- z=0", "inferno")
subplot_2d(x, z, densidad_y["H"], 0, 20, ax, 1, 0, "H y=0", "inferno")
subplot_2d(x, z, densidad_y["e"], 0, 20, ax, 1, 1, "e- y=0", "inferno")
for i in np.linspace(0, 12800, 258):
    i = int(i)
    ax[0, 0].quiver(
        x[i], y[i], velocidad_z["H"][i, 0], velocidad_z["H"][i, 1], color="k", alpha=0.5
    )
for i in np.linspace(0, 12800, 258):
    i = int(i)
    ax[0, 1].quiver(
        x[i], y[i], velocidad_z["e"][i, 0], velocidad_z["e"][i, 1], color="k", alpha=0.5
    )
for i in np.linspace(0, 12800, 258):
    i = int(i)
    ax[1, 0].quiver(
        x[i], z[i], velocidad_y["H"][i, 0], velocidad_y["H"][i, 2], color="k", alpha=0.5
    )
for i in np.linspace(0, 12800, 258):
    i = int(i)
    ax[1, 1].quiver(
        x[i], z[i], velocidad_y["e"][i, 0], velocidad_y["e"][i, 2], color="k", alpha=0.5
    )
for i in [0, 1]:
    plt.setp(ax[0, i].get_xticklabels(), visible=False)
for i in [0, 1]:
    plt.setp(ax[i, 1].get_yticklabels(), visible=False)
ax[0, 0].set_ylabel("y (RM)")
ax[1, 0].set_ylabel("z (RM)")
ax[1, 0].set_xlabel("x (RM)")
ax[1, 1].set_xlabel("x (RM)")
cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
fig.colorbar(cm.ScalarMappable(norm=Normalize(0, 20), cmap="inferno"), cax=cbar_ax)
plt.show()

# Normas
# fig, ax = plt.subplots(2, 2)
# subplot_2d(x, y, np.linalg.norm(B, axis=1), 0, 50, ax, 1, 0, "|B|")
# subplot_2d(x, y, np.linalg.norm(J * 1e3, axis=1), 0, 100, ax, 0, 1, "|J|")
# subplot_2d(x, y, np.linalg.norm(Ecv * 1e3, axis=1), 0, 3, ax, 1, 1, "|Ecv|")
# subplot_2d(x, y, np.linalg.norm(Ehall * 1e3, axis=1), 0, 5, ax, 0, 0, "|Ehall|")
# for i in [0, 1]:
#     plt.setp(ax[0, i].get_xticklabels(), visible=False)
# for i in [0, 1]:
#     plt.setp(ax[i, 1].get_yticklabels(), visible=False)
#     # plt.setp(ax[i, 2].get_yticklabels(), visible=False)
# plt.show()
# cbar_B = fig.add_axes([0.9, 0.55, 0.04, 0.35])  # [left, bottom, width, height]
# fig.colorbar(cm.ScalarMappable(norm=Normalize(-5, 5), cmap="coolwarm"), cax=cbar_ax)
# cbar_J = fig.add_axes([0.9, 0.1, 0.04, 0.35])  # [left, bottom, width, height]
# fig.colorbar(cm.ScalarMappable(norm=Normalize(-2.5, 2.5), cmap="coolwarm"), cax=cbar_ax)
# cbar_Ecv = fig.add_axes([0.9, 0.55, 0.04, 0.35])  # [left, bottom, width, height]
# fig.colorbar(cm.ScalarMappable(norm=Normalize(-5, 5), cmap="coolwarm"), cax=cbar_ax)
# cbar_Eh = fig.add_axes([0.9, 0.1, 0.04, 0.35])  # [left, bottom, width, height]
# fig.colorbar(cm.ScalarMappable(norm=Normalize(-2.5, 2.5), cmap="coolwarm"), cax=cbar_ax)
# fig.colorbar(cm.ScalarMappable(norm=Normalize(-2.5, 2.5), cmap="coolwarm"), cax=cbar_ax)
#
#
# fig, ax = plt.subplots(2, 3)
# plt.figure()
# plot_2d(x, y, beta, 0, 2)
# # plt.scatter(x, y, s=J[:, 0] * 1e3)
# plt.plot(X1, Y1, c="k")
# for i in np.linspace(0, 12800, 258):
#     i = int(i)
#     plt.quiver(x[i], y[i], B[i, 0], B[i, 2], color="r")
# # plt.quiver(1.2, 0.3, -0.0005, 0.05)
# plt.title("beta in y=0")
# plt.xlabel("x (RM)")
# plt.ylabel("z (RM)")
# plt.xlim([1, 1.6])
# plt.ylim([-1, 1])
# plt.show()
#     plt.quiver(x[i], y[i], B[i, 0], B[i, 2], color="r")
# # plt.quiver(1.2, 0.3, -0.0005, 0.05)
# plt.title("beta in y=0")
# plt.xlabel("x (RM)")
# plt.ylabel("z (RM)")
# plt.xlim([1, 1.6])
# plt.ylim([-1, 1])
# plt.show()
