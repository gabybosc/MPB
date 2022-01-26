import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sys
import scipy.interpolate
from scipy.stats import multivariate_normal
import seaborn as sns

sns.set()

path = "../../../datos/simulacion_chuanfei/"
datos = np.loadtxt(path + "y=0_HallOn_new2.gz")

sys.path.append("..")
from funciones import donde, angulo


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
x = datos[:, 0]  # RM
y = datos[:, 1]  # RM
B = datos[:, 10:13]  # nT
b1 = datos[:, 13:16]  # nT
J = datos[:, 22:25]  # ua/m2
grad_p = datos[:, 25:28]  # nPa/m
presion = {
    "e": datos[:, 16],
    "H": datos[:, 18],
    "O": datos[:, 19],
    "O2": datos[:, 20],
    "CO2": datos[:, 21],
}  # nPa
densidad = {
    "e": datos[:, 2],
    "H": datos[:, 3],
    "O": datos[:, 4],
    "O2": datos[:, 5],
    "CO2": datos[:, 6],
}  # Mp/cc
velocidad = {
    "H": datos[:, 7:10],
    "O": datos[:, 28:31],
    "O2": datos[:, 31:34],
    "CO2": datos[:, 34:37],
    "e": datos[:, 37:],
}  # km/s

P_heavy = presion["O"] + presion["O2"] + presion["CO2"]
P_mag = np.linalg.norm(b1, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_dyn = 1.67e-6 * densidad["H"] * velocidad["H"][:, 0] ** 2  # nPa
P_total = P_heavy + P_mag + presion["e"] + P_dyn + presion["H"]

n_SI = densidad["H"] * 1e6  # 1/m3
delta_v = np.array([J[i, :] * 1e-6 / (e_SI * n_SI[i]) for i in range(len(J))])  # m/s

beta = (P_heavy + presion["H"]) / P_mag

beta_str = (P_heavy + presion["H"] + P_dyn) / P_mag

densidad_heavies = densidad["O"] + densidad["O2"] + densidad["CO2"]
density_ratio = densidad["H"] / densidad_heavies
mass_ratio = densidad["H"] / (
    densidad["O"] * 16 + densidad["O2"] * 32 + densidad["CO2"] * 44
)

v_plus = np.zeros((len(x), 3))
for ion in ["H", "O", "O2", "CO2"]:
    for i in range(3):
        v_plus[:, i] += densidad[ion] * velocidad[ion][:, i] / densidad["e"]

v_SI = v_plus * 1e3  # m/s
B_SI = B * 1e-9  # T
n_SI = densidad["H"] * 1e6  # 1/m3
J_SI = J * 1e-6  # A/m2
grad_p_SI = grad_p * 1e-9

Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)
Ep = np.array([-1 / (e_SI * n_SI[i]) * grad_p_SI[i, :] for i in range(len(grad_p))])

# MPB
x0 = 0.5  # 0.65 en y=0, 0.5 en z=0
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

# plt.scatter(x, y, c=beta_str, s=35, vmin=0, vmax=2, cmap="coolwarm")
# plt.plot(X1, Y1, c="r")
# plt.xlim([1, 1.5])
# plt.ylim([-1, 1])
# plt.colorbar()
# plt.show()


# Sample from 3D Gaussian distribution

idx = [i for i in range(len(beta_str)) if np.abs(beta_str[i] - 1) < 0.05]

adentro = [i for i in range(len(x)) if x[i] < 1 and y[i] < 1]

xy = np.column_stack([x.flat, y.flat])  # Create a (N, 2) array of (x, y) pairs.

# Interpolate and generate heatmap:
grid_x, grid_y = np.mgrid[x.min() : x.max() : 1000j, y.min() : y.max() : 1000j]

z = np.linalg.norm(J, axis=1) * 1e3  # nA/m²
plt.figure()
grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method="cubic")
# [pcolormesh with missing values?](https://stackoverflow.com/a/31687006/395857)
plt.pcolormesh(
    grid_x, grid_y, ma.masked_invalid(grid_z), cmap="inferno", vmin=0, vmax=100
)
plt.plot(X1, Y1, c="C0")
# plt.title("{0} interpolation".format(method))
plt.title("|J| in z=0")
plt.xlabel("x (RM)")
plt.ylabel("y (RM)")
plt.colorbar()
plt.xlim([1, 2])
plt.ylim([-1, 1])
plt.show()

z = beta
plt.figure()
grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method="cubic")
# [pcolormesh with missing values?](https://stackoverflow.com/a/31687006/395857)
plt.pcolormesh(
    grid_x, grid_y, ma.masked_invalid(grid_z), cmap="coolwarm", vmin=0, vmax=2
)
plt.plot(X1, Y1, c="C0")
# plt.title("{0} interpolation".format(method))
plt.title("beta in z=0")
plt.xlabel("x (RM)")
plt.ylabel("y (RM)")
plt.colorbar()
plt.xlim([1, 2])
plt.ylim([-1, 1])
plt.show()

z = np.linalg.norm(B, axis=1)  # nT
plt.figure()
grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method="cubic")
# [pcolormesh with missing values?](https://stackoverflow.com/a/31687006/395857)
plt.pcolormesh(
    grid_x, grid_y, ma.masked_invalid(grid_z), cmap="inferno", vmin=0, vmax=60
)
plt.plot(X1, Y1, c="C0")
# plt.title("{0} interpolation".format(method))
plt.title("|B| in y=0")
plt.xlabel("x (RM)")
plt.ylabel("z (RM)")
plt.colorbar()
plt.xlim([1, 2])
plt.ylim([-1, 1])
plt.show()

fig, axs = plt.subplots(1, 3)
for i in range(3):
    z = B[:, i]
    grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method="cubic")
    sc = axs[i].pcolormesh(
        grid_x, grid_y, ma.masked_invalid(grid_z), cmap="inferno", vmin=0, vmax=60
    )
    axs[i].plot(X1, Y1, c="C2")
    axs[i].set_xlim([1, 2])
    axs[i].set_ylim([-1, 1])
    axs[i].set_title(f"Y=0, B {i}")
    axs[i].set_xlabel("X MSO (RM)")
    if i > 0:
        plt.setp(axs[i].get_yticklabels(), visible=False)
    axs[0].set_ylabel("Z MSO (RM)")
plt.colorbar(sc)
plt.show()

# z = densidad["H"]
# plt.figure()
# grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method="cubic")
# # [pcolormesh with missing values?](https://stackoverflow.com/a/31687006/395857)
# plt.pcolormesh(
#     grid_x, grid_y, ma.masked_invalid(grid_z), cmap="inferno", vmin=0, vmax=10
# )
# plt.plot(X1, Y1, c="C0")
# # plt.title("{0} interpolation".format(method))
# plt.title("Dependence of H+ density with SZA")
# plt.xlabel("x (RM)")
# plt.ylabel("z (RM)")
# plt.colorbar()
# plt.xlim([1, 2])
# plt.ylim([-1, 1])
# plt.show()
#
# fig, axs = plt.subplots(2, 3)
# i = 0
# for den in densidad:
#     print(den)
#     z = densidad[den]
#     grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method="cubic")
#     if i > 2:
#         axs[1, i - 3].pcolormesh(
#             grid_x, grid_y, ma.masked_invalid(grid_z), cmap="inferno", vmin=0, vmax=10
#         )
#         axs[1, i - 3].plot(X1, Y1, c="C2")
#         axs[1, i - 3].set_xlim([1, 2])
#         axs[1, i - 3].set_ylim([-1, 1])
#         axs[1, i - 3].set_title(f"Y=0, densidad de {den}")
#         if i - 3 > 0:
#             plt.setp(axs[1, i - 3].get_yticklabels(), visible=False)
#
#     else:
#         axs[0, i].pcolormesh(
#             grid_x, grid_y, ma.masked_invalid(grid_z), cmap="inferno", vmin=0, vmax=10
#         )
#         axs[0, i].plot(X1, Y1, c="C2")
#         axs[0, i].set_xlim([1, 2])
#         axs[0, i].set_ylim([-1, 1])
#         axs[0, i].set_title(f"Y=0, densidad de {den}")
#         plt.setp(axs[0, i].get_xticklabels(), visible=False)
#         if i > 0:
#             plt.setp(axs[0, i].get_yticklabels(), visible=False)
#     axs[0, 0].set_ylabel("Z MSO (RM)")
#     axs[1, 0].set_ylabel("Z MSO (RM)")
#     i += 1
# plt.show()

# sc2 = axs[1, i].scatter(
#     x, y, c=J[:, i], vmin=-0.25, vmax=0.25, s=35, cmap="coolwarm",
# )
# axs[1, i].plot(X1, Y1, c="C2")
# axs[1, i].set_xlabel("X MSE (RM)")
# axs[1, i].set_title(f"Z=0, JB{coordenadas[i]}")
# axs[1, i].set_xlim([0, 2])
# axs[1, i].set_ylim([-1, 1])
# if i > 0:
#     plt.setp(axs[1, i].get_yticklabels(), visible=False)
# else:
#     axs[1, i].set_ylabel("Y MSE (RM)")
#
# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.55, 0.05, 0.38])  # [left, bottom, width, height]
# fig.colorbar(sc, cax=cbar_ax)
#
# fig.subplots_adjust(right=0.8)
# cbar_ax2 = fig.add_axes([0.85, 0.05, 0.05, 0.38])
# fig.colorbar(sc2, cax=cbar_ax2)
