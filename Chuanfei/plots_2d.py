import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.interpolate
from scipy.stats import multivariate_normal
import seaborn as sns
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import make_axes_locatable, ImageGrid, AxesGrid
import matplotlib.cm as cm
from matplotlib.colors import Normalize

sns.set()

path = "../../../datos/simulacion_chuanfei/nueva_simu/"
datos_y0 = np.loadtxt(path + "y=0.gz")
datos_z0 = np.loadtxt(path + "z=0.gz")

# path = "../../../datos/simulacion_chuanfei/"
# datos_y0 = np.loadtxt(path + "y=0_HallOn_new2.gz")
# datos_z0 = np.loadtxt(path + "z=0_HallOn_new2.gz")

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
x = datos_y0[:, 0]  # RM
y = datos_z0[:, 1]  # RM
z = datos_y0[:, 1]  # RM

B_y0 = datos_y0[:, 10:13]  # nT
b1_y0 = datos_y0[:, 13:16]  # nT
J_y0 = datos_y0[:, 22:25]  # ua/m2

B_z0 = datos_z0[:, 10:13]  # nT
b1_z0 = datos_z0[:, 13:16]  # nT
J_z0 = datos_z0[:, 22:25]  # ua/m2

presion_y0 = {
    "e": datos_y0[:, 16],
    "H": datos_y0[:, 18],
    "O": datos_y0[:, 19],
    "O2": datos_y0[:, 20],
    "CO2": datos_y0[:, 21],
}  # nPa
densidad_y0 = {
    "e": datos_y0[:, 2],
    "H": datos_y0[:, 3],
    "O": datos_y0[:, 4],
    "O2": datos_y0[:, 5],
    "CO2": datos_y0[:, 6],
}  # Mp/cc
velocidad_y0 = {
    "H": datos_y0[:, 7:10],
    "O": datos_y0[:, 28:31],
    "O2": datos_y0[:, 31:34],
    "CO2": datos_y0[:, 34:37],
    "e": datos_y0[:, 37:],
}  # km/s

presion_z0 = {
    "e": datos_z0[:, 16],
    "H": datos_z0[:, 18],
    "O": datos_z0[:, 19],
    "O2": datos_z0[:, 20],
    "CO2": datos_z0[:, 21],
}  # nPa
densidad_z0 = {
    "e": datos_z0[:, 2],
    "H": datos_z0[:, 3],
    "O": datos_z0[:, 4],
    "O2": datos_z0[:, 5],
    "CO2": datos_z0[:, 6],
}  # Mp/cc
velocidad_z0 = {
    "H": datos_z0[:, 7:10],
    "O": datos_z0[:, 28:31],
    "O2": datos_z0[:, 31:34],
    "CO2": datos_z0[:, 34:37],
    "e": datos_z0[:, 37:],
}  # km/s

presion_z0.update(
    {
        "heavies": presion_z0["O"] + presion_z0["O2"] + presion_z0["CO2"],
        "mag_b1": np.linalg.norm(b1_z0, axis=1) ** 2 * 1e-9 / (2 * mu0),
        "mag": np.linalg.norm(B_z0, axis=1) ** 2 * 1e-9 / (2 * mu0),
        "dyn": 1.67e-6 * densidad_z0["H"] * velocidad_z0["H"][:, 0] ** 2,
    }
)

presion_z0["total"] = (
    presion_z0["heavies"]
    + presion_z0["mag_b1"]
    + presion_z0["e"]
    + presion_z0["dyn"]
    + presion_z0["H"]
)

densidad_z0["heavies"] = densidad_z0["O"] + densidad_z0["O2"] + densidad_z0["CO2"]

beta_z0 = (presion_z0["heavies"] + presion_z0["H"]) / presion_z0["mag_b1"]
beta_str_z0 = (
    presion_z0["heavies"] + presion_z0["H"] + presion_z0["dyn"]
) / presion_z0["mag_b1"]

presion_y0.update(
    {
        "heavies": presion_y0["O"] + presion_y0["O2"] + presion_y0["CO2"],
        "mag_b1": np.linalg.norm(b1_y0, axis=1) ** 2 * 1e-9 / (2 * mu0),
        "mag": np.linalg.norm(B_y0, axis=1) ** 2 * 1e-9 / (2 * mu0),
        "dyn": 1.67e-6 * densidad_y0["H"] * velocidad_y0["H"][:, 0] ** 2,
    }
)

presion_y0["total"] = (
    presion_y0["heavies"]
    + presion_y0["mag_b1"]
    + presion_y0["e"]
    + presion_y0["dyn"]
    + presion_y0["H"]
)

densidad_y0["heavies"] = densidad_y0["O"] + densidad_y0["O2"] + densidad_y0["CO2"]

beta_y0 = (presion_y0["heavies"] + presion_y0["H"]) / presion_y0["mag_b1"]
beta_str_y0 = (
    presion_y0["heavies"] + presion_y0["H"] + presion_y0["dyn"]
) / presion_y0["mag_b1"]


densidad_z0["heavies"] = densidad_z0["O"] + densidad_z0["O2"] + densidad_z0["CO2"]
density_ratio = densidad_z0["H"] / densidad_z0["heavies"]
mass_ratio = densidad_z0["H"] / (
    densidad_z0["O"] * 16 + densidad_z0["O2"] * 32 + densidad_z0["CO2"] * 44
)

v_plus = np.zeros((len(x), 3))
for ion in ["H", "O", "O2", "CO2"]:
    for i in range(3):
        v_plus[:, i] += densidad_z0[ion] * velocidad_z0[ion][:, i] / densidad_z0["e"]

# xx = x[donde(x, 1) : donde(x, 1.3)]  # recorto para los x entre 1 y 1.3
# zz = z[donde(x, 1) : donde(x, 1.3)]
# Y = [j for j in yy if np.abs(j) < 0.5]  # recorto donde y está entre -0.5 y 0.5
# X = [xx[j] for j in range(len(yy)) if np.abs(yy[j]) < 0.5]
# Z = [zz[j] for j in range(len(yy)) if np.abs(yy[j]) < 0.5]


def recortar(x, y, arr):
    a = arr[donde(x, 1.1) : donde(x, 1.3)]
    s = np.array([a[j] for j in range(len(y)) if np.abs(y[j]) < 0.1])
    return s


yy = y[donde(x, 1.1) : donde(x, 1.3)]
zz = z[donde(x, 1.1) : donde(x, 1.3)]  # yy y zz son iguales
X = recortar(x, yy, x)
Y = recortar(x, zz, y)
Z = recortar(x, yy, z)
vz0H = recortar(x, yy, velocidad_z0["H"])
vy0H = recortar(x, zz, velocidad_y0["H"])
vz0e = recortar(x, yy, velocidad_z0["e"])
vy0e = recortar(x, zz, velocidad_y0["e"])
By0 = recortar(x, zz, B_y0)
Bz0 = recortar(x, yy, B_z0)
Jy0 = recortar(x, zz, J_y0)
Jz0 = recortar(x, yy, J_z0)
rhoH_y0 = recortar(x, zz, densidad_y0["H"])
rhoe_y0 = recortar(x, zz, densidad_y0["e"])
rhoH_z0 = recortar(x, yy, densidad_z0["H"])
rhoe_z0 = recortar(x, yy, densidad_z0["e"])


dif_vel_z0 = np.array(
    [rhoe_z0[i] * vz0e[i, :] - rhoH_z0[i] * vz0H[i, :] for i in range(len(vz0e))]
)
dif_vel_y0 = np.array(
    [rhoe_y0[i] * vy0e[i, :] - rhoH_y0[i] * vy0H[i, :] for i in range(len(vy0e))]
)

# Ehall_z0 = np.array(
#     [
#         1 / (e_SI * rhoH_z0[i] * 1e6) * np.cross(Jz[i] * 1e-6, Bz[i] * 1e-9) * 1e-3
#         for i in range(len(Bz))
#     ]
# )
#
# plt.plot(X, Ehall_z0)
# plt.show()
# v_SI = v_plus * 1e3  # m/s
# B_SI = B_z0 * 1e-9  # T
# n_SI = densidad_z0["H"] * 1e6  # 1/m3
# J_SI = J_z0 * 1e-6  # A/m2
# grad_p_SI = grad_p * 1e-9
# Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
# Ehall = np.array(
#     [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
# )
# Ep = np.array([-1 / (e_SI * n_SI[i]) * grad_p_SI[i, :] for i in range(len(grad_p))])

# beta = 1
x0 = 0.5

e = 0.9

R = [1.19, 0, 0]  # la standoff distance en x que tiene beta
theta = np.linspace(0, np.pi * 2, 100)

r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(theta))

X1 = x0 + r1 * np.cos(theta)
Y1 = r1 * np.sin(theta)


# MAVEN
x0 = 0.78
e = 0.9

L0 = 0.87  # el L obtenido en mi tesis de lic
r1 = L0 / (1 + e * np.cos(theta))

X1_M = x0 + r1 * np.cos(theta)
Y1_M = r1 * np.sin(theta)


def data(i, xy, z, grid_x, grid_y, zmin, zmax, metodo="linear", colormap="coolwarm"):
    grid_z = scipy.interpolate.griddata(xy, z, (grid_x, grid_y), method=metodo)
    grid[i].pcolormesh(
        grid_x, grid_y, ma.masked_invalid(grid_z), cmap=colormap, vmin=zmin, vmax=zmax
    )
    grid[i].set_aspect("equal", "box")
    grid[i].set_xlim([2, 1])
    grid[i].set_ylim([-0.5, 0.5])
    grid[i].plot(X1, Y1, c="k", linestyle="--", label="beta=1")
    grid[i].plot(X1_M, Y1_M, c="k", linestyle="-", label="MAVEN")


xy = np.column_stack([x.flat, y.flat])  # Create a (N, 2) array of (x, y) pairs.
grid_x, grid_y = np.mgrid[x.min() : x.max() : 1000j, y.min() : y.max() : 1000j]

# betas  # ojo que hay que escribir a mano el eje z(RM)
figure = plt.figure(1, (1.0, 2.0))
grid = AxesGrid(
    figure,
    111,  # similar to subplot(142)
    nrows_ncols=(1, 2),
    axes_pad=0.5,
    share_all=True,
    label_mode="L",
    cbar_location="right",
    cbar_mode="single",
)

data(0, xy, np.log(beta_str_z0), grid_x, grid_y, -4, 4)
data(1, xy, np.log(beta_str_y0), grid_x, grid_y, -4, 4)

grid[0].set_title(r"Z=0 log($\beta*$)")
grid[1].set_title(r"Y=0 log($\beta*$)")
grid[0].set_ylabel(r"y (R$_M$)")
grid[1].set_ylabel(r"z (R$_M$)")
grid[0].set_xlabel(r"x (R$_M$)")
grid[1].set_xlabel(r"x (R$_M$)")

figure.colorbar(
    cm.ScalarMappable(norm=Normalize(-4, 4), cmap="coolwarm"), cax=grid.cbar_axes[0]
)

figure.set_size_inches(9, 6)
# plt.savefig("../../../Dropbox/Paper2/beta_2d.png", dpi=600)
plt.show()


# Campo B

fig = plt.figure(1, (2.0, 3.0))
grid = AxesGrid(
    fig,
    111,  # similar to subplot(142)
    nrows_ncols=(2, 3),
    axes_pad=0.22,
    share_all=True,
    label_mode="L",
    cbar_location="right",
    cbar_mode="single",
)

data(0, xy, B_z0[:, 0], grid_x, grid_y, -50, 50)  # r"Z=0 B$_x$", "coolwarm")
data(1, xy, B_z0[:, 1], grid_x, grid_y, -50, 50)  # r"Z=0 B$_y$", "coolwarm")
data(2, xy, B_z0[:, 2], grid_x, grid_y, -50, 50)  # r"Z=0 B$_z$", "coolwarm")
data(3, xy, B_y0[:, 0], grid_x, grid_y, -50, 50)  # r"Y=0 B$_x$", "coolwarm")
data(4, xy, B_y0[:, 1], grid_x, grid_y, -50, 50)  # r"Y=0 B$_y$", "coolwarm")
data(5, xy, B_y0[:, 2], grid_x, grid_y, -50, 50)  # r"Y=0 B$_z$", "coolwarm")
for i in [3, 4, 5]:
    grid[i].set_xlabel(r"x (R$_M$)")
grid[0].set_ylabel(r"y (R$_M$)")
grid[3].set_ylabel(r"z (R$_M$)")
grid[0].set_title(r"B$_x$ Z=0")
grid[1].set_title(r"B$_y$ Z=0")
grid[2].set_title(r"B$_z$ Z=0")
grid[3].set_title(r"B$_x$ Y=0")
grid[4].set_title(r"B$_y$ Y=0")
grid[5].set_title(r"B$_z$ Y=0")

cb = fig.colorbar(
    cm.ScalarMappable(norm=Normalize(-50, 50), cmap="coolwarm"), cax=grid.cbar_axes[0]
)
cb.ax.set_title("B (nT)")

fig.set_size_inches(9, 6)
# when saving, specify the DPI
# plt.savefig("../../../Dropbox/Paper2/B_2d.png", dpi=600)
plt.show()


# # densidad

fig = plt.figure(1, (2.0, 3.0))
grid = AxesGrid(
    fig,
    111,  # similar to subplot(142)
    nrows_ncols=(2, 3),
    axes_pad=0.22,
    share_all=True,
    label_mode="L",
    cbar_location="right",
    cbar_mode="single",
)

data(0, xy, densidad_z0["H"], grid_x, grid_y, 0, 20, colormap="inferno")
data(1, xy, densidad_z0["heavies"], grid_x, grid_y, 0, 20, colormap="inferno")
data(2, xy, densidad_z0["e"], grid_x, grid_y, 0, 20, colormap="inferno")
data(3, xy, densidad_y0["H"], grid_x, grid_y, 0, 20, colormap="inferno")
data(4, xy, densidad_y0["heavies"], grid_x, grid_y, 0, 20, colormap="inferno")
data(5, xy, densidad_y0["e"], grid_x, grid_y, 0, 20, colormap="inferno")
for i in [3, 4, 5]:
    grid[i].set_xlabel(r"x (R$_M$)")
grid[0].set_ylabel(r"y (R$_M$)")
grid[3].set_ylabel(r"z (R$_M$)")
grid[0].set_title(r"H$^+$ dens. Z=0", fontsize=10)
grid[1].set_title(r"heavy ion dens. Z=0", fontsize=10)
grid[2].set_title(r"e$^-$ dens. Z=0", fontsize=10)
grid[3].set_title(r"H$^+$ dens. Y=0", fontsize=10)
grid[4].set_title(r"heavy ion dens. Y=0", fontsize=10)
grid[5].set_title(r"e$^-$ dens. Y=0", fontsize=10)

cb = figure.colorbar(
    cm.ScalarMappable(norm=Normalize(0, 20), cmap="inferno"), cax=grid.cbar_axes[0]
)
cb.ax.set_title("dens.\n(cm⁻³)", fontsize=10)
figure.set_size_inches(9, 6)
# when saving, specify the DPI
# plt.savefig("../../../Dropbox/Paper2/dens_2d.png", dpi=600)
plt.show()


# #  en escala log
# fig, ax = plt.subplots(2, 3)
# subplot_2d(x, y, np.log(densidad_z0["H"]), 0, 20, ax, 0, 0, "Z=0 H⁺ dens.", "inferno")
# subplot_2d(
#     x, y, np.log(densidad_z0["heavies"]), 0, 20, ax, 0, 1, "Z=0 heavy ion dens.", "inferno"
# )
# subplot_2d(x, y, np.log(densidad_z0["e"]), 0, 20, ax, 0, 2, "Z=0 e⁻ dens.", "inferno")
# subplot_2d(x, z, np.log(densidad_y0["H"]), 0, 20, ax, 1, 0, "Y=0 H⁺ dens.", "inferno")
# subplot_2d(
#     x, z, np.log(densidad_y0["heavies"]), 0, 20, ax, 1, 1, "Y=0 heavy ion dens.", "inferno"
# )
# subplot_2d(x, z, np.log(densidad_y0["e"]), 0, 20, ax, 1, 2, "Y=0 e⁻ dens.", "inferno")
# for i in [0, 1, 2]:
#     for j in [0, 1]:
#         plt.setp(ax[j, i].get_yticklabels(), visible=False)
#     plt.setp(ax[0, i].get_xticklabels(), visible=False)
#     ax[0, i].set_aspect("equal", "box")
#     ax[1, i].set_aspect("equal", "box")
# ax[0, 0].set_ylabel(r"y (R$_M$)")
# ax[1, 0].set_ylabel(r"z (R$_M$)")
# ax[1, 0].set_xlabel(r"x (R$_M$)")
# ax[1, 1].set_xlabel(r"x (R$_M$)")
# cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
# cb = fig.colorbar(cm.ScalarMappable(norm=Normalize(-1, 1), cmap="inferno"), cax=cbar_ax)
# cb.ax.set_title("dens. (cm⁻³)")
# plt.show()


#
# dif_vel_z0 = vz0H - vz0e
# dif_vel_y0 = vy0H - vy0e
# dv_H = np.array([rhoH_z0[i] * vz0H[i, :] for i in range(len(rhoH_z0))])
# dv_e = np.array([rhoe_z0[i] * vz0e[i, :] for i in range(len(rhoH_z0))])
# dv_Hy = np.array([rhoH_y0[i] * vy0H[i, :] for i in range(len(rhoH_z0))])
# dv_ey = np.array([rhoe_y0[i] * vy0e[i, :] for i in range(len(rhoH_z0))])
# fig, ax = plt.subplots(2, 2)
# subplot_2d(x, y, densidad_z0["H"], 0, 20, ax, 0, 0, "H")
# subplot_2d(x, y, densidad_z0["e"], 0, 20, ax, 1, 0, "e")
# subplot_2d(x, y, densidad_y0["H"], 0, 20, ax, 0, 1, "H y=0")
# subplot_2d(x, y, densidad_y0["e"], 0, 20, ax, 1, 1, "e y=0")
# plt.plot(X1, Y1, c="k")
# for i in np.arange(0, len(X), 5):
#     ax[0, 0].quiver(
#         X[i], Y[i], dv_H[i, 0], dv_H[i, 1], scale=5000, color="k", alpha=0.5
#     )
#     ax[1, 0].quiver(
#         X[i], Y[i], dv_e[i, 0], dv_e[i, 1], scale=5000, color="k", alpha=0.5
#     )
#     ax[0, 1].quiver(
#         X[i], Y[i], dv_Hy[i, 0], dv_Hy[i, 1], scale=5000, color="k", alpha=0.5
#     )
#     ax[1, 1].quiver(
#         X[i], Y[i], dv_ey[i, 0], dv_ey[i, 1], scale=5000, color="k", alpha=0.5
#     )
# for i in [0, 1]:
#     plt.setp(ax[0, i].get_xticklabels(), visible=False)
# for i in [0, 1]:
#     plt.setp(ax[i, 1].get_yticklabels(), visible=False)
# plt.show()


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
# plt.ylabel("z (RM)")
# plt.xlim([1, 1.6])
# plt.ylim([-1, 1])
# plt.show()
# plt.ylabel("z (RM)")
# plt.xlim([1, 1.6])
# plt.ylim([-1, 1])
# plt.show()
# plt.xlim([1, 1.6])
# plt.ylim([-1, 1])
# plt.show()
# plt.ylim([-1, 1])
# plt.show()
