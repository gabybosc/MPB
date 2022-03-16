import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.interpolate

# from scipy.stats import multivariate_normal
import seaborn as sns
import numpy.ma as ma
import matplotlib.cm as cm
from matplotlib.colors import Normalize

sns.set()

path = "../../../datos/simulacion_chuanfei/"
datos_y0 = np.loadtxt(path + "y=0_HallOn_new2.gz")
datos_z0 = np.loadtxt(path + "z=0_HallOn_new2.gz")

sys.path.append("..")
from funciones import donde
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
grad_p_y0 = datos_y0[:, 25:28]  # nPa/m

B_z0 = datos_z0[:, 10:13]  # nT
b1_z0 = datos_z0[:, 13:16]  # nT
J_z0 = datos_z0[:, 22:25]  # ua/m2
grad_p_z0 = datos_z0[:, 25:28]  # nPa/m

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

P_heavy_z0 = presion_z0["O"] + presion_z0["O2"] + presion_z0["CO2"]
P_mag_z0 = np.linalg.norm(b1_z0, axis=1) ** 2 * 1e-9 / (2 * mu0)  # sin corticales
P_mag_cort_z0 = np.linalg.norm(B_z0, axis=1) ** 2 * 1e-9 / (2 * mu0)  # con corticales
P_dyn_z0 = 1.67e-6 * densidad_z0["H"] * velocidad_z0["H"][:, 0] ** 2  # nPa
P_total_z0 = P_heavy_z0 + P_mag_z0 + presion_z0["e"] + P_dyn_z0 + presion_z0["H"]

beta_z0 = (P_heavy_z0 + presion_z0["H"]) / P_mag_z0
beta_str_z0 = (P_heavy_z0 + presion_z0["H"] + P_dyn_z0) / P_mag_z0

P_heavy_y0 = presion_y0["O"] + presion_y0["O2"] + presion_y0["CO2"]
P_mag_y0 = np.linalg.norm(b1_y0, axis=1) ** 2 * 1e-9 / (2 * mu0)  # sin corticales
P_mag_cort_y0 = np.linalg.norm(B_y0, axis=1) ** 2 * 1e-9 / (2 * mu0)  # con corticales
P_dyn_y0 = 1.67e-6 * densidad_y0["H"] * velocidad_y0["H"][:, 0] ** 2  # nPa
P_total_y0 = P_heavy_y0 + P_mag_y0 + presion_y0["e"] + P_dyn_y0 + presion_y0["H"]

beta_y0 = (P_heavy_y0 + presion_y0["H"]) / P_mag_y0
beta_str_y0 = (P_heavy_y0 + presion_y0["H"] + P_dyn_y0) / P_mag_y0

# beta_cort = (P_heavy + presion_z0["H"]) / P_mag_cort
# beta_str_cort = (P_heavy + presion_z0["H"] + P_dyn) / P_mag_cort

densidad_z0_heavies = densidad_z0["O"] + densidad_z0["O2"] + densidad_z0["CO2"]
density_ratio = densidad_z0["H"] / densidad_z0_heavies
mass_ratio = densidad_z0["H"] / (
    densidad_z0["O"] * 16 + densidad_z0["O2"] * 32 + densidad_z0["CO2"] * 44
)

v_plus_z0 = np.zeros((len(x), 3))
v_plus_y0 = np.zeros((len(x), 3))
for ion in ["H", "O", "O2", "CO2"]:
    for i in range(3):
        v_plus_z0[:, i] += densidad_z0[ion] * velocidad_z0[ion][:, i] / densidad_z0["e"]
        v_plus_y0[:, i] += densidad_y0[ion] * velocidad_y0[ion][:, i] / densidad_y0["e"]

SI_z0 = {
    "v": v_plus_z0 * 1e3,
    "B": B_z0 * 1e-9,
    "n": densidad_z0["H"] * 1e6,
    "J": J_z0 * 1e-6,
    "grad_p": grad_p_z0 * 1e-9,
}
SI_y0 = {
    "v": v_plus_y0 * 1e3,
    "B": B_y0 * 1e-9,
    "n": densidad_y0["H"] * 1e6,
    "J": J_y0 * 1e-6,
    "grad_p": grad_p_y0 * 1e-9,
}

Ecv_y0 = np.array(
    [-np.cross(SI_y0["v"][i], SI_y0["B"][i]) for i in range(len(B_y0))]
)  # V/m
Ehall_y0 = np.array(
    [
        1 / (e_SI * SI_y0["n"][i]) * np.cross(SI_y0["J"][i], SI_y0["B"][i])
        for i in range(len(B_y0))
    ]
)
Ep_y0 = np.array(
    [-1 / (e_SI * SI_y0["n"][i]) * SI_y0["grad_p"][i, :] for i in range(len(grad_p_y0))]
)

Ecv_z0 = np.array(
    [-np.cross(SI_z0["v"][i], SI_z0["B"][i]) for i in range(len(B_z0))]
)  # V/m
Ehall_z0 = np.array(
    [
        1 / (e_SI * SI_z0["n"][i]) * np.cross(SI_z0["J"][i], SI_z0["B"][i])
        for i in range(len(B_z0))
    ]
)
Ep_z0 = np.array(
    [-1 / (e_SI * SI_z0["n"][i]) * SI_z0["grad_p"][i, :] for i in range(len(grad_p_z0))]
)


"""
probar rotarlo a MSE pero voy a necesitar tener los datos 3D
con x que sea igual a la MSO
El B upstream que uso que sea el mismo que le di a Chuanfei
El E va a quedar más prolijo
"""


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
    ax[i, j].set_xlim([1.1, 1.3])
    ax[i, j].set_ylim([-0.5, 0.5])
    ax[i, j].set_title(titulo)


# betas
fig, ax = plt.subplots(2, 2)
subplot_2d(x, y, beta_z0, 0, 2, ax, 0, 0, "beta z=0", "Set1")
subplot_2d(x, y, beta_str_z0, 0, 2, ax, 0, 1, "beta* z=0", "Set1")
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1]:
        ax[0, k].quiver(x[i], y[i], B_z0[i, 0], B_z0[i, 1], color="k", alpha=0.5)
subplot_2d(x, z, beta_y0, 0, 2, ax, 1, 0, "beta y=0", "coolwarm")
subplot_2d(x, z, beta_str_y0, 0, 2, ax, 1, 1, "beta* y=0", "coolwarm")
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1]:
        ax[1, k].quiver(x[i], z[i], B_y0[i, 0], B_y0[i, 2], color="k", alpha=0.5)
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

xx = x[donde(x, 1) : donde(x, 1.3)]  # recorto para los x entre 1 y 1.3
yy = y[donde(x, 1) : donde(x, 1.3)]
B_y0y = B_y0[donde(x, 1) : donde(x, 1.3)]

Y = [j for j in yy if np.abs(j) < 0.5]  # recorto donde y está entre -0.5 y 0.5
X = [xx[j] for j in range(len(yy)) if np.abs(yy[j]) < 0.5]
BB = np.array([B_y0y[j] for j in range(len(yy)) if np.abs(yy[j]) < 0.5])

fig, ax = plt.subplots(2, 2)
for i in np.arange(0, len(BB), 5):
    ax[1, 1].quiver(X[i], Y[i], BB[i, 0], BB[i, 2], color="k", alpha=0.5)
plt.show()


fig, ax = plt.subplots(2, 3)
subplot_2d(x, z, Ecv_y0[:, 0] * 1e3, -2, 2, ax, 0, 0, "Ecv_x y=0", "coolwarm")
subplot_2d(x, z, Ecv_y0[:, 1] * 1e3, -2, 2, ax, 0, 1, "Ecv_y0 y=0", "coolwarm")
subplot_2d(x, z, Ecv_y0[:, 2] * 1e3, -2, 2, ax, 0, 2, "Ecv_z0 y=0", "coolwarm")
subplot_2d(x, y, Ecv_z0[:, 0] * 1e3, -2, 2, ax, 1, 0, "Ecv_x z=0", "coolwarm")
subplot_2d(x, y, Ecv_z0[:, 1] * 1e3, -2, 2, ax, 1, 1, "Ecv_y0 z=0", "coolwarm")
subplot_2d(x, y, Ecv_z0[:, 2] * 1e3, -2, 2, ax, 1, 2, "Ecv_z0 z=0", "coolwarm")

for i in [0, 1, 2]:
    plt.setp(ax[0, i].get_xticklabels(), visible=False)
    ax[1, i].set_xlabel("x (RM)")
for i in [0, 1]:
    plt.setp(ax[i, 1].get_yticklabels(), visible=False)
    plt.setp(ax[i, 2].get_yticklabels(), visible=False)
ax[0, 0].set_ylabel("z (RM)")
ax[1, 0].set_ylabel("y (RM)")
cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
fig.colorbar(cm.ScalarMappable(norm=Normalize(-2, 2), cmap="coolwarm"), cax=cbar_ax)
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1, 2]:
        ax[0, k].quiver(x[i], z[i], Ecv_y0[i, 0], Ecv_y0[i, 2], color="k", alpha=0.5)
        ax[1, k].quiver(x[i], y[i], Ecv_z0[i, 0], Ecv_z0[i, 1], color="k", alpha=0.5)
plt.show()


fig, ax = plt.subplots(2, 3)
subplot_2d(x, z, Ehall_y0[:, 0] * 1e3, -2, 2, ax, 0, 0, "Ehall_x y=0", "coolwarm")
subplot_2d(x, z, Ehall_y0[:, 1] * 1e3, -2, 2, ax, 0, 1, "Ehall_y0 y=0", "coolwarm")
subplot_2d(x, z, Ehall_y0[:, 2] * 1e3, -2, 2, ax, 0, 2, "Ehall_z0 y=0", "coolwarm")
subplot_2d(x, y, Ehall_z0[:, 0] * 1e3, -2, 2, ax, 1, 0, "Ehall_x z=0", "coolwarm")
subplot_2d(x, y, Ehall_z0[:, 1] * 1e3, -2, 2, ax, 1, 1, "Ehall_y0 z=0", "coolwarm")
subplot_2d(x, y, Ehall_z0[:, 2] * 1e3, -2, 2, ax, 1, 2, "Ehall_z0 z=0", "coolwarm")

for i in [0, 1, 2]:
    plt.setp(ax[0, i].get_xticklabels(), visible=False)
    ax[1, i].set_xlabel("x (RM)")
for i in [0, 1]:
    plt.setp(ax[i, 1].get_yticklabels(), visible=False)
    plt.setp(ax[i, 2].get_yticklabels(), visible=False)
ax[0, 0].set_ylabel("z (RM)")
ax[1, 0].set_ylabel("y (RM)")
cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
fig.colorbar(cm.ScalarMappable(norm=Normalize(-2, 2), cmap="coolwarm"), cax=cbar_ax)
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1, 2]:
        ax[0, k].quiver(
            x[i], z[i], Ehall_y0[i, 0], Ehall_y0[i, 2], color="k", alpha=0.5
        )
        ax[1, k].quiver(
            x[i], y[i], Ehall_z0[i, 0], Ehall_z0[i, 1], color="k", alpha=0.5
        )
plt.show()


fig, ax = plt.subplots(2, 3)
subplot_2d(x, z, Ep_y0[:, 0] * 1e3, -2, 2, ax, 0, 0, "Ep_x y=0", "coolwarm")
subplot_2d(x, z, Ep_y0[:, 1] * 1e3, -2, 2, ax, 0, 1, "Ep_y0 y=0", "coolwarm")
subplot_2d(x, z, Ep_y0[:, 2] * 1e3, -2, 2, ax, 0, 2, "Ep_z0 y=0", "coolwarm")
subplot_2d(x, y, Ep_z0[:, 0] * 1e3, -2, 2, ax, 1, 0, "Ep_x z=0", "coolwarm")
subplot_2d(x, y, Ep_z0[:, 1] * 1e3, -2, 2, ax, 1, 1, "Ep_y0 z=0", "coolwarm")
subplot_2d(x, y, Ep_z0[:, 2] * 1e3, -2, 2, ax, 1, 2, "Ep_z0 z=0", "coolwarm")

for i in [0, 1, 2]:
    plt.setp(ax[0, i].get_xticklabels(), visible=False)
    ax[1, i].set_xlabel("x (RM)")
for i in [0, 1]:
    plt.setp(ax[i, 1].get_yticklabels(), visible=False)
    plt.setp(ax[i, 2].get_yticklabels(), visible=False)
ax[0, 0].set_ylabel("z (RM)")
ax[1, 0].set_ylabel("y (RM)")
cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
fig.colorbar(cm.ScalarMappable(norm=Normalize(-2, 2), cmap="coolwarm"), cax=cbar_ax)
for i in np.linspace(0, 12800, 258):
    i = int(i)
    for k in [0, 1, 2]:
        ax[0, k].quiver(x[i], z[i], Ep_y0[i, 0], Ep_y0[i, 2], color="k", alpha=0.5)
        ax[1, k].quiver(x[i], y[i], Ep_z0[i, 0], Ep_z0[i, 1], color="k", alpha=0.5)
plt.show()
