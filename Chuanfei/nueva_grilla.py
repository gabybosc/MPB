import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.interpolate
from scipy.stats import multivariate_normal
import seaborn as sns
import numpy.ma as ma
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from pykrige.ok import OrdinaryKriging

sns.set()

path = "../../../datos/simulacion_chuanfei/"
datos_y0 = np.loadtxt(path + "y=0_HallOn_new2.gz")
datos_z0 = np.loadtxt(path + "z=0_HallOn_new2.gz")

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


def recortar(x, y, arr):
    a = arr[donde(x, 1.1) : donde(x, 1.3)]
    s = np.array([a[j] for j in range(len(y)) if np.abs(y[j]) < 0.1])
    return s


yy = y[donde(x, 1.1) : donde(x, 1.3)]
zz = z[donde(x, 1.1) : donde(x, 1.3)]
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

# defino la grilla


def regrid(X, Y, z, Nx=50, Ny=50, model="linear"):
    """
    X, Y son las coordenadas de la grilla, z es el valor en cada punto,
    Nx, Ny es el nuevo tamaño de grilla
    """
    OK = OrdinaryKriging(X, Y, z, variogram_model=model)
    dx = (max(X) - min(X)) / (Nx - 1)
    dy = (max(Y) - min(Y)) / (Ny - 1)
    xgr = min(X) + dx * np.arange(0, Nx)
    ygr = min(Y) + dy * np.arange(0, Ny)
    out, ss = OK.execute("grid", xgr, ygr)  # out es la grilla, ss es la varianza
    return out


v_regrid_xl = regrid(X, Y, vz0H[:, 0], model="linear")
v_regrid_xg = regrid(X, Y, vz0H[:, 0], model="gaussian")
v_regrid_yl = regrid(X, Y, vz0H[:, 1], model="linear")
v_regrid_yg = regrid(X, Y, vz0H[:, 1], model="gaussian")
rho_regrid = regrid(X, Y, rhoH_z0, model="gaussian")

plt.figure()
plt.imshow(v_regrid_yl)
plt.figure()
plt.imshow(v_regrid_yg)
plt.show()


nueva_grilla = np.column_stack(
    [x.flat, y.flat]
)  # Create a (N, 2) array of (x, y) pairs.
new_grid_x, new_grid_y = np.mgrid[
    min(X) : max(X) : 50j, min(Y) : max(Y) : 50j
]  # grilla de 50x50
fig, ax = plt.subplots(1, 1)
cf = ax.pcolormesh(new_grid_x, new_grid_y, np.transpose(rho_regrid))
for i in np.arange(0, 50):
    ax.quiver(
        new_grid_x[i],
        new_grid_y[i],
        v_regrid_xl.flatten()[i],
        v_regrid_yg.flatten()[i],
        scale=5000,
        color="k",
        alpha=0.5,
    )
fig.colorbar(cf, cmap="inferno")

xy = np.column_stack([X.flat, Y.flat])  # Create a (N, 2) array of (x, y) pairs.
# Interpolate and generate heatmap:
grid_x, grid_y = np.mgrid[min(X) : max(X) : 1000j, min(Y) : max(Y) : 1000j]

# interpolation method can be linear or cubic
grid_z = scipy.interpolate.griddata(xy, rhoH_z0, (grid_x, grid_y), method="linear")
fig, ax = plt.subplots(1, 1)
cf = ax.pcolormesh(grid_x, grid_y, ma.masked_invalid(grid_z),)
# ax.set_title(titulo)
for i in np.arange(0, len(X), 5):
    # for i in range(len(X)):
    ax.quiver(X[i], Y[i], vz0H[i, 0], vz0H[i, 1], scale=5000, color="k", alpha=0.5)
ax.set_xlim([1.1, 1.3])
ax.set_ylim([-0.1, 0.1])
fig.colorbar(cf, cmap="inferno")
plt.show()
