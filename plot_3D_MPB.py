import numpy as np
import matplotlib.pyplot as plt
from funciones_plot import equal_axes

from mpl_toolkits.mplot3d import Axes3D  # de acá importo la proyección 3D

"""Plotea nuestra MPB y las fuerzas, órbitas, o lo que quiera"""

R = [1.082, -0.064, 0.515]
normal = [0.920, -0.302, 0.251]
EHall = [26.370, -8.439, 6.846]
Ecv = [-0.098, 0.1629, 1.3616]  # es el Ecv en la MS

L = 0.96
x0 = 0.78
e = 0.9
theta = np.linspace(0, np.pi * 3 / 4, 100)
phi = np.linspace(0, 2 * np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(THETA))
# r1 = L / (1 + e * np.cos(THETA))

X1 = x0 + r1 * np.cos(THETA)
Y1 = r1 * np.sin(THETA) * np.cos(PHI)
Z1 = r1 * np.sin(THETA) * np.sin(PHI)

# ahora plotea
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection="3d")
ax.set_xlabel(r"$X_{MSO} (R_m)$", fontsize=14)
ax.set_ylabel(r"$Y_{MSO} (R_m)$", fontsize=14)
ax.set_zlabel(r"$Z_{MSO} (R_m)$", fontsize=14)
ax.plot_surface(
    X1,
    Y1,
    Z1,
    rstride=4,
    cstride=4,
    alpha=0.5,
    # edgecolor="gray",
    cmap=plt.get_cmap("Blues_r"),
)
ax.scatter(R[0], R[1], R[2], label="MAVEN", color="k", s=40)

u, v = np.mgrid[0 : 2 * np.pi : 20j, 0 : np.pi : 10j]
ax.plot_wireframe(
    np.cos(u) * np.sin(v),
    np.sin(u) * np.sin(v),
    np.cos(v),
    color="#c1440e",
    linewidth=0.5,
)

ax.quiver(
    R[0],
    R[1],
    R[2],
    EHall[0],
    EHall[1],
    EHall[2],
    color="k",
    length=0.05,
    label="E Hall",
)

ax.quiver(
    R[0],
    R[1],
    R[2],
    Ecv[0],
    Ecv[1],
    Ecv[2],
    color="r",
    length=0.5,
    label="E cv",
)

plt.legend()
equal_axes(ax, X1, Y1, Z1)
plt.show()


theta = np.linspace(0, np.pi * 3 / 4, 100)
phi = np.linspace(0, 2 * np.pi, 100)

r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(theta))
# r1 = L / (1 + e * np.cos(THETA))

X = x0 + r1 * np.cos(theta)
Y = r1 * np.sin(theta) * np.cos(phi)
Z = r1 * np.sin(theta) * np.sin(phi)


Ryz = np.sqrt(R[1]**2 + R[2]**2)
Ecv_yz = np.sqrt(Ecv[1]**2 + Ecv[2]**2)
EHall_yz = np.sqrt(EHall[1]**2 + EHall[2]**2)


fig, ax = plt.subplots()
ax.plot()
ax.axis("equal")
ax.set_xlim(0, 2)
ax.set_ylim(0, 1.5)
circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
ax.add_artist(circle)

ax.plot(X, np.sqrt(Y**2 + Z**2))
ax.quiver(R[0], Ryz, Ecv[0], Ecv_yz, color='red', scale=20, label='Ecv')
ax.quiver(R[0], Ryz, EHall[0], EHall_yz, scale=100, label='Ehall')

ax.set_title("March 16th, 2016", fontsize=16)
ax.set_xlabel(r"$X_{MSO}$ ($R_M$)", fontsize=14)
ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)", fontsize=14)
plt.legend()
plt.show()
