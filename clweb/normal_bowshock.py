import numpy as np
import sys
from MVA_sin_spreadsheet import ajuste, normal_coplanar
from importar_datos import importar_mag, importar_lpw
import matplotlib.pyplot as plt

sys.path.append("..")

from funciones import (
    angulo,
    ancho_mpb,
    corrientes,
    find_nearest,
    find_nearest_final,
    find_nearest_inicial,
    fechas,
    tiempos,
    donde,
    UTC_to_hdec,
)
from funciones_metodos import normal_fit, ajuste_conico

"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
"""

year, month, day, doy = 2016, "03", 16, 76
ti, tf = 17.5, 18.5

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

"""
Podemos probar cambiando los intervalos up/down y ver si varía mucho la normal
coplanar.
"""

inicio_up = donde(t, UTC_to_hdec("17:44:00"))
fin_up = donde(t, UTC_to_hdec("17:50:00"))
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0)  # nT

inicio_down = donde(t, UTC_to_hdec("18:06:00"))
fin_down = donde(t, UTC_to_hdec("18:10:00"))
B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0)  # nT

n_coplanar = -normal_coplanar(B_upstream, B_downstream)

theta_bn = angulo(n_coplanar, B_upstream) * 180 / np.pi

x0 = 0.64
e = 1.03

index = np.arange(fin_up, inicio_down, 100)
nn = [normal_fit(posicion, idx, 0.64, 1.03)[0] for idx in index]

angulo(nn, B_upstream) * 180 / np.pi


theta = np.linspace(0, np.pi * 3 / 4, 100)
phi = np.linspace(0, 2 * np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

R = posicion[53400, :] / 3390
r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(THETA))


"""
Hay algo raro en este fit igual, porque no parece ser normal a la superficie...
e = 1.03 me dice que es una hipérbola, no una parábola. Entonces tengo que
cambiar el cálculo de asc, bsc, csc.
"""

# ahora plotea
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection="3d")
ax.set_xlabel(r"$X_{MSO} (R_m)$")
ax.set_ylabel(r"$Y_{MSO} (R_m)$")
ax.set_zlabel(r"$Z_{MSO} (R_m)$")
ax.scatter(R[0], R[1], R[2], label="MAVEN", color="k", s=40)
X1 = x0 + r1 * np.cos(THETA)
Y1 = r1 * np.sin(THETA) * np.cos(PHI)
Z1 = r1 * np.sin(THETA) * np.sin(PHI)
ax.plot_surface(
    X1,
    Y1,
    Z1,
    rstride=4,
    cstride=4,
    alpha=0.5,
    edgecolor="none",
    cmap=plt.get_cmap("Blues_r"),
)

asc = L0 / (1 - e ** 2)  # semieje mayor
bsc = np.sqrt(asc * L0)
csc = e * asc - x0  # donde está centrada. Hay que ver el signo

norm_vignes = np.array(
    [(R[0] + csc) * 2 / asc ** 2, R[1] * 2 / bsc ** 2, R[2] * 2 / bsc ** 2]
)
norm_vignes = -norm_vignes / np.linalg.norm(norm_vignes)

ax.quiver(
    R[0],
    R[1],
    R[2],
    norm_vignes[0],
    norm_vignes[1],
    norm_vignes[2],
    color="b",
    length=0.5,
    label="Normal del ajuste",
)
ax.quiver(
    R[0],
    R[1],
    R[2],
    n_coplanar[0],
    n_coplanar[1],
    n_coplanar[2],
    color="k",
    length=0.5,
    label="Normal coplanar",
)

ax.legend()

plt.show()
