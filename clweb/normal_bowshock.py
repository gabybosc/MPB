import numpy as np
import sys
from MVA_sin_spreadsheet import ajuste, normal_coplanar
from importar_datos import importar_mag, importar_swicfa
import matplotlib.pyplot as plt
import matplotlib.dates as md

sys.path.append("..")

from funciones import (
    angulo,
    ancho_mpb,
    corrientes,
    datenum,
    find_nearest,
    find_nearest_final,
    find_nearest_inicial,
    fechas,
    tiempos,
    donde,
    UTC_to_hdec,
)

"""
calcula la normal coplanar y del fit al bowshock
"""

year, month, day, doy = 2016, "03", 16, 76
ti, tf = 17.5, 18.5

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
t_c, t_f, density_c, density_f = importar_swicfa(year, month, day, ti, tf)

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

rampa = donde(t, UTC_to_hdec("18:03:30"))

theta_coplanar = angulo(n_coplanar, B_upstream) * 180 / np.pi

x0 = 0.64
e = 1.03

index = np.arange(fin_up, inicio_down, 100)


def normal_fit(posicion, x0=0.64, e=1.03):
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    R = posicion / 3390
    r0 = R - np.array([x0, 0, 0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))

    asc = L0 / (e ** 2 - 1)
    bsc = L0 / (e ** 2 - 1) ** (1 / 2)
    csc = x0 + L0 * e / (e ** 2 - 1)

    norm_vignes = np.array(
        [(R[0] + csc) * 2 / asc ** 2, R[1] * 2 / bsc ** 2, R[2] * 2 / bsc ** 2]
    )
    norm_vignes = norm_vignes / np.linalg.norm(norm_vignes)

    if norm_vignes[0] < 0:
        norm_vignes = -norm_vignes

    return norm_vignes


nn = [normal_fit(posicion[idx]) for idx in index]

aa = [180 - angulo(nn[idx], B_upstream) * 180 / 3.14 for idx in range(len(nn))]


theta = np.linspace(0, np.pi * 3 / 4, 100)
phi = np.linspace(0, 2 * np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

R = posicion[rampa, :] / 3390
r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))

norm_vignes = normal_fit(posicion[rampa])

print("theta_B vignes =", angulo(norm_vignes, B_upstream) * 180 / np.pi)

r1 = L0 / (1 + e * np.cos(theta))
x1 = x0 + r1 * np.cos(theta)
y1 = r1 * np.sin(theta) * np.cos(phi)
z1 = r1 * np.sin(theta) * np.sin(phi)

plt.plot(x1, np.sqrt(y1 ** 2 + z1 ** 2))
plt.quiver(
    R[0],
    np.sqrt(R[1] ** 2 + R[2] ** 2),
    norm_vignes[0],
    np.sqrt(norm_vignes[1] ** 2 + norm_vignes[2] ** 2),
)
plt.xlabel("x")
plt.ylabel("sqrt(y²+z²)")
plt.show()

"""
Figura de B, densidad, theta Bn
"""

tiempo_mag = np.array([np.datetime64(datenum(year, float(month), day, x)) for x in t])
tiempo_swica = np.array(
    [np.datetime64(datenum(year, float(month), day, x)) for x in t_c]
)
tiempo_swifa = np.array(
    [np.datetime64(datenum(year, float(month), day, x)) for x in t_f]
)
t_angulo = np.array([tiempo_mag[idx] for idx in index])

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.subplot2grid((3, 1), (0, 0))
ax1.xaxis.set_major_formatter(xfmt)
ax1.set_title(f"MAVEN MAG SWIA {year}-{month}-{day}")
ax1.plot(tiempo_mag, np.linalg.norm(B, axis=1), c="C0")
ax1.set_ylabel("|B| (nT)")

ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
ax2.semilogy(tiempo_swica, density_c, c="C0", label="SWICA")
ax2.semilogy(tiempo_swifa, density_f, c="C1", label="SWIFA")
ax2.set_ylabel("SWIA ion \n density (cm⁻³)")
ax2.legend(loc="upper left")

ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax2)
ax3.xaxis.set_major_formatter(xfmt)
ax3.plot(t_angulo, aa)
# ax3.set_ylim(ymin=0.8)
ax3.set_ylabel(r"$\theta_{Bn}$ (º)")
ax3.set_xlabel("Time (UTC)")

for ax in [ax2, ax1]:
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_xlim(tiempo_mag[25000], tiempo_mag[-25000])

figure = plt.gcf()  # get current figure
figure.set_size_inches(8, 9)
# when saving, specify the DPI
plt.show()


def alpha(x, y):
    """
    Calculates angle between two vectors based on its internal product.
    The angle is given in degrees and is restricted to the first cuadrant.

    Parameters:
        x,y: arrays

    Returns:
        a: float
           angle between x and y, in degrees
    """

    a = abs(
        np.arccos((np.dot(x, y)) / (np.linalg.norm(x) * np.linalg.norm(y)))
        * 180
        / (np.pi)
    )

    # restrict to first cuadrant
    if a > 90:
        a = abs(180 - a)

    return a


# normales formulas teo coplanaridad


def norm_coplanar(Bd, Bu, Vd, Vu):
    """
    Calcula todos los tipos de normales coplanares que hay (Analysis Methods Cap10).
    nB falla en 0 y 90
    nV es una formula aprox que vale para Mach muy grandes para angulos cerca de 0 y 90.
    """

    nB = (np.cross(np.cross(Bd, Bu), (Bd - Bu))) / (
        np.linalg.norm(np.cross(np.cross(Bd, Bu), (Bd - Bu)))
    )
    if (
        nB[0] < 0
    ):  # en SR MSO el versor x apunta hacia el sol, entonces para tener normal externa nx>0
        nB = -nB

    nBuV = (np.cross(np.cross(Bu, (Vd - Vu)), (Bd - Bu))) / (
        np.linalg.norm(np.cross(np.cross(Bu, (Vd - Vu)), (Bd - Bu)))
    )
    if nBuV[0] < 0:
        nBuV = -nBuV

    nBdV = (np.cross(np.cross(Bd, (Vd - Vu)), (Bd - Bu))) / (
        np.linalg.norm(np.cross(np.cross(Bd, (Vd - Vu)), (Bd - Bu)))
    )
    if nBdV[0] < 0:
        nBdV = -nBdV

    nBduV = (np.cross(np.cross((Bd - Bu), (Vd - Vu)), (Bd - Bu))) / (
        np.linalg.norm(np.cross(np.cross((Bd - Bu), (Vd - Vu)), (Bd - Bu)))
    )
    if nBduV[0] < 0:
        nBduV = -nBduV

    nV = (Vd - Vu) / np.linalg.norm(Vd - Vu)
    if nV[0] < 0:
        nV = -nV

    return nB, nBuV, nBdV, nBduV, nV
    nBduV = (np.cross(np.cross((Bd - Bu), (Vd - Vu)), (Bd - Bu))) / (
        np.linalg.norm(np.cross(np.cross((Bd - Bu), (Vd - Vu)), (Bd - Bu)))
    )
    if nBduV[0] < 0:
        nBduV = -nBduV

    nV = (Vd - Vu) / np.linalg.norm(Vd - Vu)
    if nV[0] < 0:
        nV = -nV

    return nB, nBuV, nBdV, nBduV, nV
