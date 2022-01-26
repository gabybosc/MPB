from importar_datos import importar_mag, importar_swica, importar_lpw, importar_STATIC
from matplotlib.widgets import MultiCursor
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.dates as md
import matplotlib.pyplot as plt
from cycler import cycler
import numpy as np
import pandas as pd
import sys

sys.path.append("..")

from funciones_plot import equal_axes, onpick1
from funciones import donde, datenum

path = "../../../datos/simulacion_chuanfei/"
datos_enteros = np.loadtxt(path + "sat_trajectory_HallOn_new2.sat", skiprows=2)

"""
Mismos análisis pero para los datos sobre la trayectoria de la nave en función
del tiempo
Si tengo dudas de cómo elegí los datos, al final del script hay una función
que hace los recortes y explica todo
"""

"""
it year mo dy hr mn sc msc X Y (0 a 9)
Z Rho Ux Uy Uz Bx By Bz Pe P (10 a 19)
HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP (20 a 29)
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP (30 a 39)
b1x b1y b1z e jx jy jz (40 a 46)
"""

# Datos de la simulación
pi = 1000
pf = 1250
mu0 = 4e-7 * np.pi  # T m / A

datos = datos_enteros[pi:pf]

x = datos[:, 8]
y_full = datos[:, 9]
z_full = datos[:, 10]

B = datos[:, 15:18]  # nT
b1 = datos[:, 40:43]
J = datos[:, -3:]  # uA/m2

presion = {
    "e-": datos[:, 18],
    "H+": datos[:, 24],
    "O+": datos[:, 34],
    "O2+": datos[:, 29],
    "CO2+": datos[:, 39],
}  # nPa

densities = {
    "e-": datos[:, 11],
    "H+": datos[:, 20],
    "O+": datos[:, 30],
    "O2+": datos[:, 25],
    "CO2+": datos[:, 35],
}  # mp/cc

velocidad = {
    "H+": datos[:, 12:15],
    "O+": datos[:, 31:34],
    "O2+": datos[:, 26:29],
    "CO2+": datos[:, 36:39],
}  # km/s


# Datos de MAVEN
mag, t, B_mag, posicion = importar_mag(2016, "03", 16, 17.7, 18.5)
STATIC, t_static, H_density, O_density, O2_density, CO2_density = importar_STATIC(
    2016, "03", 16, 17.7, 18.5
)

# Datos del análisis de MAVEN
R = [1.082, -0.064, 0.515]
normal = [0.920, -0.302, 0.251]
j_maven = 282  # nA/m²
v_x = -13000  # km/h la velocidad de MAVEN en x


t1, t2, t3, t4 = 18.2167, 18.2204, 18.235, 18.2476
t_up = t1 - 0.015
t_down = t4 + 0.015

# la explicación de esto está en la función al final de todo

r_simu = np.transpose([x, y_full, z_full]) * 3390
zi = donde(posicion[:, 2], r_simu[0, 2])
zf = donde(posicion[:, 2], r_simu[-1, 2])

posicion_cut = posicion[zi:zf]
t_cut = t[zi:zf]
B_cut = B_mag[zi:zf]
t_simu = np.linspace(t_cut[0], t_cut[-1], len(r_simu))


# los valores estos los elegi mirando los gráficos de la función ancho
ti_simu = t_simu[donde(x, 1.15)]
tf_simu = t_simu[donde(x, 1.01)]

ii = donde(t_simu, ti_simu)
jj = donde(t_simu, tf_simu)

lpw, t_lpw, e_density, flag = importar_lpw(2016, "03", 16, t_cut[0], t_cut[-1])
swia, t_swia, proton_density, sw_vel = importar_swica(
    2016, "03", 16, t_cut[0], t_cut[-1]
)

year = 2016
month = 3
day = 16

tiempo_mag = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_cut])
tiempo_simu = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_simu])
tiempo_swia = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swia])
tiempo_lpw = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_lpw])
tiempo_static = np.array(
    [np.datetime64(datenum(year, month, day, x)) for x in t_static]
)


idx_flag = [i for i in range(len(flag)) if flag[i] > 50]

normal = np.array([0.920, -0.302, 0.251])
ancho_mpb = np.dot(r_simu[jj] - r_simu[ii], normal)

posicion_cut = posicion[zi:zf]
t_cut = t[zi:zf]
B_cut = B_mag[zi:zf]
t_simu = np.linspace(t_cut[0], t_cut[-1], len(r_simu))


tt = donde(t_cut, t1)
ff = donde(t_cut, t4)
# ancho = np.dot(posicion_cut[tt] - posicion_cut[ff], normal)

j_media = np.mean(np.linalg.norm(J, axis=1)[ii:jj]) * 1e3

ii = donde(t_simu, ti_simu)
jj = donde(t_simu, tf_simu)
i_menos = donde(t_simu, ti_simu - 0.0125)
j_mas = donde(t_simu, tf_simu + 0.0125)

Bdown = np.mean(B[jj:j_mas], axis=0) * 1e-9
Bup = np.mean(B[i_menos:ii], axis=0) * 1e-9
J_salto = 1 / (mu0 * ancho_mpb * 1e3) * np.cross(normal, Bup - Bdown) * 1e9

print(f"el ancho de la MPB de la simu es {np.abs(ancho_mpb):.3g} km")
print(f"la corriente media en la MPB de la simu es {np.abs(j_media):.3g} nA/m²")
print(
    f"la corriente por salto en la MPB de la simu es {np.linalg.norm(J_salto):.3g} nA/m²"
)

# long inercial
paso = 20

density_mean = [
    np.mean(densities["H+"][i : i + paso]) for i in range(len(densities["H+"]) - paso)
]

ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km

print(
    f"la longitud inercial de iones en la zona Upstream es {np.mean(ion_length[ii:jj]):.3g} km"
)

# recorto los puntos de e_density que son malos
densidades_malas = [A for A in range(len(e_density)) if e_density[A] < 1]
e_density[densidades_malas] = np.nan  # van a aparecer en blanco estos puntos
"""
Ahora vienen los plots
"""
plt.rcParams.update({"font.size": 12})
plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()

ax1 = plt.subplot2grid((5, 1), (0, 0))
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(tiempo_mag, B_cut[:, 0], label="Bx MSO")
ax1.plot(tiempo_mag, B_cut[:, 1], label="By MSO")
ax1.plot(tiempo_mag, B_cut[:, 2], label="Bz MSO")
ax1.plot(tiempo_simu, B[:, 0], label="Bx MSO simulation")
ax1.plot(tiempo_simu, B[:, 1], label="By MSO simulation")
ax1.plot(tiempo_simu, B[:, 2], label="Bz MSO simulation")
ax1.set_ylabel("B components (nT)")
ax1.legend(loc="upper left")
ax1.set_title(f"MAVEN MAG LPW SWIA {year}-{month}-{day}")

ax2 = plt.subplot2grid((5, 1), (1, 0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, np.linalg.norm(B_cut, axis=1))
plt.plot(tiempo_simu, np.linalg.norm(B, axis=1))
plt.ylabel("|B| (nT)")

ax3 = plt.subplot2grid((5, 1), (2, 0), sharex=ax1)
ax3.xaxis.set_major_formatter(xfmt)
ax3.plot(tiempo_swia, np.linalg.norm(sw_vel, axis=1))
ax3.plot(tiempo_simu, np.linalg.norm(velocidad["H+"], axis=1))
ax3.set_ylabel("SW proton \n velocity (km/s)")

ax4 = plt.subplot2grid((5, 1), (3, 0), sharex=ax1)
ax4.xaxis.set_major_formatter(xfmt)
ax4.semilogy(tiempo_lpw, e_density)
ax4.semilogy(tiempo_simu, densities["e-"])
ax4.set_ylabel("Electron \n density (cm⁻³)")

ax5 = plt.subplot2grid((5, 1), (4, 0), sharex=ax1)
ax5.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_swia, proton_density)
plt.plot(tiempo_simu, densities["H+"])
ax5.set_ylabel("SW Proton \n density (cm⁻³)")
ax5.set_xlabel("Time (UTC)")
ax5.xaxis.set_label_coords(-0.05, -0.05)

for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.axvspan(
        xmin=tiempo_mag[donde(t_cut, ti_simu)],
        xmax=tiempo_mag[donde(t_cut, tf_simu)],
        facecolor="k",
        alpha=0.2,
    )
    ax.axvspan(
        xmin=tiempo_mag[donde(t_cut, t1)],
        xmax=tiempo_mag[donde(t_cut, t4)],
        facecolor="k",
        alpha=0.5,
    )

    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()

ax2.legend(["MAVEN", "Simulation", "MPB sim", "MPB MAVEN"], loc="upper left")
for ax in [ax1, ax2, ax3, ax4]:
    plt.setp(ax.get_xticklabels(), visible=False)

plt.show()


"""
plot de Corrientes, sólo simu
"""
plt.figure()
plt.plot(t_simu, J * 1e3)
plt.plot(t_simu, np.linalg.norm(J, axis=1) * 1e3)
plt.axvspan(xmin=ti_simu, xmax=tf_simu, color="m", alpha=0.2)
plt.axvspan(xmin=t1, xmax=t4, color="k", alpha=0.2)
plt.axvline(x=t2, color="k")
plt.axvline(x=t3, color="k")
# plt.axhline(y=-34, color="C0", linestyle="--")
# plt.axhline(y=-233, color="C1", linestyle="--")
# plt.axhline(y=-153, color="C2", linestyle="--")
# plt.axhline(y=238, color="C3", linestyle="--")
plt.legend(["x", "y", "z", "norm", "simu", "maven"])
plt.ylabel("J (nA/m2)")
plt.xlabel("t (hdec)")
plt.title("Corriente en volumen. La línea punteada es el cálculo de MAVEN")
plt.grid()
plt.ylim(ymin=-250, ymax=250)
plt.show()

"""
plot y cálculo de campos eléctricos, sólo simu. No tenemos gradp
"""
v_plus = np.zeros((len(t_simu), 3))
for ion in ["H+", "O+", "O2+", "CO2+"]:
    for i in range(3):
        v_plus[:, i] += densities[ion] * velocidad[ion][:, i] / densities["e-"]

v_SI = v_plus * 1e3  # m/s
B_SI = B * 1e-9  # T
n_SI = densities["H+"] * 1e6  # 1/m3
J_SI = J * 1e-6  # A/m2
e_SI = 1.6e-19  # C

Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)

plt.figure()
ax4 = plt.subplot2grid((2, 1), (0, 0))
ax5 = plt.subplot2grid((2, 1), (1, 0))

ax4.set_title("Electric fields")
ax4.plot(t_simu, Ecv * 1e3, ".")
ax4.set_ylabel("Ecv (mV/m)")
plt.setp(ax4.get_xticklabels(), visible=False)

ax5.plot(t_simu, Ehall * 1e3, ".")
ax5.set_ylabel("Ehall (mV/m)")
plt.setp(ax5.get_xticklabels(), visible=False)

for ax in [ax4, ax5]:
    ax.set_ylim([-5, 5])
    ax.grid()
    ax.axvspan(xmin=ti_simu, xmax=tf_simu, facecolor="k", alpha=0.2)

ax5.set_xlabel("x (RM)")
ax5.legend(["x", "y", "z", "MPB"])


plt.show()


fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()
ax1 = plt.subplot2grid((1, 1), (0, 0))
ax1.xaxis.set_major_formatter(xfmt)
ax1.semilogy(tiempo_static, H_density, ".", label="H")
ax1.semilogy(tiempo_static, O_density, ".", label="O")
ax1.semilogy(tiempo_static, O2_density, ".", label="O2")
ax1.semilogy(tiempo_static, CO2_density, ".", label="CO2")
ax1.semilogy(tiempo_simu, densities["H+"], c="C0", label="H simu")
ax1.semilogy(tiempo_simu, densities["O+"], c="C1", label="O simu")
ax1.semilogy(tiempo_simu, densities["O2+"], c="C2", label="O2 simu")
ax1.semilogy(tiempo_simu, densities["CO2+"], c="C3", label="CO2 simu")
ax1.legend(loc="upper left")
ax1.axvspan(
    xmin=tiempo_mag[donde(t_cut, ti_simu)],
    xmax=tiempo_mag[donde(t_cut, tf_simu)],
    facecolor="k",
    alpha=0.2,
)
ax1.set_ylim(ymin=0.1, ymax=1e5)
ax1.grid()
plt.show()


def check(posicion, B_mag):
    """Es para checkear que está bien el recorte"""

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

    plt.plot(x, y_full, z_full, label="Órbita", color="C1")
    plt.plot(
        posicion[:, 0] / 3390,
        posicion[:, 1] / 3390,
        posicion[:, 2] / 3390,
        label="new",
        color="C2",
        linestyle="--",
    )
    equal_axes(ax, X1, Y1, Z1)
    plt.show()

    """
    Ahora que sabemos que se superponen, pero que la simu es más cortita,
    hay que encontrar los puntos en los que coinciden inicial y final y
    recortar MAVEN ahí. Luego, hay que escribir eso en función de t.
    """

    r_simu = np.transpose([x, y_full, z_full]) * 3390

    xi = donde(posicion[:, 0], r_simu[0, 0])
    yi = donde(posicion[:, 1], r_simu[0, 1])
    zi = donde(posicion[:, 2], r_simu[0, 2])

    xf = donde(posicion[:, 0], r_simu[-1, 0])
    yf = donde(posicion[:, 1], r_simu[-1, 1])
    zf = donde(posicion[:, 2], r_simu[-1, 2])

    print("iniciales", xi, yi, zi, "finales", xf, yf, zf)

    # como son todos valores parecidos/iguales, podemos recortar ahí
    posicion_cut = posicion[zi:zf]
    t_cut = t[zi:zf]
    B_cut = B_mag[zi:zf]

    """
    Ahora vamos a ver si alcanza con igualar el t_inicial y el t_final.
    """

    t_simu = np.linspace(t_cut[0], t_cut[-1], len(r_simu))

    plt.figure()
    plt.plot(t_simu, r_simu)
    plt.plot(t_cut, posicion_cut, "--")
    plt.show()

    # se superponen perfectamente!!
    # se superponen perfectamente!!


def ancho(x, b1):
    dBx_dx = np.gradient(b1[:, 0], np.abs(x[0] - x[2]) * 3390e3)
    dBy_dx = np.gradient(b1[:, 1], np.abs(x[0] - x[2]) * 3390e3)
    dBz_dx = np.gradient(b1[:, 2], np.abs(x[0] - x[2]) * 3390e3)
    B_norm = np.linalg.norm(b1, axis=1)

    dBdx = np.transpose(np.vstack((dBx_dx, dBy_dx, dBz_dx)))

    # plt.figure()
    # plt.plot(x, dBdx, ".")
    # plt.plot(x, np.linalg.norm(dBdx, axis=1), ".")
    # plt.legend(["x", "y", "z", "norm"])
    # plt.title("Derivada de B")
    # plt.show()

    fig = plt.figure()
    fig.subplots()
    plt.title("Campo magnético y su derivada en x")

    ax1 = plt.subplot2grid((4, 1), (0, 0))
    plt.plot(x, B_norm, ".")
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel("|B|")
    ax1.grid()

    ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
    ax2.plot(x, np.linalg.norm(dBdx, axis=1), ".")
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.set_ylabel("dB/dx")
    ax2.set_xlabel("x (RM)")
    ax2.grid()

    ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
    ax3.plot(x, np.linalg.norm(Ehall * 1e3, axis=1), ".")
    ax3.set_ylabel("Ehall")
    ax3.set_xlabel("x (RM)")
    ax3.set_ylim([0, 5])
    ax3.grid()

    ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)
    ax4.plot(x, np.linalg.norm(J * 1e3, axis=1), ".")
    ax4.set_ylabel("J")
    ax4.set_xlabel("x (RM)")
    ax4.set_ylim([0, 50])
    ax4.grid()

    multi = MultiCursor(fig.canvas, (ax1, ax2, ax3, ax4), color="black", lw=1)
    plt.show()


paso = 20

density_mean = [
    np.mean(densities["H+"][i : i + paso]) for i in range(len(densities["H+"]) - paso)
]

ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km

b1_norm = np.linalg.norm(b1, axis=1)
x_cut = x[ii:jj]
B_cut = b1_norm[ii:jj]

coef = np.polynomial.polynomial.polyfit(x_cut, B_cut, deg=1)

# MPR = np.mean(b1_norm[donde(x, 1.16) : donde(x, 1.174)])

MPB_inicio = donde(coef[0] + x * coef[1], 1.0)
MPB_fin = donde(x, 1.11)

f, (ax1, ax) = plt.subplots(2, 1, sharex=True)
ax1.plot(x, b1_norm, ".")
ax1.plot(x, coef[0] + x * coef[1], label="fit")
# ax1.axhline(y=MPR, c="k")
ax1.axhline(y=b1_norm[donde(x, 1.22)], c="k")
ax1.set_title("MPB thickness")
ax1.set_ylabel("|B|")
ax1.set_ylim([0, 50])
ax1.grid()

ax2 = ax.twinx()
ax.plot(x[:-paso], ion_length, ".")
ax.set_ylabel("proton inertial length (km)")
ax2.set_ylabel("proton inertial length (RM)")

ax.set_ylim([50, 250])
ax2.set_ylim([50 / 3390, 250 / 3390])
# set an invisible artist to twin axes
# to prevent falling back to initial values on rescale events
ax2.plot([], [])
ax.set_xlabel("X (RM)")
# ax.set_xlim([1.15, 1.6])

ax.grid()
plt.show()
