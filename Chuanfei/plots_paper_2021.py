import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from mpl_toolkits.mplot3d import Axes3D
import sys
from cycler import cycler
from importar_datos import importar_mag, importar_swica, importar_lpw, importar_STATIC

sys.path.append("..")

from funciones_plot import equal_axes, onpick1
from funciones import (
    donde,
    datenum,
    fechas,
    find_nearest_inicial,
    find_nearest_final,
    tiempos,
)

plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#cdcdcd"],
)

path = "../../../datos/simulacion_chuanfei/"
datos = np.loadtxt(path + "ejex_new2_+.gz")  # high resolution

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
e_SI = 1.6e-19  # C
g = 3.7  # Mars surface gravity, m/s²

"""
r1 r2 rho Hrho Orho O2rho CO2rho H_vx H_vy H_vz
Bx By Bz b1x b1y b1z Pe P HP OP
O2P CO2P jx jy jz gradx grady gradz O_vx O_vy
O_vz O2_vx O2_vy O2_vz CO2_vx CO2_vy CO2_vz
"""

x = datos[:, 0]  # RM
z = datos[:, 1]  # RM
B = datos[:, 10:13]  # nT
b1 = datos[:, 13:16]  # nT
J = datos[:, 22:25]  # ua/m2
grad_p = datos[:, 25:28]  # nPa/m
presion = {
    "e-": datos[:, 16],
    "H+": datos[:, 18],
    "O+": datos[:, 19],
    "O2+": datos[:, 20],
    "CO2+": datos[:, 21],
}  # nPa
densities = {
    "e-": datos[:, 2],
    "H+": datos[:, 3],
    "O+": datos[:, 4],
    "O2+": datos[:, 5],
    "CO2+": datos[:, 6],
}  # Mp/cc
velocities = {
    "H+": datos[:, 7:10],
    "O+": datos[:, 28:31],
    "O2+": datos[:, 31:34],
    "CO2+": datos[:, 34:],
}  # km/s

B_norm = np.linalg.norm(B, axis=1)

v_plus = np.zeros((len(x), 3))
for ion in ["H+", "O+", "O2+", "CO2+"]:
    for i in range(3):
        v_plus[:, i] += densities[ion] * velocities[ion][:, i] / densities["e-"]

v_SI = v_plus * 1e3  # m/s
B_SI = B * 1e-9  # T
n_SI = densities["H+"] * 1e6  # 1/m3
J_SI = J * 1e-6  # A/m2
grad_p_SI = grad_p * 1e-9

Ecv = np.array([-np.cross(v_SI[i], B_SI[i]) for i in range(len(B))])  # V/m
Ehall = np.array(
    [1 / (e_SI * n_SI[i]) * np.cross(J_SI[i], B_SI[i]) for i in range(len(B))]
)
Ep = np.array([-1 / (e_SI * n_SI[i]) * grad_p_SI[i, :] for i in range(len(grad_p))])

P_heavy = presion["O+"] + presion["O2+"] + presion["CO2+"]
P_B = np.linalg.norm(b1, axis=1) ** 2 * 1e-9 / (2 * mu0)
P_ram = 1.67e-6 * densities["H+"] * velocities["H+"][:, 0] ** 2  # nPa
P_total = P_heavy + P_B + presion["e-"] + P_ram + presion["H+"]

rho_heavies = np.zeros(len(x))
for ion in ["O+", "O2+", "CO2+"]:
    rho_heavies += densities[ion]


"""Ancho MPB: hay una función al fondo que uso para elegir los límites que
tomé acá"""

x_cut = x[donde(x, 1.174) : donde(x, 1.22)]
B_cut = B_norm[donde(x, 1.174) : donde(x, 1.22)]

coef = np.polynomial.polynomial.polyfit(x_cut, B_cut, deg=1)

MPR = np.mean(B_norm[donde(x, 1.16) : donde(x, 1.174)])

MPB_inicio = donde(coef[0] + x * coef[1], MPR)
MPB_fin = donde(x, 1.22)

ancho_mpb = (x[MPB_fin] - x[MPB_inicio]) * 3390  # km

# valores medios de los campos en la MPB
EH_medio = np.mean(np.linalg.norm(Ehall[MPB_inicio:MPB_fin], axis=1)) * 1e3
Ecv_medio = np.mean(np.linalg.norm(Ecv[MPB_inicio:MPB_fin], axis=1)) * 1e3
Ep_medio = np.mean(np.linalg.norm(Ep[MPB_inicio:MPB_fin], axis=1)) * 1e3

# campos up y down
ti_up = 1.16
tf_up = 1.174
ti_down = 1.22
tf_down = 1.234
Bup = np.mean(B[donde(x, ti_up) : donde(x, tf_up)], axis=0) * 1e-9
Bdown = np.mean(B[donde(x, ti_down) : donde(x, tf_down)], axis=0) * 1e-9

# longitud inercial de protones
paso = 20

density_mean = [
    np.mean(densities["H+"][i : i + paso]) for i in range(len(densities["H+"]) - paso)
]

ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km

fig, ax = plt.subplots()
ax2 = ax.twinx()

"""
Giroradio térmico en la MPB
"""

v_th = np.sqrt(
    np.mean(presion["H+"][MPB_inicio:MPB_fin], axis=0)
    * 1e-21
    / (np.mean(densities["H+"][MPB_inicio:MPB_fin], axis=0) * mp)
)  # km/s

B_medio = np.mean(B[MPB_inicio:MPB_fin], axis=0)

rg = mp * v_th / (e_SI * np.linalg.norm(B_medio)) * 1e9  # km

# punto a punto
v_th = np.sqrt(presion["H+"] * 1e-21 / (densities["H+"] * mp))  # km/s

rg = [
    mp * v_th[i] / (e_SI * np.linalg.norm(B[i, :])) * 1e9 for i in range(len(B))
]  # km

sw_i = donde(x, 1.6)
sw_f = donde(x, 1.8)
BS = donde(x, 1.55)
ms_i = donde(x, 1.3)
ms_f = donde(x, 1.5)
mpr_i = donde(x, 1.05)
mpr_f = donde(x, 1.15)

B_sw = np.mean(B[sw_i:sw_f, :], axis=0)
B_bs = B[BS, :]
B_ms = np.mean(B[ms_i:ms_f, :], axis=0)
B_mpb = np.mean(B[MPB_inicio:MPB_fin, :], axis=0)
B_mpr = np.mean(B[mpr_i:mpr_f, :], axis=0)

# los vectores en 2D: hago el sqrt (y²+z²)
B_mpr_2d = [B_mpr[0], np.sqrt(B_mpr[1] ** 2 + B_mpr[2] ** 2)]
B_mpb_2d = [B_mpb[0], np.sqrt(B_mpb[1] ** 2 + B_mpb[2] ** 2)]
B_ms_2d = [B_ms[0], np.sqrt(B_ms[1] ** 2 + B_ms[2] ** 2)]
B_bs_2d = [B_bs[0], np.sqrt(B_bs[1] ** 2 + B_bs[2] ** 2)]
B_sw_2d = [B_sw[0], np.sqrt(B_sw[1] ** 2 + B_sw[2] ** 2)]

plt.plot(x, B)
plt.axvspan(xmin=x[MPB_inicio], xmax=x[MPB_fin], alpha=0.2)
plt.quiver(1.1, 40, B_mpr_2d[0], B_mpr_2d[1], scale=500)
plt.quiver(1.2, 40, B_mpb_2d[0], B_mpb_2d[1], scale=250)
plt.quiver(1.4, 10, B_ms_2d[0], B_ms_2d[1], scale=100)
plt.quiver(1.55, 10, B_bs_2d[0], B_bs_2d[1], scale=50)
plt.quiver(1.7, 10, B_sw_2d[0], B_sw_2d[1], scale=40)
plt.legend(["Bx", "By", "Bz"])
plt.ylabel("|B| (nT)")
plt.xlabel("x MSO")
plt.title("las flechas son el campo en Bx, sqrt(By²+Bz²)")
plt.ylim([-10, 70])
plt.grid()
plt.show()

# plots

plt.rcParams.update({"font.size": 11})
plt.figure()
ax1 = plt.subplot2grid((4, 1), (0, 0))
ax4 = plt.subplot2grid((4, 1), (1, 0))
ax5 = plt.subplot2grid((4, 1), (2, 0))
ax6 = plt.subplot2grid((4, 1), (3, 0))

plt.subplots_adjust(hspace=0.0)

for ax in [ax1, ax4, ax5, ax6]:
    ax.set_xlim([1.15, 1.3])
    ax.axvspan(xmin=ti_up, xmax=x[MPB_inicio], facecolor="#428AE0", alpha=0.4)  # down
    ax.axvspan(
        xmin=x[MPB_inicio], xmax=x[MPB_fin], facecolor="#79B953", alpha=0.4
    )  # mpb
    ax.axvspan(xmin=x[MPB_fin], xmax=tf_down, facecolor="#FE6779", alpha=0.4)  # up

    ax.grid()
ax1.legend(["Downstream", "MPB", "Upstream"], loc="upper right")

ax1.plot(x, np.linalg.norm(B, axis=1), c="C3")
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("|B| (nT)")
ax1.set_ylim([0, 50])
ax1.set_title(
    "Magnetic field, electric current and electric field\n from simulation over the Mars-Sun line"
)

ax4.plot(x, J[:, 0] * 1e3, label=r"J$_x$")
ax4.plot(x, J[:, 1] * 1e3, label=r"J$_y$")
ax4.plot(x, J[:, 2] * 1e3, label=r"J$_z$")
ax4.plot(x, np.linalg.norm(J, axis=1) * 1e3, label="|J|")
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.set_ylabel("J (nA/m²)")
ax4.set_ylim([-120, 120])
ax4.legend(loc="upper right")

ax5.plot(x, np.linalg.norm(Ecv, axis=1) * 1e3, label=r"E$_{cv}$", c="C3")
ax5.plot(x, np.linalg.norm(Ehall, axis=1) * 1e3, label=r"E$_H$", c="C4")
ax5.plot(x, np.linalg.norm(Ep, axis=1) * 1e3, label=r"E$_p$", c="C5")
ax5.set_ylabel("E (mV/m)")
ax5.set_ylim([-0.5, 5])
plt.setp(ax5.get_xticklabels(), visible=False)
ax5.legend(loc="upper right", prop={"size": 10})

ax6.plot(x, Ehall * 1e3)
ax6.set_ylabel(r"E$_H$ (mV/m)")
ax6.set_xlabel(r"X ($R_M$)")
ax6.legend(["x", "y", "z"], loc="upper right", prop={"size": 10})
ax6.set_ylim([-1.1, 4.5])


figure = plt.gcf()  # get current figure
figure.set_size_inches(8, 9)
# when saving, specify the DPI
plt.savefig("../../../Dropbox/Paper2/simulación_ejex.png", dpi=300)


"""
Trayectoria temporal
"""
path = "../../../datos/simulacion_chuanfei/"
datos_enteros = np.loadtxt(path + "sat_trajectory_HallOn_new2.sat", skiprows=2)

# Datos de la simulación
pi = 1000
pf = 1250
mu0 = 4e-7 * np.pi  # T m / A

datos_tray = datos_enteros[pi:pf]

x_tray = datos_tray[:, 8]
y_tray = datos_tray[:, 9]
z_tray = datos_tray[:, 10]

B_tray = datos_tray[:, 15:18]  # nT
b1_tray = datos_tray[:, 40:43]
J_tray = datos_tray[:, -3:]  # uA/m2

presion_tray = {
    "e-": datos_tray[:, 18],
    "H+": datos_tray[:, 24],
    "O+": datos_tray[:, 34],
    "O2+": datos_tray[:, 29],
    "CO2+": datos_tray[:, 39],
}  # nPa

densities_tray = {
    "e-": datos_tray[:, 11],
    "H+": datos_tray[:, 20],
    "O+": datos_tray[:, 30],
    "O2+": datos_tray[:, 25],
    "CO2+": datos_tray[:, 35],
}  # mp/cc

velocidad_tray = {
    "H+": datos_tray[:, 12:15],
    "O+": datos_tray[:, 31:34],
    "O2+": datos_tray[:, 26:29],
    "CO2+": datos_tray[:, 36:39],
}  # km/s


# Datos de MAVEN
mag, t, B_mag, posicion = importar_mag(2016, "03", 16, 17.7, 18.5)

# Datos del análisis de MAVEN
R = [1.082, -0.064, 0.515]
normal = [0.920, -0.302, 0.251]
j_maven = 282  # nA/m²
v_x = -13000  # km/h la velocidad de MAVEN en x


t1, t2, t3, t4 = 18.2167, 18.2204, 18.235, 18.2476
t_up = t1 - 0.015
t_down = t4 + 0.015

# la explicación de esto está en la función al final de todo

r_simu = np.transpose([x_tray, y_tray, z_tray]) * 3390
zi = donde(posicion[:, 2], r_simu[0, 2])
zf = donde(posicion[:, 2], r_simu[-1, 2])

posicion_cut = posicion[zi:zf]
t_cut = t[zi:zf]
B_cut = B_mag[zi:zf]
t_simu = np.linspace(t_cut[0], t_cut[-1], len(r_simu))

lpw, t_lpw, e_density, flag = importar_lpw(2016, "03", 16, t_cut[0], t_cut[-1])
swia, t_swia, proton_density, sw_vel = importar_swica(
    2016, "03", 16, t_cut[0], t_cut[-1]
)
STATIC, t_static, H_density, O_density, O2_density, CO2_density = importar_STATIC(
    2016, "03", 16, 17.7, 18.5
)

# los valores estos los elegi mirando los gráficos de la función ancho
ti_simu = t_simu[donde(x_tray, 1.15)]
tf_simu = t_simu[donde(x_tray, 1.01)]

ii = donde(t_simu, ti_simu)
jj = donde(t_simu, tf_simu)

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
# idx_flag = [i for i in range(len(flag)) if flag[i] > 50]

normal = np.array([0.920, -0.302, 0.251])
ancho_mpb = np.dot(r_simu[jj] - r_simu[ii], normal)

tt = donde(t_cut, t1)
ff = donde(t_cut, t4)
# ancho = np.dot(posicion_cut[tt] - posicion_cut[ff], normal)

j_media = np.mean(J[ii:jj]) * 1e3

i_menos = donde(t_simu, ti_simu - 0.0125)
j_mas = donde(t_simu, tf_simu + 0.0125)

Bup = np.mean(B_tray[jj:j_mas], axis=0) * 1e-9
Bdown = np.mean(B_tray[i_menos:ii], axis=0) * 1e-9
J_salto = 1 / (mu0 * ancho_mpb * 1e3) * np.cross(normal, Bup - Bdown) * 1e9

densidades_malas = [A for A in range(len(e_density)) if e_density[A] < 1]
e_density[densidades_malas] = np.nan

"""
Ahora vienen los plots
"""
plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

plt.rcParams.update({"font.size": 12})

fig = plt.figure()
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax2 = plt.gca()

ax2 = plt.subplot2grid((4, 1), (0, 0))
ax2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, np.linalg.norm(B_cut, axis=1), c="C0")
plt.plot(tiempo_simu, np.linalg.norm(B_tray, axis=1), c="C1")
plt.ylabel("|B| (nT)")
ax2.set_title(f"MAVEN MAG LPW SWIA STATIC {year}-{month}-{day}")

ax1 = plt.subplot2grid((4, 1), (3, 0), sharex=ax2)
ax1.xaxis.set_major_formatter(xfmt)
ax1.semilogy(tiempo_static, H_density, ".", label=r"H$^+$")
ax1.semilogy(tiempo_static, O_density, ".", label=r"O$^+$")
ax1.semilogy(tiempo_static, O2_density, ".", label=r"O$_2^+$")
ax1.semilogy(tiempo_static, CO2_density, ".", label=r"CO$_2^+$")
ax1.semilogy(tiempo_simu, densities_tray["H+"], c="C0")
ax1.semilogy(tiempo_simu, densities_tray["O+"], c="C1")
ax1.semilogy(tiempo_simu, densities_tray["O2+"], c="C2")
ax1.semilogy(tiempo_simu, densities_tray["CO2+"], c="C3")
ax1.legend(loc="upper left")
ax1.set_ylim(ymin=0.1, ymax=1e5)
ax1.set_ylabel("H⁺ and heavies\ndensities (cm⁻³)")
ax1.set_xlabel("Time (UTC)")

ax4 = plt.subplot2grid((4, 1), (1, 0), sharex=ax2)
ax4.xaxis.set_major_formatter(xfmt)
ax4.semilogy(tiempo_lpw, e_density, c="C0")
ax4.semilogy(tiempo_simu, densities_tray["e-"], c="C1")
ax4.set_ylabel("Electron \n density (cm⁻³)")

ax5 = plt.subplot2grid((4, 1), (2, 0), sharex=ax2)
ax5.xaxis.set_major_formatter(xfmt)
plt.semilogy(tiempo_swia, proton_density, c="C0")
plt.semilogy(tiempo_simu, densities_tray["H+"], c="C1")
ax5.set_ylim(ymin=1)
ax5.set_ylabel("SW ion \n density (cm⁻³)")
# ax5.xaxis.set_label_coords(-0.05, -0.05)

for ax in [ax1, ax2, ax4, ax5]:
    ax.axvspan(
        xmin=tiempo_mag[donde(t_cut, ti_simu)],
        xmax=tiempo_mag[donde(t_cut, tf_simu)],
        facecolor="#79B953",
        alpha=0.4,
    )
    ax.axvspan(
        xmin=tiempo_mag[donde(t_cut, t1)],
        xmax=tiempo_mag[donde(t_cut, t4)],
        facecolor="k",
        alpha=0.5,
    )

    ax.set_xlim(tiempo_mag[12000], tiempo_mag[-17000])
    ax.grid()

ax2.legend(["MAVEN", "Simulation", "MPB simulation", "MPB MAVEN"], loc="upper left")
for ax in [ax2, ax5, ax4]:
    plt.setp(ax.get_xticklabels(), visible=False)

figure = plt.gcf()  # get current figure
figure.set_size_inches(8, 9)
# when saving, specify the DPI
plt.savefig("../../../Dropbox/Paper2/MAVEN_vs_simu_log.png", dpi=600)
plt.show()


# el plot de la órbita


# def datos_fijos(year, month, day, ti, tf):
#     # path a los datos desde la laptop
#     path = f"../../../datos/clweb/{year}-{month}-{day}/"
#     mag = np.loadtxt(path + "MAG.asc")
#
#     M = len(mag[:, 0])  # el numero de datos
#
#     hh = mag[:, 3]
#     mm = mag[:, 4]
#     ss = mag[:, 5]
#
#     t = hh + mm / 60 + ss / 3600  # hdec
#
#     posicion = np.zeros((M, 3))
#     for i in range(9, 12):
#         posicion[:, i - 9] = mag[:, i] / 3390
#
#     inicio = np.where(t == find_nearest_inicial(t, ti))[0][0]
#     fin = np.where(t == find_nearest_final(t, tf))[0][0]
#
#     posicion_cut = posicion[inicio:fin, :]
#     t_cut = t[inicio:fin]
#
#     return (t_cut, posicion_cut, year, month, day)
#
#
# t, posicion_cut, year, month, day = datos_fijos(2016, "03", 16, 17.5, 18.4)
#
#
# def BS_MPB(L, e, x0):
#     THETA = np.linspace(0, np.pi * 3 / 4, 100)
#     PHI = np.linspace(0, 2 * np.pi, 100)
#
#     r1 = L / (1 + e * np.cos(THETA))
#
#     X1 = x0 + r1 * np.cos(THETA)
#     Y1 = r1 * np.sin(THETA) * np.cos(PHI)
#     Z1 = r1 * np.sin(THETA) * np.sin(PHI)
#     yz = np.sqrt(Y1 ** 2 + Z1 ** 2)
#     return (X1, yz)
#
#
# def marte(x_bs, yz_bs, x_mpb, yz_mpb):
#     fig, ax = plt.subplots()
#     fig.set_size_inches(7, 7)
#     ax.plot()
#     ax.plot(x_bs, yz_bs, color="#FF1493", linestyle="-.")
#     ax.plot(x_mpb, yz_mpb, color="#79B953", linestyle="-.")
#     ax.axis("equal")
#     ax.set_xlim(0, 2)
#     ax.set_ylim(0, 2)
#     circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
#     ax.add_artist(circle)
#     ax.set_title("MAVEN MSO coordinates", fontsize=16)
#     ax.set_xlabel(r"$X_{MSO}$ ($R_M$)", fontsize=14)
#     ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)", fontsize=14)
#
#
# def orbitas(posicion, year, month, day, puntos):
#     x = posicion[:, 0]
#     y = posicion[:, 1]
#     z = posicion[:, 2]
#     proyeccion = np.sqrt(y ** 2 + z ** 2)
#     fig = plt.gcf()
#     ax = fig.gca()
#     ax.plot(x, proyeccion, label=f"{day}-{month}-{year}", linewidth=0.5, c="k")
#     ax.scatter(x[puntos], proyeccion[puntos], c="k")
#
#
# x_bs, yz_bs = BS_MPB(2.04, 1.03, 0.64)
# x_mpb, yz_mpb = BS_MPB(0.96, 0.9, 0.78)
#
# # t, posicion_cut, year, month, day = datos()
# puntos = [donde(t, aa) for aa in [18.33, 18.16, 18, 17.83, 17.66]]
#
# marte(x_bs, yz_bs, x_mpb, yz_mpb)
# orbitas(posicion_cut, year, month, day, puntos)
#
# plt.show(block=False)
# # plt.savefig(f"../../../Dropbox/AGU2021/orbita.png", dpi=300)
