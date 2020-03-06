import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mpl
from importar_datos import importar_mag, importar_lpw, importar_swea, importar_swia
from cycler import cycler
import sys

sys.path.append("..")
from funciones import find_nearest, datenum, Bpara_Bperp, tiempos, UTC_to_hdec

np.set_printoptions(precision=4)

"""
Este script plotea mag, swea, swia y lpw en la región de interés
Es la fig principal del grl.

"""

year, month, day, doy = 2016, "03", 16, 76
ti, tf = 17.85, 18.4

# path = f'../../../datos/clweb/{year}-{month}-{day}/'
path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la desktop.
datos = np.loadtxt("../outputs/t1t2t3t4.txt")

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
t1, t2, t3, t4 = 18.2167, 18.2204, 18.235, 18.2476
t1t2t3t4 = np.array([t1, t2, t3, t4])

t_up = t1 - 0.015
t_down = t4 + 0.015

zoom_inicial = np.where(t == find_nearest(t, t1 - 0.05))[0][0]
zoom_final = np.where(t == find_nearest(t, t4 + 0.05))[0][0]

Bnorm = np.linalg.norm(B, axis=1)

mag_low = np.loadtxt(path + "mag_1s.sts", skiprows=160)
tlow = mag_low[:, 6]  # el dia decimal
tlow = (tlow - int(doy)) * 24  # para que me de sobre la cantidad de horas

Mlow = np.size(tlow)  # el numero de datos
# el campo
Blow = np.zeros((Mlow, 3))
for i in range(7, 10):
    Blow[:, i - 7] = mag_low[:, i]


B_para, B_perp_norm, t_plot = Bpara_Bperp(Blow, tlow, t[0], t[-1])

# ######### SWEA

swea, t_swea, energias, energy, JE = importar_swea(year, month, day, ti, tf)

# ####### densidad SWIA

swia, t_swia, density = importar_swia(year, month, day, ti, tf)

# ###### densidad electrones
lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)

# ############ tiempos UTC
year = int(year)
month = int(month)
day = int(day)

tiempo_mag = np.array(
    [np.datetime64(datenum(year, month, day, x)) for x in t]
)  # datenum es una función mía
tiempo_swea = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swea])
tiempo_swia = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swia])
tiempo_lpw = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_lpw])
tiempo_low = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_plot])

tm1 = np.where(t == find_nearest(t, t1))
tm2 = np.where(t == find_nearest(t, t2))
tm3 = np.where(t == find_nearest(t, t3))
tm4 = np.where(t == find_nearest(t, t4))
tm_up = np.where(t == find_nearest(t, t_up))
tm_down = np.where(t == find_nearest(t, t_down))
tmva = np.where(t == find_nearest(t, UTC_to_hdec("18:13:33")))
tbs = np.where(t == find_nearest(t, UTC_to_hdec("18:02:00")))
tmpr = np.where(t == find_nearest(t, UTC_to_hdec("18:19:00")))

tiempo_lim = [tiempo_mag[tm1], tiempo_mag[tm2], tiempo_mag[tm3], tiempo_mag[tm4]]


#####plot grl zoom

mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)
fig = plt.figure(
    1, figsize=(20, 10)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

axz1 = plt.gca()

axz1.xaxis.set_major_formatter(xfmt)
axz1 = plt.subplot2grid((2, 1), (0, 0))
axz1.plot(tiempo_mag, B[:, 0], label="Bx MSO")
axz1.plot(tiempo_mag, B[:, 1], label="By MSO")
axz1.plot(tiempo_mag, B[:, 2], label="Bz MSO")
plt.setp(axz1.get_xticklabels(), visible=False)
axz1.set_ylabel("B components (nT)")
axz1.set_title(f"MAVEN MAG 2016 March 16th")

axz2 = plt.subplot2grid((2, 1), (1, 0), sharex=axz1)
axz2.xaxis.set_major_formatter(xfmt)
axz2.plot(tiempo_mag, Bnorm)
axz2.set_ylabel("|B| (nT)")

axz2.set_xlabel("Time (UTC)")


axz1.axvspan(
    xmin=tiempo_mag[tm_up][0], xmax=tiempo_mag[tm1][0], facecolor="#581845", alpha=0.5
)
axz1.axvspan(
    xmin=tiempo_mag[tmva][0], xmax=tiempo_mag[tm4][0], facecolor="#FFC300", alpha=0.75
)
axz1.axvspan(
    xmin=tiempo_mag[tm4][0], xmax=tiempo_mag[tm_down][0], facecolor="#C70039", alpha=0.5
)

axz2.axvspan(
    xmin=tiempo_mag[tm_up][0],
    xmax=tiempo_mag[tm1][0],
    facecolor="#581845",
    alpha=0.5,
    label="Upstream",
)
axz2.axvspan(
    xmin=tiempo_mag[tmva][0],
    xmax=tiempo_mag[tm4][0],
    facecolor="#FFC300",
    alpha=0.75,
    label="MVA",
)
axz2.axvspan(
    xmin=tiempo_mag[tm4][0],
    xmax=tiempo_mag[tm_down][0],
    facecolor="#C70039",
    alpha=0.5,
    label="Downstream",
)


for ax in [axz1, axz2]:
    ax.set_xlim(tiempo_mag[zoom_inicial], tiempo_mag[zoom_final])
    ax.grid()
    ax.legend()
    for xc in tiempo_lim:
        ax.axvline(x=xc, color="k", linewidth=1.5)

# plt.tight_layout()
plt.show(block=False)


###### funciones para el plot que se repiten
def regiones(ax, tiempo_mag, tm1, tm4, tm_up, tm_down, tbs, tmpr):
    ax.axvspan(
        xmin=tiempo_mag[tbs][0], xmax=tiempo_mag[tm1][0], facecolor="#FE6779", alpha=0.6
    )  # Magnetosheath
    ax.axvspan(
        xmin=tiempo_mag[tm1][0], xmax=tiempo_mag[tm4][0], facecolor="#79B953", alpha=0.6
    )  # MPB
    ax.axvspan(
        xmin=tiempo_mag[tm4][0],
        xmax=tiempo_mag[tmpr][0],
        facecolor="#428AE0",
        alpha=0.5,
    )  # MPR


#########La fig de la tesis pero sin LPW
mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)
fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()

ax1.xaxis.set_major_formatter(xfmt)
ax1 = plt.subplot2grid((5, 1), (0, 0))
ax1.plot(tiempo_mag, B[:, 0], label="Bx MSO")
ax1.plot(tiempo_mag, B[:, 1], label="By MSO")
ax1.plot(tiempo_mag, B[:, 2], label="Bz MSO")
ax1.set_ylabel("B components (nT)")
ax1.set_title(f"MAVEN MAG SWEA SWIA {year}-{month}-{day}")

ax2 = plt.subplot2grid((5, 1), (1, 0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, Bnorm)
plt.ylabel("|B| (nT)")

ax3 = plt.subplot2grid((5, 1), (2, 0), sharex=ax1)
ax3.xaxis.set_major_formatter(xfmt)
ax3.plot(tiempo_low, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_low, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Relative variation \n of B")

ax4 = plt.subplot2grid((5, 1), (3, 0), sharex=ax1)
ax4.xaxis.set_major_formatter(xfmt)
ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[0]
    JEplot = JE[index]
    plt.semilogy(tiempo_swea, JEplot, label=f"{energia} eV")

ax5 = plt.subplot2grid((5, 1), (4, 0), sharex=ax1)
ax5.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_swia, density)
ax5.set_ylabel("SW proton \n density (cm⁻³)")
ax5.set_xlabel("Time (UTC)")


for ax in [ax1, ax2, ax3, ax4]:
    regiones(ax, tiempo_mag, tm1, tm4, tm_up, tm_down, tbs, tmpr)
    plt.setp(ax.get_xticklabels(), visible=False)
ax5.axvspan(
    xmin=tiempo_mag[tbs][0],
    xmax=tiempo_mag[tm1][0],
    facecolor="#FE6779",
    alpha=0.6,
    label="Magnetosheath",
)
ax5.axvspan(
    xmin=tiempo_mag[tm1][0],
    xmax=tiempo_mag[tm4][0],
    facecolor="#79B953",
    alpha=0.6,
    label="MPB",
)
ax5.axvspan(
    xmin=tiempo_mag[tm4][0],
    xmax=tiempo_mag[tmpr][0],
    facecolor="#428AE0",
    alpha=0.5,
    label="MPR",
)


for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()
    ax.legend()
    ax.axvline(x=tiempo_mag[tbs][0], color="c", linewidth=1.5)
    for xc in tiempo_lim:
        ax.axvline(x=xc, color="k", linewidth=1.5)


# plt.tight_layout()
plt.show()
