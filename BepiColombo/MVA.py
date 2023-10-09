import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from leer_datos import importar_bepi
from funciones_bepi import tiempos_UTC
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import sys

sys.path.append("../")

from funciones import (
    SZA,
    error,
    donde,
    Mij,
    Bpara_Bperp,
    SZA,
)

# from funciones_metodos import (
#     plot_bootstrap,
#     bootstrap,
# )
from funciones_plot import hodograma


"""
Hace el MVA y nada más. Es para correr rápido buscando la hoja de corriente.
"""


np.set_printoptions(precision=4)

year, month, day = 2021, "08", 10  # fechas()
ti_MVA, tf_MVA = 13.911111, 13.922222
t1, t2, t3, t4 = [13.8974, 13.9078, 13.9283, 13.9469]


t, B, posicion = importar_bepi(t1 - 0.5, t4 + 0.5)

# ti_MVA = t2  # comentar si quiero elegir los propios
# tf_MVA = t3

inicio = donde(t, ti_MVA)
fin = donde(t, tf_MVA)

#################
B_cut = B[inicio:fin]

M_ij = Mij(B_cut)

# ahora quiero los autovectores y autovalores
[lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

# Los ordeno de mayor a menor
idx = lamb.argsort()[::-1]
lamb = lamb[idx]
x = x[:, idx]
# ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
x1 = x[:, 0]
x2 = x[:, 1]
x3 = x[:, 2]

av = np.concatenate([x1, x2, x3])

if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
    x3 = -x3
if any(np.cross(x1, x2) - x3) > 0.01:
    print("Cambio el signo de x1 para que los av formen terna derecha")
    x1 = -x1

# las proyecciones
B1 = np.dot(B_cut, x1)
B2 = np.dot(B_cut, x2)
B3 = np.dot(B_cut, x3)


# el B medio
B_medio_vectorial = np.mean(B_cut, axis=0)
altitud = np.mean(posicion[inicio:fin, -1])
SZAngle = SZA(posicion[inicio:fin], 0)

B_norm_medio = np.linalg.norm(B_medio_vectorial)

hodograma(B1, B2, B3)

# el error
phi, delta_B3 = error(lamb, B_cut, x)

###############

n_fit = np.array([0.8938, 0.4485])
# B3_fit = np.dot(B_cut, n_fit)  # esto sirve para una normal 3D
###############
# Bootstrap

# N_boot = 1000
# normal_boot, phi_boot, delta_B3_boot, out, out_phi = bootstrap(N_boot, B_cut)

# muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

# B3_boot = np.dot(B_cut, normal_boot)

# #######
# # Errores
# if phi[2, 1] > phi[2, 0]:
#     error_normal = phi[2, 1] * 57.2958
# else:
#     error_normal = phi[2, 0] * 57.2958
#     # quiero ver si el error más grande es phi31 o phi32

# if sigma31 > sigma32:
#     error_boot = sigma31
# else:
#     error_boot = sigma32

n = np.array([x3[0], np.sqrt(x3[1] ** 2 + x3[2] ** 2)])
angulo_mva = np.arccos(np.clip(np.dot(n_fit, n), -1.0, 1.0))


print(f"SZA = {SZAngle:.3g}º y altitud = {int(altitud)}km")
print(f"MVA entre los tiempos {ti_MVA} y {tf_MVA}")
print(f"Cociente de lambdas = {lamb[1]/lamb[2]:.4g}")
print(
    f"El ángulo entre las normales 2D de MVA y del fit es {angulo_mva * 180/np.pi:.3g}º"
)

ti = t1 - 0.15
tf = t4 + 0.15

B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, ti, tf)

t, B, pos = importar_bepi(13.5, 14.1)
Bnorm = np.linalg.norm(B, axis=1)
pos_RV = pos / 6050


Bpara, Bperp, tpara = Bpara_Bperp(B, t, 13.5, 14.1)

yy = 2021
mm = 8
dd = 10
tiempo_mag = tiempos_UTC(yy, mm, dd, t)
tiempo_paraperp = tiempos_UTC(yy, mm, dd, tpara)

outs = [13.8807, 13.8933, 13.9094, 13.929]
MPB = tiempos_UTC(yy, mm, dd, outs)

fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")

ax1 = plt.gca()
ax1 = plt.subplot2grid((3, 1), (0, 0))
ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
for ax in [ax1, ax2, ax3]:
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()

ax1.plot(tiempo_mag, Bnorm, linewidth=0.5)
ax1.set_ylabel("|B| (nT)")
ax1.set_title(f"Bepi-Colombo MAG 2021-08-10")

ax2.plot(tiempo_mag, B[:, 0], label="Bx VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 1], label="By VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 2], label="Bz VSO", linewidth=0.5)
ax2.set_ylabel("B components (nT)")

ax3.plot(tiempo_paraperp, Bpara, linewidth=0.5, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_paraperp, Bperp, "-.", linewidth=0.5, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Relative variation \n of B")
ax3.set_xlabel("Tiempo (UTC)")
ax3.set_ylim([-0.1, 1])


for ax in [ax2, ax3]:
    ax.legend(loc="upper left")
for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
for ax in [ax1, ax2, ax3]:
    ax.axvspan(xmin=MPB[1], xmax=MPB[2], facecolor="#79B953", alpha=0.5)
    ax.axvspan(xmin=MPB[0], xmax=MPB[1], facecolor="#cdcdcd", alpha=0.7)
    ax.axvspan(xmin=MPB[2], xmax=MPB[3], facecolor="#cdcdcd", alpha=0.7)
    plt.axvline(x=ti_MVA, color="r", linewidth=1)
    plt.axvline(x=tf_MVA, color="r", linewidth=1)
    # en un radio de 10 min de la MPB
    ax.set_xlim([MPB[0] - np.timedelta64(10, "m"), MPB[-1] + np.timedelta64(10, "m")])


plt.show()
