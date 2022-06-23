import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_t1t2t3t4, importar_mag, importar_swea, importar_swia
from funciones import (
    fechas,
    tiempos,
    SZA,
    error,
    donde,
    Mij,
    Bpara_Bperp,
    SZA,
)
from funciones_metodos import (
    ajuste_conico,
    plot_bootstrap,
    bootstrap,
)
from funciones_plot import hodograma
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm


"""
Hace el MVA y nada más. Es para correr rápido buscando la hoja de corriente.
"""


np.set_printoptions(precision=4)

year, month, day, doy = 2016, "03", 20, 80  # fechas()
# ti_MVA = input("hora aproximada\n")
ti_MVA, tf_MVA = tiempos()
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day, int(ti_MVA))


mag, t, B, posicion = importar_mag(year, month, day, t1 - 0.5, t4 + 0.5)

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

plt.show(block=False)
# el error
phi, delta_B3 = error(lamb, B_cut, x)

###############
# fit
orbita = posicion[donde(t, t1 - 1) : donde(t, t4 + 1)] / 3390  # radios marcianos

index = donde(
    t[inicio:fin], (t2 + t3) / 2
)  # el tiempo en el medio de la hoja de corriente
x0 = 0.78
e = 0.9
normal_fit, X1, Y1, Z1, R, L0 = ajuste_conico(posicion[inicio:fin], index, orbita, x3)

B3_fit = np.dot(B_cut, normal_fit)

###############
# Bootstrap

N_boot = 1000
normal_boot, phi_boot, delta_B3_boot, out, out_phi = bootstrap(N_boot, B_cut)

muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

B3_boot = np.dot(B_cut, normal_boot)

#######
# Errores
if phi[2, 1] > phi[2, 0]:
    error_normal = phi[2, 1] * 57.2958
else:
    error_normal = phi[2, 0] * 57.2958
    # quiero ver si el error más grande es phi31 o phi32

if sigma31 > sigma32:
    error_boot = sigma31
else:
    error_boot = sigma32

angulo_mva = np.arccos(np.clip(np.dot(normal_fit, x3), -1.0, 1.0))


print(f"SZA = {SZAngle:.3g}º y altitud = {int(altitud)}km")
print(f"MVA entre los tiempos {ti_MVA} y {tf_MVA}")
print(f"Cociente de lambdas = {lamb[1]/lamb[2]:.4g}")
print(f"El ángulo entre las normales de MVA y del fit es {angulo_mva * 180/np.pi:.3g}º")

ti = t1 - 0.15
tf = t4 + 0.15

B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, t1, t4)

swea, t_swea, energia, flux_plot = importar_swea(year, month, day, t1 - 0.5, t4 + 0.5)
swia, t_swia, density, temp, vel = importar_swia(year, month, day, t1 - 0.5, t4 + 0.5)

plt.clf()  # clear figure
fig = plt.figure(
    4, constrained_layout=True
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.93, bottom=0.07, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
fig.set_size_inches(15, 10)  # con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((3, 2), (0, 0))
plt.plot(t_plot, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
plt.plot(t_plot, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / B")
plt.setp(ax1.get_xticklabels(), visible=False)
for xc in [t1, t2, t3, t4]:
    plt.axvline(x=xc, color="k", linewidth=1)
plt.axvline(x=ti_MVA, color="r", linewidth=1)
plt.axvline(x=tf_MVA, color="r", linewidth=1)
ax1.set_ylabel(r"|$\Delta B$|/ B")
ax1.grid()
ax1.legend()

ax4 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
ax4.plot(t, B)
for xc in [t1, t2, t3, t4]:
    plt.axvline(x=xc, color="k", linewidth=1)
plt.axvline(x=ti_MVA, color="r", linewidth=1)
plt.axvline(x=tf_MVA, color="r", linewidth=1)
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.set_ylabel("Bx, By, Bz (nT)")
ax4.legend(["Bx", "By", "Bz"])
ax4.grid()

ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
plt.plot(t, np.linalg.norm(B, axis=1))
ax3.grid()
for xc in [t1, t2, t3, t4]:
    plt.axvline(x=xc, color="k", linewidth=1)
plt.axvline(x=ti_MVA, color="r", linewidth=1)
plt.axvline(x=tf_MVA, color="r", linewidth=1)
ax3.set_ylabel("|B| (nT)")
ax3.set_xlabel("Tiempo (hdec)")

ax5 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
ax5.set_ylabel("Energia", picker=True)  # , bbox=dict(facecolor='red'))
plt.setp(ax5.get_xticklabels(), visible=False)
for xc in [t1, t2, t3, t4]:
    plt.axvline(x=xc, color="k", linewidth=1)
plt.axvline(x=ti_MVA, color="r", linewidth=1)
plt.axvline(x=tf_MVA, color="r", linewidth=1)
im = plt.imshow(
    flux_plot,
    aspect="auto",
    origin="lower",
    extent=(t_swea[0], t_swea[-1], energia[-1], energia[0]),
    cmap="inferno",
    norm=LogNorm(vmin=1e4, vmax=1e9),
)
divider = make_axes_locatable(ax5)
cax = divider.append_axes("top", size="7%", pad="1%")
cb = plt.colorbar(im, cax=cax, orientation="horizontal")
cax.xaxis.set_ticks_position("top")


ax7 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
plt.setp(ax7.get_xticklabels(), visible=False)
ax7.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
plt.plot(t_swia, density)
plt.axvline(x=ti_MVA, color="r", linewidth=1)
plt.axvline(x=tf_MVA, color="r", linewidth=1)
ax7.grid()
for xc in [t1, t2, t3, t4]:
    plt.axvline(x=xc, color="k", linewidth=1)

plt.suptitle(f"MAVEN {year}-{doy}")


plt.show(block=False)
