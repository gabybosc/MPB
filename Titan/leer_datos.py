import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.widgets import MultiCursor
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.append("../..")
from funciones_plot import onpick1, hodograma
from funciones import (
    donde,
    Bpara_Bperp,
    Mij,
    autovectores,
    error,
    find_nearest,
    angulo,
    ancho_mpb,
    corrientes,
)

"""
tdec, Bx, By, Bz, modulo B,pos x, pos y, pos z, distancia km
"""

path = "../../../datos/Titan/t96_tswis_1s.ascii"
datos = np.loadtxt(path)
tiempo = datos[:, 0]
i = donde(tiempo, 22)
f = donde(tiempo, 26)

t, B, Bnorm, posicion = datos[i:f, 0], datos[i:f, 1:4], datos[i:f, 4], datos[i:f, 5:8]
B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, t[0] + 0.2, t[-1] - 0.2)

# happy = False
# while not happy:
#     val = []
#     while len(val) < 4:
#         plt.clf()  # clear figure
#         fig = plt.figure(
#             1, constrained_layout=True
#         )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
#         fig.subplots_adjust(
#             top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
#         )
#         plt.title("Spacebar when ready to click:")
#
#         ax1 = plt.subplot2grid((3, 1), (0, 0))
#         plt.plot(t_plot, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
#         plt.plot(
#             t_plot, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / B"
#         )
#         plt.setp(ax1.get_xticklabels(), visible=False)
#         ax1.set_ylabel(r"|$\Delta B$|/ B")
#         ax1.grid()
#         ax1.legend()
#
#         ax4 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
#         ax4.plot(t, B)
#         plt.setp(ax4.get_xticklabels(), visible=False)
#         ax4.set_ylabel("Bx, By, Bz (nT)")
#         ax4.legend(["Bx", "By", "Bz"])
#         ax4.grid()
#
#         ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
#         plt.plot(t, Bnorm)
#         ax3.grid()
#         ax3.set_ylabel("|B| (nT)")
#         ax3.set_xlabel("Tiempo (hdec)")
#         fig.canvas.mpl_connect("pick_event", onpick1)
#         multi = MultiCursor(fig.canvas, (ax1, ax3, ax4), color="black", lw=1)
#
#         zoom_ok = False
#         print("\nSpacebar when ready to click:\n")
#         while not zoom_ok:
#             zoom_ok = plt.waitforbuttonpress(-1)
#         print("Click to select MPB: ")
#         val = np.asarray(plt.ginput(4))[:, 0]
#         print("Selected values: ", val)
#         outs = sorted(val)
#
#     print("Happy? Keyboard click for yes, mouse click for no.")
#     happy = plt.waitforbuttonpress()
t1, t2, t3, t4 = 24.54175182, 24.55058123, 24.57651763, 24.58203602


def MVA(t, B, posicion):
    M = len(t)

    n_p = int(M / 2)

    M_ij = Mij(B)

    avec, lamb = autovectores(M_ij)

    print("la normal del MVA es ", avec[2])

    # las proyecciones
    B1 = np.dot(B, avec[0])
    B2 = np.dot(B, avec[1])
    B3 = np.dot(B, avec[2])

    # el B medio
    B_medio_vectorial = np.mean(B, axis=0)
    altitud = np.linalg.norm(posicion, axis=1) - 3390  # km
    altitud_media = np.mean(altitud)

    SZA = angulo(posicion[n_p, :], [1, 0, 0]) * 180 / np.pi
    print(f"altitud = {altitud_media}, SZA = {SZA}")

    print("cociente de lambdas = ", lamb[1] / lamb[2])
    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    print(f"El B medio es {B_norm_medio}")
    hodograma(B1, B2, B3)

    # el error
    phi, delta_B3 = error(lamb, B, avec[2])
    print("MVA terminado")
    return avec[2], B, t, posicion


def corte(t, ti, tf, vector):
    inicio = donde(t, ti)
    fin = donde(t, tf)
    if vector.ndim > 1:
        vector_cut = np.nanmean(vector[inicio:fin, :], axis=0)
    else:
        vector_cut = np.nanmean(vector[inicio:fin])
    return vector_cut


ti = donde(t, 24.5516666666)
tf = donde(t, 24.56111111)
x3, B_cut, t_cut, posicion_cut = MVA(t[ti:tf], B[ti:tf], posicion[ti:tf])

B_upstream = corte(t, t1 - 0.015, t1, B)
B_downstream = corte(t, t4, t4 + 0.015, B)

omega = angulo(B_upstream, B_downstream)

v_punto = np.zeros((len(posicion_cut) - 1, 3))
norma_v = np.zeros(len(posicion_cut) - 1)
for i in range(len(v_punto)):
    v_punto[i, :] = posicion_cut[i + 1, :] - posicion_cut[i]
    # en km/s, tiene resolución de 1Hz
    norma_v[i] = np.linalg.norm(v_punto[i, :])

# la velocidad promedio
v_media = np.mean(v_punto, axis=0)
print("v media", v_media)

x14, x23 = ancho_mpb(t1, t2, t3, t4, x3, v_media)

J_s_MVA, J_v_MVA = corrientes(x3, B_upstream, B_downstream, x23)

print(f"J sup = {np.linalg.norm(J_s_MVA)}, J vol = {np.linalg.norm(J_v_MVA)}")
inicio_down = donde(t, t1 - 0.015)
# fuerza_mva = fuerza(J_v_MVA, B[inicio_down, :])
