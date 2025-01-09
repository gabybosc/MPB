import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.widgets import MultiCursor
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from _update_parametros import hoja_MVA_update, hoja_MVA_analisis, hoja_param

sys.path.append("..")
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
    UTC_to_hdec,
    SZA,
    altitude,
)

"""
tdec, Bx, By, Bz, modulo B,pos x, pos y, pos z, distancia km
"""

"""
Uso los datos de 1s para la posición porque están en TSWIS
"""
# path = "../../../datos/Titan/t96_tswis_1s.ascii"
# datos = np.loadtxt(path)
# tiempo = datos[:, 0]
# i = donde(tiempo, 23)
# f = donde(tiempo, 26)
# t_low, posicion = datos[i:f, 0], datos[i:f, 5:8]
"""
Uso los datos de alta resolución para el MVA, pero están en KSO
"""
# path_t_hires = "../../../datos/Titan/t96_kso_hires.txt"  # solo para el tiempo
# path_B_hires = "../../../datos/Titan/t96_kso_hires_filt.gz"  # filtrado, para B
# tiempo = np.genfromtxt(path_t_hires, skip_header=1, dtype="str", usecols=[1])
# t_hires = np.array([UTC_to_hdec(x) for x in tiempo])
#
# B = np.loadtxt(path_B_hires)
# Bnorm = np.linalg.norm(B, axis=1)
# B_para, B_perp_norm, t_plot = Bpara_Bperp(
#     B, t_hires, t_hires[0] + 0.2, t_hires[-1] - 0.2
# )

# np.save("../../../datos/Titan/t_low.npy", t_low)
# np.save("../../../datos/Titan/posicion.npy", posicion)
# np.save("../../../datos/Titan/tiempo.npy", tiempo)
# np.save("../../../datos/Titan/t_hires.npy", t_hires)
# np.save("../../../datos/Titan/B.npy", B)
# np.save("../../../datos/Titan/Bnorm.npy", Bnorm)
# np.save("../../../datos/Titan/B_para.npy", B_para)
# np.save("../../../datos/Titan/B_perp_norm.npy", B_perp_norm)
# np.save("../../../datos/Titan/t_plot.npy", t_plot)
t_low = np.load("../../../datos/Titan/t_low.npy")
tiempo = np.load("../../../datos/Titan/tiempo.npy")
t_hires = np.load("../../../datos/Titan/t_hires.npy")
t_plot = np.load("../../../datos/Titan/t_plot.npy")
posicion = np.load("../../../datos/Titan/posicion.npy")
B = np.load("../../../datos/Titan/B.npy")
Bnorm = np.load("../../../datos/Titan/Bnorm.npy")
B_para = np.load("../../../datos/Titan/B_para.npy")
B_perp_norm = np.load("../../../datos/Titan/B_perp_norm.npy")

t_mva_i = "00:33:05"  # 33:01
t_mva_f = "00:33:26"  # 33:42

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
t1, t2, t3, t4 = 0.54175182, 0.55058123, 0.57651763, 0.58203602


def importar_gdocs():
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive.file",
        "https://www.googleapis.com/auth/drive",
    ]

    creds = ServiceAccountCredentials.from_json_keyfile_name("../mpb_api.json", scope)

    client = gspread.authorize(creds)

    hoja_parametros = client.open("MPB").worksheet("Parametros")
    hoja_MVA = client.open("MPB").worksheet("MVA")

    return hoja_parametros, hoja_MVA


def MVA(t, B, posicion):
    # Calcular la matriz de covarianza de B
    M_ij = np.cov(B.T)
    # avec, lamb = autovectores(M_ij)
    lamb, x = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

    # Los ordeno de mayor a menor
    idx = lamb.argsort()[::-1]
    lamb = lamb[idx]
    x = x[:, idx]
    # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
    x1 = x[:, 0]
    x2 = x[:, 1]
    x3 = x[:, 2]

    if x3[0] > 0:  # si la normal aputna para afuera me la da vuelta xq estoy en TSWIS
        x3 = -x3
    if any(np.cross(x1, x2) - x3) > 0.01:
        # print("Cambio el signo de x1 para que los av formen terna derecha")
        print(x1, np.cross(x1, x2) - x3)
        x1 = -x1
        print(x1)

    avec = [x1, x2, x3]

    print("La normal del MVA es:", avec[2])

    # Calcular las proyecciones
    B1, B2, B3 = [np.dot(B, avec[i]) for i in range(3)]

    # Calcular el campo medio y la altitud
    B_medio_vectorial = np.mean(B, axis=0)
    altitud = np.linalg.norm(posicion, axis=1) - 2570  # km
    altitud_media = np.mean(altitud)

    # Ángulo solar cenital y altitud
    print(
        f"Altitud media = {altitud_media} km, SZA = {SZA(posicion, len(posicion) // 2) - 90}°"
    )

    # Mostrar eigenvalores y cociente
    print(f"l1 = {lamb[0]}, l2 = {lamb[1]}, l3 = {lamb[2]}")
    print("Cociente de lambdas:", lamb[1] / lamb[2])

    # Calcular la norma de B y la componente normal
    B_norm_medio = np.linalg.norm(B_medio_vectorial)
    print(
        f"El campo magnético medio es {B_medio_vectorial}, su norma es {B_norm_medio:.2f}"
    )
    print(f"La componente normal media de B es {np.mean(B3):.2f}")

    # Cociente entre la componente normal y la norma del campo
    print(f"Cociente <B_3>/|<B>| es {np.mean(B3) / B_norm_medio:.2f}")

    # Graficar el hodograma
    hodograma(B1, B2, B3)

    # Calcular el ángulo entre el campo medio y la normal
    angulo_B = angulo(B_medio_vectorial, avec[2]) * 180 / np.pi
    print(
        f"El ángulo entre el vector de campo magnético medio y la normal es {angulo_B:.2f}°"
    )

    # Calcular el error
    phi, delta_B3 = error(lamb, B, avec[2])
    if phi[2, 1] > phi[2, 0]:
        error_normal = phi[2, 1] * 180 / np.pi
    else:
        error_normal = phi[2, 0] * 180 / np.pi
    hoja_parametros, hoja_mva = importar_gdocs()
    nr = 70
    hoja_MVA_update(
        hoja_mva, nr, t, lamb, avec, error_normal, B3, delta_B3, B_norm_medio, angulo_B
    )
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


# 30-nov-2013: central: 24:33:15, radio = 10s

# ti = donde(t, UTC_to_hdec("24:33:05"))
# tf = donde(t, UTC_to_hdec("24:33:25"))
ti, ti_hires = donde(t_low, 24 + UTC_to_hdec(t_mva_i)), donde(
    t_hires, UTC_to_hdec(t_mva_i)
)
tf, tf_hires = donde(t_low, 24 + UTC_to_hdec(t_mva_f)), donde(
    t_hires, UTC_to_hdec(t_mva_f)
)

alt_i = np.linalg.norm(posicion[donde(t_low, 24 + UTC_to_hdec(t_mva_i)), :]) - 2570
alt_f = np.linalg.norm(posicion[donde(t_low, 24 + UTC_to_hdec(t_mva_f)), :]) - 2570
print(f"Altitud al inicio y fin del MVA: inicio = {alt_i}, fin = {alt_f}")
print("SZA inicio", SZA(posicion, donde(t_low, 24 + UTC_to_hdec(t_mva_i))) - 90)
print("SZA fin", SZA(posicion, donde(t_low, 24 + UTC_to_hdec(t_mva_f))) - 90)

x3, B_cut, t_cut, posicion_cut = MVA(
    t_hires[ti_hires:tf_hires], B[ti_hires:tf_hires], posicion[ti:tf]
)

B_upstream = corte(t_hires, t1 - 0.015, t1, B)
B_downstream = corte(t_hires, t4, t4 + 0.015, B)

omega = angulo(B_upstream, B_downstream)
print(f"omega = {omega * 180 / np.pi}º")

v_punto = np.zeros((len(posicion_cut) - 1, 3))
norma_v = np.zeros(len(posicion_cut) - 1)
for i in range(len(v_punto)):
    v_punto[i, :] = posicion_cut[i + 1, :] - posicion_cut[i]
    # en km/s, tiene resolución de 1Hz
    norma_v[i] = np.linalg.norm(v_punto[i, :])

# la velocidad promedio
v_media = np.mean(v_punto, axis=0)
angulo_v = angulo(v_media, x3) * 180 / np.pi
print("v media", v_media)
print(f"El ángulo entre v media y la normal es {angulo_v}")

x14, x23 = ancho_mpb(t1, t2, t3, t4, x3, v_media)
print(f"El ancho x14 es {x14}km y el ancho x23={x23} km")

J_s, J_v = corrientes(x3, B_upstream, B_downstream, x23)

print(
    f"J sup = {J_s}, |Js| = {np.linalg.norm(J_s)}, J vol {J_v}, |Jv| = {np.linalg.norm(J_v)}"
)
inicio_down = donde(t_hires, t1 - 0.015)
# fuerza_mva = fuerza(J_v_MVA, B[inicio_down, :])

# altitude(posicion, 2574)
print("SZA para t1, t2, t3, t4:\n")
for ts in [t1, t2, t3, t4]:
    print(ts, SZA(posicion, donde(t_low, ts)))

hoja_parametros, hoja_mva = importar_gdocs()
nr = 70

hoja_param(
    hoja_parametros,
    nr,
    [t1, t2, t3, t4],
    SZA(posicion, donde(t_low, t2)),
    v_media,
    omega,
    B_upstream,
    B_downstream,
)
hoja_MVA_analisis(hoja_mva, nr, x14, x23, J_s, J_v, fuerza=0)
