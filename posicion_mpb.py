import numpy as np
import matplotlib.pyplot as plt
from funciones import Bpara_Bperp, fechas, hdec_to_UTC, UTC_to_hdec, donde
from importar_datos import importar_mag_1s
import pickle


"""
Grafica la distribución de los cruces de MPB y BS del grupo 4 (por ahora)
"""

np.set_printoptions(precision=4)


def crear_txt_posiciones():
    "esto lo corro una sola vez y después ya tengo el txt con las posiciones"
    catalogo = np.genfromtxt("outputs/grupo4.txt", dtype="str")

    lst = []
    for cat in catalogo:
        year, month, day = cat[0].split("-")

        t_bs = UTC_to_hdec(cat[1])
        t_mpb = UTC_to_hdec(cat[2])

        if t_bs < t_mpb:
            ti = t_bs - 1
            tf = t_mpb + 1
        else:
            ti = t_mpb - 1
            tf = t_bs + 1
        if ti < 0:
            ti = 0
        if tf > 24:
            tf = 24
        mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
        pos_mpb = posicion[donde(t, t_mpb)]
        pos_bs = posicion[donde(t, t_bs)]
        lst.append([year, month, day, t_bs, pos_bs, t_mpb, pos_mpb])

    with open("outputs/etiquetas_grupo4", "wb") as fp:
        pickle.dump(lst, fp)


def marte():
    fig, ax = plt.subplots()
    ax.plot()
    ax.axis("equal")
    ax.set_xlim(-3, 3)
    ax.set_ylim(0, 3)
    circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
    ax.add_artist(circle)
    ax.set_title("MAVEN MSO coordinates", fontsize=16)
    ax.set_xlabel(r"$X_{MSO}$ ($R_M$)", fontsize=14)
    ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)", fontsize=14)


def frontera(pos):
    x = np.zeros(len(pos))
    yz = np.zeros(len(pos))
    for i in range(len(pos)):
        x[i] = pos[i][0] / 3390
        yz[i] = np.sqrt(pos[i][1] ** 2 + pos[i][2] ** 2) / 3390

    return x, yz


# crear_txt_posiciones()

with open("outputs/etiquetas_grupo4", "rb") as fp:
    cruces = pickle.load(fp)


pos_bs = [cruces[i][4] for i in range(len(cruces))]
pos_mpb = [cruces[i][6] for i in range(len(cruces))]

marte()
x_bs, yz_bs = frontera(pos_bs)
x_mpb, yz_mpb = frontera(pos_mpb)
plt.scatter(x_bs, yz_bs)
plt.scatter(x_mpb, yz_mpb)
plt.show()
