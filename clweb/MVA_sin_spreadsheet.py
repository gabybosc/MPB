import numpy as np
import sys
from importar_datos import importar_mag

sys.path.append("..")
from funciones import error, find_nearest, Mij, angulo, autovectores, donde
from funciones_metodos import normal_fit
from funciones_plot import hodograma

"""

Para datos de mag del clweb.

Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o doy-año)
y los tiempos entre los que voy a realizar el MVA.
Eventualmente podría simplemente encontrar todos los cruces que quiero
y decirle que lea directamente de algún lugar eso. (es lo que hace MVA_automatico)
Antes de correr este programa hay que haber usado plot_seleccionar_puntos, para tener los cuatro tiempos elegidos.
Toma los datos de alta resolución y les aplica un filtro pasabajos con
ventana Butterworth con frecuencia de corte de 0.1 Hz de orden 3.
A los datos filtrados les aplica el MVA.
Nos devuelve los lambdas, el cociente de lambdas, el omega, los autovectores y la normal.
Nos da el valor medio de B, de la altitud y el SZA.
Devuelve el ancho de la MPB y la corriente que pasa. Calcula también la fuerza de lorentz y el campo de hall.
Grafica el hodograma, el ajuste de vignes, y la comparación de las normales obtenidas por los distintos métodos.
"""

np.set_printoptions(precision=4)


def MVA(year, month, day, ti_MVA, tf_MVA):

    mag, t, B, posicion = importar_mag(year, month, day, ti_MVA, tf_MVA)
    # ya los importa cortados a los datos, entonces no hace falta que haga el cut yo

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
    phi, delta_B3 = error(lamb, B, x)
    print("MVA terminado")
    return avec[2], B, t, posicion


def ajuste(year, month, day, doy, ti_MVA, tf_MVA):

    datos_tiempo = np.loadtxt("../outputs/t1t2t3t4.txt")
    idx_d = np.where(int(doy) == datos_tiempo[:, 1].astype(int))[0]
    idx_h = np.where(int(ti_MVA) == datos_tiempo[:, 2].astype(int))[0]
    idx = np.intersect1d(idx_d, idx_h)[0]
    t1 = datos_tiempo[idx, 2]
    t2 = datos_tiempo[idx, 3]
    t3 = datos_tiempo[idx, 4]
    t4 = datos_tiempo[idx, 5]

    mag, t, B, posicion = importar_mag(year, month, day, ti_MVA, tf_MVA)

    t_nave = find_nearest(t, (t2 + t3) / 2)
    # el tiempo en el medio de la hoja de corriente
    index = donde(t, t_nave)
    # x0 = 0.78
    # e = 0.9
    normal_ajuste, L0 = normal_fit(posicion, index)

    B3_fit = np.dot(B, normal_ajuste)

    print(f"La normal del fit es {normal_ajuste}")
    # print(f"B3 del fit es {B3_fit}")
    print("Ajuste terminado")

    return normal_ajuste, t1, t2, t3, t4


def normal_coplanar(B_upstream, B_downstream):
    # Anda mal cuando el ángulo entre B y la normal es 0 o 90.
    deltaB = B_downstream - B_upstream

    coplanar = np.cross(np.cross(B_downstream, B_upstream), deltaB)

    normal = coplanar / np.linalg.norm(coplanar, axis=0)

    return normal
