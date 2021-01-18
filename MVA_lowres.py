import numpy as np
from funciones import (
    error,
    find_nearest,
    find_nearest_final,
    find_nearest_inicial,
    Mij,
    fechas,
    donde,
)
from funciones_MVA import (
    ajuste_conico,
    plot_bootstrap,
    bootstrap,
)
from importar_datos import importar_t1t2t3t4, importar_mag_1s, importar_lpw


"""
Para datos de MAg de baja resolución

Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o dia_del_año-año) y los cuatro tiempos t1 t2 t3 t4.
Eventualmente podría simplemente encontrar todos los cruces que quiero y decirle que lea directamente de algún lugar eso.
Antes de correr este programa hay que haber usado plot_seleccionar_puntos, para tener los cuatro tiempos elegidos.
Nos devuelve los lambdas, el cociente de lambdas, el omega, los autovectores y la normal para el MVA, el ajuste y el bootstrap.
Nos da el valor medio de B, de la altitud y el SZA. Grafica el hodograma.

Guarda los datos en una spreadsheet de google
"""


np.set_printoptions(precision=4)

year, month, day, doy = fechas()
hour = int(input("Hora en HH"))
t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day, doy, hour)

mag, t, B, posicion = importar_mag_1s(year, month, day, t1, t4)
lpw, t_lpw, e_density = importar_lpw(year, month, day, t1, t4)

t_1 = donde(t, t1)
t_2 = donde(t, t2)
t_3 = donde(t, t3)
t_4 = donde(t, t4)

#################

# ahora empieza el MVA con los datos que elegí
B_cut = B[t_2 : t_3 + 1, :]
posicion_cut = posicion[t_2 : t_3 + 1, :]

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
altitud = np.mean(np.linalg.norm(posicion_cut, axis=1) - 3390)
SZA = (
    np.arccos(
        np.clip(
            np.dot(posicion_cut[0, :] / np.linalg.norm(posicion_cut[0, :]), [1, 0, 0]),
            -1.0,
            1.0,
        )
    )
    * 180
    / np.pi
)

B_norm_medio = np.linalg.norm(B_medio_vectorial)

# hodograma(B1, B2, B3, 'nT', 'MAVEN MAG MVA ')

# el error
phi, delta_B3 = error(lamb, B_cut, len(B_cut), x)

###############
####fit
orbita = (
    posicion[
        np.where(t == find_nearest_inicial(t, t1 - 1))[0][0] : np.where(
            t == find_nearest_final(t, t4 + 1)
        )[0][0],
        :,
    ]
    / 3390
)  # radios marcianos

t_nave = find_nearest(t, (t2 + t3) / 2)  # el tiempo en el medio de la hoja de corriente
index = np.where(t == t_nave)[0][0]
x0 = 0.78
e = 0.9
normal_fit, X1, Y1, Z1, R, L0 = ajuste_conico(posicion, index, orbita, x3)

B3_fit = np.dot(B_cut, normal_fit)

#############
##Bootstrap
N_boot = 1000
normal_boot, phi_boot, delta_B3_boot, out, out_phi = bootstrap(
    N_boot, B_cut, len(B_cut)
)

muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

B3_boot = np.dot(B_cut, normal_boot)

##########
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
