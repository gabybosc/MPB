import numpy as np
import os
import matplotlib.pyplot as plt
from importar_datos import importar_swia
from funciones import find_nearest, fechas, donde, tiempos

np.set_printoptions(precision=4)

"""
Calcula el E convectivo. Igual debería asegurarme de que esté funcionando bien porque tengo mis dudas.
Usa los datos de SWIA de pds y los de MAG de clweb.
"""

# ##########DATOS
year, month, day, doy = fechas()
ti, tf = tiempos("magnetofunda + mpb")


def importar_mag(year, month, day, ti, tf):
    # path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
    path = f"../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la desktop.
    # Estaría bueno ponerle un if para que detecte en cuál estoy.
    if os.path.isfile(path + "mag_filtrado.txt"):
        mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    else:
        mag = np.loadtxt(path + "MAG.asc")
        B = mag[:, 6:9]

    hh = mag[:, 3]
    mm = mag[:, 4]
    ss = mag[:, 5]

    t = hh + mm / 60 + ss / 3600  # hdec

    posicion = np.zeros((len(t), 3))
    for i in range(9, 12):
        posicion[:, i - 9] = mag[:, i]

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    posicion_cut = posicion[inicio:fin]
    return mag, t_cut, B_cut, posicion_cut


swia, t_swia, density, temperature, vel_mso = importar_swia(year, month, day, ti, tf)
mag, t_mag, B, posicion = importar_mag(year, month, day, ti, tf)


# quiero diezmar el tiempo y el campo para que sean los mismos que tiene swia
idx = np.zeros(len(t_swia))
for i in range(len(t_swia)):
    idx[i] = donde(t_mag, t_swia[i])
idx = idx.astype(int)

tmag_diezmado = t_mag[idx]  # lo diezmó

B_cut = B[idx]
posicion_cut = posicion[idx]

ti_up, tf_up = tiempos("región upstream")

inicio_up = donde(t_swia, ti_up)
fin_up = donde(t_swia, tf_up)

inicio_up_mag = donde(tmag_diezmado, ti_up)
fin_up_mag = donde(tmag_diezmado, tf_up)

####################


B_avg = np.empty((len(idx), 3))
v_maven = np.empty((len(idx), 3))
v_planet = np.empty((len(idx), 3))


for i in range(len(idx)):
    B_avg[i, :] = np.mean(B_cut[i : i + 30, :], axis=0)
    v_maven[i, :] = (
        (posicion_cut[inicio_up_mag + 1, :] - posicion_cut[inicio_up_mag, :])
        / (tmag_diezmado[inicio_up_mag + 1] - tmag_diezmado[inicio_up_mag])
        / 3600
    )  # en km/s
for i in range(len(idx)):
    v_planet[i, :] = np.nanmean(
        vel_mso[i : i + 30, :] + v_maven[i : i + 30, :], axis=0
    )  # km/s

E_convective = np.cross(-v_planet * 1e3, B_avg * 1e-9)
E_convective_normalized = np.linalg.norm(E_convective, axis=1)  # V/m

print(f"El Ecv medio es {np.mean(E_convective, axis=0)}")

plt.plot(t_swia, E_convective * 1e3)
plt.legend(["x", "y", "z"])
plt.xlabel("tiempo hdec")
plt.ylabel("Ecv mV/m")

ti_funda = donde(tmag_diezmado, 18.06)
tf_funda = donde(tmag_diezmado, 18.2193)
print(
    "El E convectivo medio en la magnetofunda es {0:1.3g} mV/m".format(
        np.mean(E_convective_normalized[ti_funda:tf_funda] * 1e3)
    )
)  # en V/m


plt.show()
