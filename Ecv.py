import numpy as np
import os
import matplotlib.pyplot as plt
from importar_datos import importar_swia
from funciones import fechas, donde, tiempos, diezmar

np.set_printoptions(precision=4)

"""
Calcula el E convectivo. Igual debería asegurarme de que esté funcionando bien
porque tengo mis dudas.
Usa los datos de SWIA de pds y los de MAG de clweb.
"""

# ##########DATOS
year, month, day, doy = fechas()
ti_ms, tf_ms = 18.06, 18.2193  # tiempos("magnetofunda")
ti_up, tf_up = 18.3, 18.5  # tiempos("región upstream")


def importar_mag(year, month, day, ti, tf):
    path = (
        f"../../datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la laptop
    )
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

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    posicion_cut = posicion[inicio:fin]
    return mag, t_cut, B_cut, posicion_cut


swia, t_swia, density, temperature, vel_mso = importar_swia(
    year, month, day, ti_ms, tf_up
)
mag, t_mag, B, posicion = importar_mag(year, month, day, ti_ms, tf_up)


# quiero diezmar el tiempo y el campo para que sean los mismos que tiene swia
idx = diezmar(t_mag, t_swia)

tmag_diezmado = t_mag[idx]
B_cut = B[idx]
posicion_cut = posicion[idx]

inicio_up = donde(t_swia, ti_up)  # tiene que dar lo mismo si uso tswia o tmag
fin_up = donde(t_swia, tf_up)

ti_funda = donde(tmag_diezmado, ti_ms)
tf_funda = donde(tmag_diezmado, tf_ms)
####################


B_avg = np.empty((len(idx), 3))
v_maven = np.empty((len(idx), 3))
v_planet = np.empty((len(idx), 3))


for i in range(len(idx)):
    B_avg[i, :] = np.mean(B_cut[i : i + 30, :], axis=0)
    v_maven[i, :] = (
        (posicion_cut[inicio_up + 1, :] - posicion_cut[inicio_up, :])
        / (tmag_diezmado[inicio_up + 1] - tmag_diezmado[inicio_up])
        / 3600
    )  # en km/s
for i in range(len(idx)):
    v_planet[i, :] = np.nanmean(
        vel_mso[i : i + 30, :] + v_maven[i : i + 30, :], axis=0
    )  # km/s

E_convective = np.cross(-v_planet * 1e3, B_avg * 1e-9)
E_convective_normalized = np.linalg.norm(E_convective, axis=1)  # V/m

print(f"El Ecv medio es {np.mean(E_convective, axis=0)}")
print(
    "El E convectivo medio en la magnetofunda es {0:1.3g} mV/m".format(
        np.mean(E_convective_normalized[ti_funda:tf_funda] * 1e3)
    )
)  # en mV/m

plt.plot(t_swia, E_convective * 1e3)
plt.plot(t_swia, np.linalg.norm(E_convective * 1e3, axis=1))
plt.axvspan(18.21, 18.25, color="red", alpha=0.2)
plt.legend(["x", "y", "z", "|B|", "MPB"])
plt.xlabel("tiempo hdec")
plt.ylabel("Ecv mV/m")

plt.show()
