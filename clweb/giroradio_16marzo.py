from sys import exit
import numpy as np
from importar_datos import importar_mag, importar_swia, importar_fila
import matplotlib.pyplot as plt

np.set_printoptions(precision=4)
import sys

sys.path.append("..")
from funciones import fechas, donde

"""
La long inercial depende de la densidad. Entonces hace tres cosas:
Para el SW calcula usando SWIFA entre 600 y 2000 eV para excluir a los alfas
Para la MS calcula usando SWICA entre 1 y 5000 (con alfa) y después entre 600 y 1400 eV (sin alfa)
En la MS no es tan fácil discriminar alfa de protones entonces no es claro
"""

# #########CONSTANTES
mp = 1.67e-27  # masa del proton en kg
kB = 1.38e-23  # cte de Boltzmann en J/K
q_e = 1.602e-19  # carga del electrón en C

# ##########DATOS
year, month, day, doy = 2016, "03", 16, "076"
ti = 17
tf = 18.5

normal = [0.920, -0.302, 0.251]
t1, t2, t3, t4 = 18.2167, 18.2204, 18.235, 18.2476


def importar_swica(ti, tf):
    path = f"../../../datos/clweb/2016-03-16/momentos/SWICA/"

    dens_full = np.loadtxt(path + "dens_full.asc")
    dens = np.loadtxt(path + "dens_corte.asc")
    temp_full = np.loadtxt(path + "temp_full.asc")
    temp = np.loadtxt(path + "temp_corte.asc")
    vel_full = np.loadtxt(path + "vel_full.asc")
    vel = np.loadtxt(path + "vel_corte.asc")

    t = dens[:, 3] + dens[:, 4] / 60 + dens[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_full = dens_full[inicio:fin, 6]
    density = dens[inicio:fin, 6]
    temperature_full = temp_full[inicio:fin, 6]
    temperature = temp[inicio:fin, 6]
    velocity_full = vel_full[inicio:fin, 6:]
    velocity = vel[inicio:fin, 6:]

    return (
        t_cut,
        density_full,
        density,
        temperature_full,
        temperature,
        velocity_full,
        velocity,
    )


def importar_swifa(ti, tf):
    path = f"../../../datos/clweb/2016-03-16/momentos/SWIFA/"

    dens_full = np.loadtxt(path + "dens_full.asc")
    dens = np.loadtxt(path + "dens_corte.asc")
    temp_full = np.loadtxt(path + "temp_full.asc")
    temp = np.loadtxt(path + "temp_corte.asc")
    vel_full = np.loadtxt(path + "vel_full.asc")
    vel = np.loadtxt(path + "vel_corte.asc")

    t = dens[:, 3] + dens[:, 4] / 60 + dens[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_full = dens_full[inicio:fin, 6]
    density = dens[inicio:fin, 6]
    temperature_full = temp_full[inicio:fin, 6]
    temperature = temp[inicio:fin, 6]
    velocity_full = vel_full[inicio:fin, 6:]
    velocity = vel[inicio:fin, 6:]

    return (
        t_cut,
        density_full,
        density,
        temperature_full,
        temperature,
        velocity_full,
        velocity,
    )


mag, t_mag_entero, B_entero, posicion = importar_mag(year, month, day, ti, tf)
(
    t_swia_entero,
    density_full,
    density,
    temperature_full,
    temperature,
    velocity_full,
    velocity,
) = importar_swica(ti, tf)
# t_swifa_entero, density_full, density, temperature_full, temperature, velocity_full, velocity = importar_swifa(ti, tf)
#  plt.plot(t_swia_entero, density)
# plt.show()

sw_i = 17.5
sw_f = 17.8

ms_i = 18.08
ms_f = 18.21

plt.plot(t_swia_entero, density, label="recorte")
plt.plot(t_swia_entero, density_full, label="full")
plt.axvline(x=sw_i, c="k")
plt.axvline(x=sw_f, c="k")
plt.axvline(x=ms_f)
plt.axvline(x=ms_f)
plt.legend()
plt.show()

"""
giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
where vperp is the component of the velocity perpendicular to the
direction of the magnetic field and B is the strength of the magnetic field
"""

# primero diezmamos
def diezmar(arr_grande, arr_chico):
    """recorta el arr_grande para que tenga el largo del chico"""
    paso = int(len(arr_grande) / len(arr_chico))
    r = arr_grande[::paso]
    return r


def normalizar(vector):
    r = vector / np.linalg.norm(vector)
    return r


def proyeccion(v1, v2):
    """proyecta v1 sobre v2"""
    prod = np.dot(v1, v2)
    N = np.zeros((len(prod), len(v2)))

    for i in range(len(N)):
        N[i, :] = prod[i] * v2
    return N


def giroradio(B, velocidad):
    # ahora que están recortadas las convierto en lo que busco: la intensidad de B y la v perp
    B_medio = np.mean(B, axis=0)
    B_medio_normalizado = normalizar(B_medio)

    # proyecto la velocidad sobre B
    N = proyeccion(velocidad, B_medio_normalizado)

    v_perp = np.mean(velocidad - N, axis=0)

    # el giroradio entonces:
    gf = mp / (q_e * np.linalg.norm(B_medio))  # girofreq
    rg = gf * np.linalg.norm(v_perp)
    rg_min = (
        mp
        * min(np.linalg.norm(velocidad - N, axis=1))
        / (q_e * np.linalg.norm(B_medio))
    )
    rg_max = (
        mp
        * max(np.linalg.norm(velocidad - N, axis=1))
        / (q_e * np.linalg.norm(B_medio))
    )

    return rg


B = diezmar(B_entero, density)
t_mag = diezmar(t_mag_entero, density)
mag_i = donde(t_mag, ti)
mag_f = donde(t_mag, tf)
swia_i = donde(t_swia_entero, ti)
swia_f = donde(t_swia_entero, tf)
# recorto
B_cut = B[mag_i:mag_f] * 1e-9  # T
velocidad = velocity[swia_i:swia_f, :3]


rg = giroradio(t_mag, t_swia_entero, ti, tf, B, velocity)
rg_full = giroradio(t_mag, t_swia_entero, ti, tf, B, velocity_full)

print(f"gr MS 1-5000 eV = {rg_full} km")
print(f"gr MS 600-1400 eV = {rg} km")


def giroradio_termico(B, temp):

    B_avg = np.array(
        [np.mean(B[i : i + 30, :], axis=0) * 1e-5 for i in range(len(B) - 1)]
    )

    thermal_gyroradius = np.zeros(swia_f - swia_i)
    for i in range(swia_f - swia_i):
        thermal_gyroradius[i] = (
            1.02e02
            * np.sqrt(temp[swia_i + i])
            / np.linalg.norm(B_avg[swia_i + i, :])
            * 1e-5
        )  # km


th_gr = giroradio_termico()
print(
    f"thermal ion gyroradius mean = {np.nanmean(th_gr, axis=0):1.3g} km"
)  # nanmean ignora los nans


"""
Longitud inercial
"""


def long_inercial(density_full, swia_i, swia_f):
    density_mean = np.zeros(swia_f - swia_i)  # upstream
    paso = 20  # cada paso son 4 segundos.
    if swia_i - paso > 0:  # si no se cumple, va a calcularlo mal
        for i in range(swia_f - swia_i):
            density_mean[i] = np.mean(
                density_full[swia_i + i - paso : swia_i + i]
            )  # toma desde atrás del ti así no se mete en la MPB nunca

        ion_length = 2.28e07 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
        # ion_min = 2.28e07 / np.sqrt(max(density_mean)) * 1e-5  # km
        # ion_max = 2.28e07 / np.sqrt(min(density_mean)) * 1e-5  # km
        return ion_length
