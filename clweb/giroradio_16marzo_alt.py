from sys import exit
import numpy as np
from importar_datos import importar_mag, importar_swia, importar_fila
import matplotlib.pyplot as plt

np.set_printoptions(precision=4)
import sys

sys.path.append("..")
from funciones import fechas, donde

"""
VERSIÓN ALTERNATIVA donde primero recorta y después diezma (no cambia mucho los resultados)
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
    return rg


def giroradio_termico(B, temp):
    """sería mejor haciendo un barrido, pero vamos a empezar por lo fácil"""
    B_avg = np.mean(B, axis=0) * 1e-5  # gauss

    thermal_gyroradius = 1.02e02 * np.sqrt(np.mean(temp)) / np.linalg.norm(B_avg)  # cm
    return thermal_gyroradius * 1e-5  # km


def long_inercial(density):
    density_mean = np.zeros(len(density))
    paso = 20  # cada paso son 4 segundos.
    for i in range(len(density)):
        density_mean[i] = np.mean(density[i : +i + paso])
        # lo ideal sería tomar desde atrás del ti así no se mete en la MPB nunca

    ion_length = 2.28e07 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
    return ion_length


mag, t_mag_entero, B_entero, posicion = importar_mag(year, month, day, ti, tf)

# plt.plot(t_swia, density, label="recorte")
# plt.plot(t_swia, density_full, label="full")
# plt.axvline(x=sw_i, c="k")
# plt.axvline(x=sw_f, c="k")
# plt.axvline(x=ms_f)
# plt.axvline(x=ms_f)
# plt.legend()
# plt.title("SWICA")
# plt.show()

# tomamos los límites de inicio y fin de la MS para ambos instrumentos
def limites(inp):
    if inp == "MS":
        (
            t_swia,
            density_full,
            density,
            temperature_full,
            temperature,
            velocity_full,
            velocity,
        ) = importar_swica(ti, tf)
        inicio = 18.08
        fin = 18.21
    elif inp == "SW":
        (
            t_swia,
            density_full,
            density,
            temperature_full,
            temperature,
            velocity_full,
            velocity,
        ) = importar_swifa(ti, tf)
        inicio = 17.5
        fin = 17.8
    mag_i = donde(t_mag_entero, inicio)
    mag_f = donde(t_mag_entero, fin)
    swia_i = donde(t_swia, inicio)
    swia_f = donde(t_swia, fin)
    # recorto en la MS
    # tengo el B (único) y dos valores de velocidad, dependiendo del recorte de la DF
    # lo mismo para la temperatura y densidad
    B = B_entero[mag_i:mag_f] * 1e-9  # T
    vel_cut = velocity[swia_i:swia_f, :3]
    vel_full_cut = velocity_full[swia_i:swia_f, :3]
    temp_cut = temperature[swia_i:swia_f]
    temp_full_cut = temperature_full[swia_i:swia_f]
    density_cut = density[swia_i:swia_f]
    density_full_cut = density_full[swia_i:swia_f]

    # lo diezmamos para que tengan la misma cantidad de puntos
    B_cut = diezmar(B, density)

    return (
        B_cut,
        vel_cut,
        vel_full_cut,
        temp_cut,
        temp_full_cut,
        density_cut,
        density_full_cut,
    )


# para la magnetofunda
(
    B_cut,
    vel_cut,
    vel_full_cut,
    temp_cut,
    temp_full_cut,
    density_cut,
    density_full_cut,
) = limites("MS")


# ahora sí, calculamos
rg = giroradio(B_cut, vel_cut)
rg_full = giroradio(B_cut, vel_full_cut)

print(f"gr MS 1-5000 eV = {rg_full} km")
print(f"gr MS recorte = {rg} km")

th_gr = giroradio_termico(B_cut, temp_cut)
th_gr_full = giroradio_termico(B_cut, temp_full_cut)
print(f"thermal gr MS 1-5000 eV = {np.nanmean(th_gr_full, axis=0):1.3g} km")
print(f"thermal gr MS recorte = {np.nanmean(th_gr, axis=0):1.3g} km")
# nanmean ignora los nans

proton_length = long_inercial(density_cut)
proton_length_full = long_inercial(density_full_cut)

print(f"long inercial MS 1-5000 eV = {proton_length:1.3g} km")
print(f"long inercial MS recorte = {proton_length_full:1.3g} km")


# para el SW

(
    B_cut,
    vel_cut,
    vel_full_cut,
    temp_cut,
    temp_full_cut,
    density_cut,
    density_full_cut,
) = limites("SW")
# ahora sí, calculamos
rg = giroradio(B_cut, vel_cut)
rg_full = giroradio(B_cut, vel_full_cut)

print(f"gr SW 1-5000 eV = {rg_full} km")
print(f"gr SW recorte = {rg} km")

th_gr = giroradio_termico(B_cut, temp_cut)
th_gr_full = giroradio_termico(B_cut, temp_full_cut)
print(f"thermal gr SW 1-5000 eV = {np.nanmean(th_gr_full, axis=0):1.3g} km")
print(f"thermal gr SW recorte = {np.nanmean(th_gr, axis=0):1.3g} km")
# nanmean ignora los nans

proton_length = long_inercial(density_cut)
proton_length_full = long_inercial(density_full_cut)

print(f"long inercial SW 1-5000 eV = {proton_length:1.3g} km")
print(f"long inercial SW recorte = {proton_length_full:1.3g} km")
