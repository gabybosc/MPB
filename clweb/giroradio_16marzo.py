from sys import exit
import numpy as np
from importar_datos import importar_mag, importar_swia, importar_fila

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


def importar_swia(ti, tf):
    path = f"../../../datos/clweb/2016-03-16/"

    swica_full = np.loadtxt(path + "densidad_swica_1-5000.asc")
    swica = np.loadtxt(path + "densidad_swica_600-1400.asc")
    swifa = np.loadtxt(path + "densidad_swifa_600-2000.asc")

    t = swica[:, 3] + swica[:, 4] / 60 + swica[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_full = swica_full[inicio:fin, 6]
    density_swica = swica[inicio:fin, 6]
    density_swifa = swifa[inicio:fin, 6]

    return t_cut, density_full, density_swica, density_swifa


# def importar_temp(ti, tf):
#     path = "../../datos/temp_swica_ms.asc"
#     swica = np.loadtxt(path)

#     t = swica[:, 3] + swica[:, 4] / 60 + swica[:, 5] / 3600

#     inicio = donde(t, ti)
#     fin = donde(t, tf)

#     t_cut = t[inicio:fin]
#     temp = swica[inicio:fin, -1]

#     return t_cut, temp


mag, t_mag_entero, B_entero, posicion = importar_mag(year, month, day, ti, tf)
t_swia_entero, density_full, density_swica, density_swifa = importar_swia(ti, tf)
# plt.plot(t_swia_entero, density)
# plt.show()

sw_i = 17.5
sw_f = 17.8

ms_i = 18.08
ms_f = 18.21

"""
giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
where vperp is the component of the velocity perpendicular to the
direction of the magnetic field and B is the strength of the magnetic field
"""

# primero diezmamos
paso = int(len(B_entero) / len(density_swica))
B = B_entero[::paso]
t_mag = t_mag_entero[::paso]


for k in [(sw_i, sw_f), (ms_i, ms_f)]:

    # selecciono las regiones:
    mag_i = donde(t_mag, k[0])
    mag_f = donde(t_mag, k[1])
    swia_i = donde(t_swia_entero, k[0])
    swia_f = donde(t_swia_entero, k[1])
    # tt, temperature = importar_temp(ms_i, ms_f)

    # # recorto
    # B_cut = B[mag_i:mag_f] * 1e-9  # T
    # velocidad = vel_mso_xyz[swia_i:swia_f]

    # # ahora que están recortadas las convierto en lo que busco: la intensidad de B y la v perp
    # B_medio = np.mean(B_cut, axis=0)
    # B_medio_normalizado = B_medio / np.linalg.norm(B_medio)

    # # proyecto la velocidad sobre B
    # dot = np.dot(velocidad, B_medio_normalizado)
    # N = np.zeros((len(dot), len(B_medio)))

    # for i in range(len(N)):
    #     N[i, :] = dot[i] * B_medio_normalizado

    # v_perp = np.mean(velocidad - N, axis=0)

    # # el giroradio entonces:
    # gf = mp / (q_e * np.linalg.norm(B_medio))  # girofreq
    # rg = gf * np.linalg.norm(v_perp)
    # rg_min = (
    #     mp
    #     * min(np.linalg.norm(velocidad - N, axis=1))
    #     / (q_e * np.linalg.norm(B_medio))
    # )
    # rg_max = (
    #     mp
    #     * max(np.linalg.norm(velocidad - N, axis=1))
    #     / (q_e * np.linalg.norm(B_medio))
    # )

    # B_avg = np.empty((len(B_cut), 3))
    # B_avg_normalized = np.empty((len(B_cut), 3))
    # temp_para_xyz = np.empty((len(B_cut), 3))

    # for i in range(len(B_cut) - 1):
    #     B_avg[i, :] = np.mean(B_cut[i : i + 30, :], axis=0) * 1e-5  # lo paso a gauss
    #     B_avg_normalized[i, :] = B_avg[i, :] / np.linalg.norm(
    #         B_avg[i, :]
    #     )  # adimensional
    #     temp_para_xyz[i, :] = (
    #         np.dot(B_avg_normalized[i, :], temperature[i, :]) * B_avg_normalized[i, :]
    #     )  # eV

    # temp_para = np.linalg.norm(temp_para_xyz, axis=1)  # eV
    # temp_perp = np.linalg.norm(temperature - temp_para_xyz, axis=1)  # eV

    # thermal_gyroradius = np.empty(swia_f - swia_i)
    # for i in range(swia_f - swia_i):
    #     thermal_gyroradius[i] = (
    #         1.02e02
    #         * np.sqrt(temp_perp[swia_i + i])
    #         / np.linalg.norm(B_avg[swia_i + i, :])
    #         * 1e-5
    #     )  # km

    # print(
    #     f"thermal ion gyroradius mean = {np.nanmean(thermal_gyroradius, axis=0):1.3g} km"
    # )  # nanmean ignora los nans

    # if k[0] == sw_i:
    #     print(
    #         f"El radio de Larmor de iones SW es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
    #     )
    # if k[0] == ms_i:
    #     print(
    #         f"El radio de Larmor de iones MS es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
    #     )

    # if k[0] == sw_i:
    #     print(f"La girofrecuencia de iones SW es {gf:1.3g} s⁻¹")
    # if k[0] == ms_i:
    #     print(f"La girofrecuencia de iones MS es {gf:1.3g} s⁻¹")

    # """
    # Proyectando en la normal:
    # """
    # # normal = np.array([0.920,-0.302,0.251])
    # v_normal = np.dot(np.mean(velocidad, axis=0), normal)

    # # el giroradio entonces:
    # rg_normal = mp * np.linalg.norm(v_normal) / (q_e * np.linalg.norm(B_medio))

    # if k[0] == sw_i:
    #     print(
    #         f"El radio de Larmor de iones SW con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
    #     )
    # if k[0] == ms_i:
    #     print(
    #         f"El radio de Larmor de iones MS con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
    #     )

    """
    Longitud inercial
    """
    density_mean = np.zeros(swia_f - swia_i)  # upstream
    paso = 20  # cada paso son 4 segundos.
    if k[0] == sw_i:
        if swia_i - paso > 0:  # si no se cumple, va a calcularlo mal
            for i in range(swia_f - swia_i):
                density_mean[i] = np.mean(
                    density_swifa[swia_i + i - paso : swia_i + i]
                )  # toma desde atrás del ti así no se mete en la MPB nunca

            ion_length = 2.28e07 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
            ion_min = 2.28e07 / np.sqrt(max(density_mean)) * 1e-5  # km
            ion_max = 2.28e07 / np.sqrt(min(density_mean)) * 1e-5  # km
        print(
            f"La longitud inercial de protones (sin alfa) SW es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g})"
        )

    if k[0] == ms_i:
        if swia_i - paso > 0:  # si no se cumple, va a calcularlo mal
            for i in range(swia_f - swia_i):
                density_mean[i] = np.mean(
                    density_full[swia_i + i - paso : swia_i + i]
                )  # toma desde atrás del ti así no se mete en la MPB nunca

            ion_length = 2.28e07 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
            ion_min = 2.28e07 / np.sqrt(max(density_mean)) * 1e-5  # km
            ion_max = 2.28e07 / np.sqrt(min(density_mean)) * 1e-5  # km
        print(
            f"La longitud inercial de iones (protones + alfa) MS es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g})"
        )
        if swia_i - paso > 0:  # si no se cumple, va a calcularlo mal
            for i in range(swia_f - swia_i):
                density_mean[i] = np.mean(
                    density_swica[swia_i + i - paso : swia_i + i]
                )  # toma desde atrás del ti así no se mete en la MPB nunca

            ion_length = 2.28e07 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
            ion_min = 2.28e07 / np.sqrt(max(density_mean)) * 1e-5  # km
            ion_max = 2.28e07 / np.sqrt(min(density_mean)) * 1e-5  # km
        print(
            f"La longitud inercial de protones (sin alfa) MS es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g})"
        )
