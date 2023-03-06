from sys import exit
import numpy as np
from funciones import fechas, donde
from importar_datos import importar_mag_1s, importar_swea, importar_fila

np.set_printoptions(precision=4)

"""
Calcula el giroradio y la long inercial de electrones. El giroradio lo calcula usando
v_perp al campo B y también usando la velocidad proyectada en la normal.
"""

# #########CONSTANTES
me = 9.1e-31  # masa del electron en kg
kB = 1.38e-23  # cte de Boltzmann en J/K
q_e = 1.602e-19  # carga del electrón en C

# ##########DATOS
year, month, day, doy = 2016, "03", 16, "076"
ti = 17
tf = 19

normal = [0.920, -0.302, 0.251]
t1, t2, t3, t4 = 18.2167, 18.2204, 18.235, 18.2476

mag, t_mag_entero, B_entero, posicion = importar_mag_1s(year, month, day, ti, tf)


def importar_lpw(ti, tf):
    path = "../../datos/lpw_density.asc"
    lpw = np.loadtxt(path)

    t = lpw[:, 3] + lpw[:, 4] / 60 + lpw[:, 5] / 3600

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    e_density_cut = lpw[inicio:fin, 6]

    return lpw, t_cut, e_density_cut


def importar_swea(ti, tf):
    path_den = "../../datos/swea_density.asc"
    path_vel = "../../datos/swea_vel.asc"
    den = np.loadtxt(path_den)
    vel = np.loadtxt(path_vel)

    t = den[:, 3] + den[:, 4] / 60 + den[:, 5] / 3600

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    e_density_cut = den[inicio:fin, 6]
    vel_cut = vel[inicio:fin, 6:-1]

    return t_cut, e_density_cut, vel_cut


t_swea_entero, density, vel_mso_xyz = importar_swea(ti, tf)
lpw, t_lpw, density_lpw = importar_lpw(ti, tf)

# plt.plot(t_swea_entero, density)
# plt.show()

sw_i = 17.5
sw_f = 17.8

ms_i = 18.06
ms_f = 18.21

"""
giroradio: rg = me * v_perp / (q_e * B)  en la región upstream
where vperp is the component of the velocity perpendicular to the
direction of the magnetic field and B is the strength of the magnetic field
"""

# primero diezmamos
paso = int(len(B_entero) / len(vel_mso_xyz))
B = B_entero[::paso]
t_mag = t_mag_entero[::paso]

for k in [(sw_i, sw_f), (ms_i, ms_f)]:

    # selecciono las regiones:
    mag_i = donde(t_mag, k[0])
    mag_f = donde(t_mag, k[1])
    swea_i = donde(t_swea_entero, k[0])
    swea_f = donde(t_swea_entero, k[1])
    lpw_i = donde(t_lpw, k[0])
    lpw_f = donde(t_lpw, k[1])

    # recorto
    B_cut = B[mag_i:mag_f] * 1e-9  # T
    velocidad = vel_mso_xyz[swea_i:swea_f]

    # ahora que están recortadas las convierto en lo que busco: la intensidad de B y la v perp
    B_medio = np.mean(B_cut, axis=0)
    B_medio_normalizado = B_medio / np.linalg.norm(B_medio)

    # proyecto la velocidad sobre B
    dot = np.dot(velocidad, B_medio_normalizado)
    N = np.zeros((len(dot), len(B_medio)))

    for i in range(len(N)):
        N[i, :] = dot[i] * B_medio_normalizado

    v_perp = np.mean(velocidad - N, axis=0)

    # el giroradio entonces:
    gf = me / (q_e * np.linalg.norm(B_medio))  # girofreq
    rg = gf * np.linalg.norm(v_perp)
    rg_min = (
        me
        * min(np.linalg.norm(velocidad - N, axis=1))
        / (q_e * np.linalg.norm(B_medio))
    )
    rg_max = (
        me
        * max(np.linalg.norm(velocidad - N, axis=1))
        / (q_e * np.linalg.norm(B_medio))
    )

    if k[0] == sw_i:
        print(
            f"El radio de Larmor de electrones SW es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
        )
    if k[0] == ms_i:
        print(
            f"El radio de Larmor de electrones MS es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
        )

    if k[0] == sw_i:
        print(f"La girofrecuencia de electrones SW es {gf:1.3g} s⁻¹")
    if k[0] == ms_i:
        print(f"La girofrecuencia de electrones MS es {gf:1.3g} s⁻¹")

    """
    Proyectando en la normal:
    """
    # normal = np.array([0.920,-0.302,0.251])
    v_normal = np.dot(np.mean(velocidad, axis=0), normal)

    # el giroradio entonces:
    rg_normal = me * np.linalg.norm(v_normal) / (q_e * np.linalg.norm(B_medio))

    if k[0] == sw_i:
        print(
            f"El radio de Larmor de electrones SW con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
        )
    if k[0] == ms_i:
        print(
            f"El radio de Larmor de electrones MS con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
        )

    """
    Longitud inercial
    """
    density_mean = np.zeros(swea_f - swea_i)  # upstream
    paso = 20  # cada paso son 4 segundos.
    if swea_i - paso > 0:  # si no se cumple, va a calcularlo mal
        for i in range(swea_f - swea_i):
            density_mean[i] = np.mean(
                density[swea_i + i - paso : swea_i + i]
            )  # toma desde atrás del ti así no se mete en la MPB nunca
        ion_length = 5.31e5 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
        ion_min = 5.31e5 / np.sqrt(max(density_mean)) * 1e-5  # km
        ion_max = 5.31e5 / np.sqrt(min(density_mean)) * 1e-5  # km

    if k[0] == sw_i:
        print(
            f"La longitud inercial de electrones usando SWEA SW es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g}"
        )
    if k[0] == ms_i:
        print(
            f"La longitud inercial de electrones usando SWEA MS es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g}"
        )

    density_mean = np.zeros(lpw_f - lpw_i)  # upstream
    paso = 20  # cada paso son 4 segundos.
    if lpw_i - paso > 0:  # si no se cumple, va a calcularlo mal
        for i in range(lpw_f - lpw_i):
            density_mean[i] = np.mean(
                density_lpw[lpw_i + i - paso : lpw_i + i]
            )  # toma desde atrás del ti así no se mete en la MPB nunca

        ion_length = 5.31e5 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
        ion_min = 5.31e5 / np.sqrt(max(density_mean)) * 1e-5  # km
        ion_max = 5.31e5 / np.sqrt(min(density_mean)) * 1e-5  # km

    if k[0] == sw_i:
        print(
            f"La longitud inercial de electrones usando LPW SW es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g}"
        )
    if k[0] == ms_i:
        print(
            f"La longitud inercial de electrones usando LPW MS es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g}"
        )


# ##########
# # guarda en la spreadsheet

# # hoja_MVA.update_acell(f"AA{fila}", f"{ion_length:1.3g}")
# # hoja_MVA.update_acell(f"AB{fila}", f"{rg:1.3g}")

# # hoja_Bootstrap.update_acell(f"Q{fila}", f"{ion_length:1.3g}")
# # hoja_Bootstrap.update_acell(f"R{fila}", f"{rg:1.3g}")

# # hoja_Ajuste.update_acell(f"Q{fila}", f"{ion_length:1.3g}")
# # hoja_Ajuste.update_acell(f"R{fila}", f"{rg:1.3g}")


# # en el SW
# inicio_sw = donde(t_swea_entero, 17.7)
# fin_sw = donde(t_swea_entero, 17.8)

# density_sw = np.zeros(fin_sw - inicio_sw)  # upstream
# paso = 20  # cada paso son 4 segundos.
# for i in range(fin_sw - inicio_sw):
#     density_sw[i] = np.mean(
#         density[inicio_sw + i : fin_sw + i]
#     )  # toma desde atrás del ti así no se mete en la MPB nunca

# ion_length = 2.28e07 / np.sqrt(np.mean(density_sw)) * 1e-5  # km
# ion_min = 2.28e07 / np.sqrt(max(density_sw)) * 1e-5  # km
# ion_max = 2.28e07 / np.sqrt(min(density_sw)) * 1e-5  # km
# print(
#     f"La longitud inercial de iones es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g})"
# )
