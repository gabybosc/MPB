import numpy as np
from funciones import donde, long_inercial_iones, giroradio, giroradio_termico
from importar_datos import importar_mag_1s, importar_swia

np.set_printoptions(precision=4)

"""
Calcula el giroradio y la long inercial de iones. El giroradio lo calcula usando
v_perp al campo B y también usando la velocidad proyectada en la normal.
"""

# #########CONSTANTES
mp = 1.67e-27  # masa del proton en kg
kB = 1.38e-23  # cte de Boltzmann en J/K
q_e = 1.602e-19  # carga del electrón en C

# ##########DATOS
year, month, day, doy = 2016, "03", 16, "076"
ti = 17
tf = 19

normal = [0.920, -0.302, 0.251]
t1, t2, t3, t4 = 18.2167, 18.2204, 18.235, 18.2476

mag, t_mag_entero, B_entero, posicion = importar_mag_1s(year, month, day, ti, tf)
swia, t_swia_entero, density, temperature_all, vel_mso_xyz = importar_swia(
    year, month, day, ti, tf
)

# def importar_temp(ti, tf):
#     path = "../../datos/temp_swica_ms.asc"
#     swica = np.loadtxt(path)
#
#     t = swica[:, 3] + swica[:, 4] / 60 + swica[:, 5] / 3600
#
#     inicio = donde(t, ti)
#     fin = donde(t, tf)
#
#     t_cut = t[inicio:fin]
#     temp = swica[inicio:fin, -1]
#
#     return t_cut, temp


# plt.plot(t_swia_entero, density)
# plt.show()

sw_i = 17.5
sw_f = 17.8

ms_i = 18.06
ms_f = 18.21

"""
giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
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
    swia_i = donde(t_swia_entero, k[0])
    swia_f = donde(t_swia_entero, k[1])
    # tt, temperature = importar_temp(ms_i, ms_f)

    # recorto
    B_cut = B[mag_i:mag_f]
    velocidad = vel_mso_xyz[swia_i:swia_f]
    temperature = temperature_all[swia_i:swia_f]

    # ahora que están recortadas las convierto en lo que busco: la intensidad de B y la v perp
    B_medio = np.mean(B_cut, axis=0)
    B_medio_normalizado = B_medio / np.linalg.norm(B_medio)

    # proyecto la velocidad sobre B
    dot = np.dot(velocidad, B_medio_normalizado)
    N = np.zeros((len(dot), len(B_medio)))

    for i in range(len(N)):
        N[i, :] = dot[i] * B_medio_normalizado

    v_perp = np.mean(velocidad - N, axis=0)

    gf, rg, rg_min, rg_max = giroradio(B_cut, velocidad)
    if k[0] == sw_i:
        print(f"La girofrecuencia de iones SW es {gf:1.3g} s⁻¹")
    if k[0] == ms_i:
        print(f"La girofrecuencia de iones MS es {gf:1.3g} s⁻¹")

    if k[0] == sw_i:
        print(
            f"El radio de Larmor de iones SW es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
        )
    if k[0] == ms_i:
        print(
            f"El radio de Larmor de iones MS es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
        )
    gth = giroradio_termico(B_cut, temperature)
    if k[0] == sw_i:
        print(f"El giroradio térmico de iones SW es {gth:1.3g} km")
    if k[0] == ms_i:
        print(f"El giroradio térmico de iones MS es {gth:1.3g} km")
    """
    Proyectando en la normal:
    """
    # normal = np.array([0.920,-0.302,0.251])
    v_normal = np.dot(np.mean(velocidad, axis=0), normal)

    # el giroradio entonces:
    rg_normal = mp * np.linalg.norm(v_normal) / (q_e * np.linalg.norm(B_medio * 1e-9))

    if k[0] == sw_i:
        print(
            f"El radio de Larmor de iones SW con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
        )
    if k[0] == ms_i:
        print(
            f"El radio de Larmor de iones MS con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
        )

    """
    Longitud inercial
    """
    ion_length, ion_min, ion_max = long_inercial_iones(density[swia_i:swia_f])
    if k[0] == sw_i:
        print(
            f"La longitud inercial de iones SW es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g}"
        )
    if k[0] == ms_i:
        print(
            f"La longitud inercial de iones MS es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g}"
        )
