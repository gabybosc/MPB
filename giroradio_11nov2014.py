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
# year, month, day = 2014, 11, 11
# ti = 4
# tf = 4.2
year, month, day = 2014, 11, 12
ti = 12.1
tf = 12.3

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

"""
giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
where vperp is the component of the velocity perpendicular to the
direction of the magnetic field and B is the strength of the magnetic field
"""

# primero diezmamos
paso = int(len(B_entero) / len(vel_mso_xyz))
B = B_entero[::paso]
t_mag = t_mag_entero[::paso]

# ahora que están recortadas las convierto en lo que busco: la intensidad de B y la v perp
B_medio = np.mean(B, axis=0)
B_medio_normalizado = B_medio / np.linalg.norm(B_medio)

# proyecto la velocidad sobre B
dot = np.dot(vel_mso_xyz, B_medio_normalizado)
N = np.zeros((len(dot), len(B_medio)))

for i in range(len(N)):
    N[i, :] = dot[i] * B_medio_normalizado

v_perp = np.mean(vel_mso_xyz - N, axis=0)

gf, rg, rg_min, rg_max = giroradio(B, vel_mso_xyz)
print(f"La frecuencia angular de iones SW es {gf:1.3g} rad s⁻¹")
print(
    f"El radio de Larmor de iones SW es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
)
gth = giroradio_termico(B, temperature_all)
print(f"El giroradio térmico de iones SW es {gth:1.3g} km")

"""
Longitud inercial
"""
ion_length, ion_min, ion_max = long_inercial_iones(density)
print(
    f"La longitud inercial de iones SW es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g}"
)
