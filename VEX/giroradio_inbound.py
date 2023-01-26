import numpy as np
from importar_datos import importar_MAG_pds, importar_ELS_clweb, importar_fila
import sys

sys.path.append("..")

from funciones import fechas, donde

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
# year, month, day, doy = fechas()
# ti = int(input("Hora del cruce (HH)\n"))
# tf = ti + 1
year = 2011
month = "04"
day = 30
doy = 120
ti = 2
tf = 4

t, B_entero, posicion = importar_MAG_pds(year, doy, ti, tf)
IMA = np.loadtxt(f"../../../datos/clweb/VEX_IMA_vel.asc")
hh = IMA[:, 3]
mm = IMA[:, 4]
ss = IMA[:, 5]
t_IMA = hh + mm / 60 + ss / 3600
vel = IMA[:, 6:]
densidad = np.loadtxt(f"../../../datos/clweb/VEX_IMA_density.asc")[:, -1]

t1 = 2.778339444
"""
giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
where vperp is the component of the velocity perpendicular to the
direction of the magnetic field and B is the strength of the magnetic field
"""
# selecciono la región upstream:
inicio_mag = donde(t, t1 - 0.1)
fin_mag = donde(t, t1 - 0.05)
inicio_swia = donde(t_IMA, t1 - 0.1)
fin_swia = donde(t_IMA, t1 - 0.05)

t_mag = t[inicio_mag:fin_mag]
B = B_entero[inicio_mag:fin_mag] * 1e-9  # T

t_sw = t_IMA[inicio_swia:fin_swia]
velocidad = vel[inicio_swia:fin_swia]

# para poder proyectar v_para voy a tener que diezmar B
paso = int(len(B) / len(velocidad))
B_diezmado = B[::paso]

# ahora que están recortadas las convierto en lo que busco: la intensidad de B y la v perp
B_medio = np.mean(B_diezmado, axis=0)
B_medio_normalizado = B_medio / np.linalg.norm(B_medio)

dot = np.dot(velocidad, B_medio_normalizado)
N = np.zeros((len(dot), len(B_medio)))

for i in range(len(N)):
    N[i, :] = dot[i] * B_medio_normalizado

v_perp = np.mean(velocidad - N, axis=0)

# el giroradio entonces:
rg = mp * np.linalg.norm(v_perp) / (q_e * np.linalg.norm(B_medio))
rg_min = (
    mp * min(np.linalg.norm(velocidad - N, axis=1)) / (q_e * np.linalg.norm(B_medio))
)
rg_max = (
    mp * max(np.linalg.norm(velocidad - N, axis=1)) / (q_e * np.linalg.norm(B_medio))
)

print(f"El radio de Larmor es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})")

"""
Proyectando en la normal:
"""
normal = np.array([0.391, -0.129, 0.911])
# 0.391	-0.129	0.911
v_normal = np.dot(np.mean(velocidad, axis=0), normal)

# el giroradio entonces:
rg_normal = mp * np.linalg.norm(v_normal) / (q_e * np.linalg.norm(B_medio))

print(
    f"El radio de Larmor con la velocidad proyectada en la normal es {rg_normal:1.3g} km"
)


"""
Longitud inercial
"""
density_mean = np.zeros(fin_swia - inicio_swia)  # upstream
paso = 20  # cada paso son 4 segundos.
if inicio_swia - paso > 0:  # si no se cumple, va a calcularlo mal
    for i in range(fin_swia - inicio_swia):
        density_mean[i] = np.mean(
            densidad[inicio_swia + i - paso : inicio_swia + i]
        )  # toma desde atrás del ti así no se mete en la MPB nunca

    ion_length = 2.28e07 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
    ion_min = 2.28e07 / np.sqrt(max(density_mean)) * 1e-5  # km
    ion_max = 2.28e07 / np.sqrt(min(density_mean)) * 1e-5  # km
    print(
        f"La longitud inercial de iones es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g})"
    )
