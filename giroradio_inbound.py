from sys import exit
import numpy as np
from funciones import fechas, donde
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
year, month, day, doy = fechas()
ti = int(input("Hora del cruce (HH)\n"))
tf = ti + 1
in_out = input("Inbound? (y/n)\n")
if in_out == "n":
    print("este no sirve!")
    exit()
lst = input("La normal como nn,nn,nn\n").split(",")
normal = np.array([float(i) for i in lst])

datos = np.loadtxt("outputs/t1t2t3t4.txt")
for j in range(len(datos)):
    if (
        datos[j, 0] == float(year)
        and datos[j, 1] == float(doy)
        and int(datos[j, 2]) == int(ti)
    ):
        i = j
    else:
        print("No tengo el ancho de la MPB para esa fecha")
        exit()

t1 = datos[i, 2]
t2 = datos[i, 3]
t3 = datos[i, 4]
t4 = datos[i, 5]

mag, t_mag_entero, B_entero, posicion = importar_mag_1s(year, month, day, ti, tf)
swia, t_swia_entero, density, temperature, vel_mso_xyz = importar_swia(
    year, month, day, ti, tf
)

"""
giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
where vperp is the component of the velocity perpendicular to the
direction of the magnetic field and B is the strength of the magnetic field
"""
# selecciono la región upstream:
inicio_mag = donde(t_mag_entero, t1 - 0.1)
fin_mag = donde(t_mag_entero, t1 - 0.05)
inicio_swia = donde(t_swia_entero, t1 - 0.1)
fin_swia = donde(t_swia_entero, t1 - 0.05)

t_mag = t_mag_entero[inicio_mag:fin_mag]
B = B_entero[inicio_mag:fin_mag] * 1e-9  # T

t_swia = t_swia_entero[inicio_swia:fin_swia]
velocidad = vel_mso_xyz[inicio_swia:fin_swia]

# para poder proyectar v_para voy a tener que diezmar B
paso = int(len(B) / len(velocidad))
B_diezmado = B[::paso]

# ahora que están recortadas las convierto en lo que busco: la intensidad de B y la v perp
B_norm = np.linalg.norm(B, axis=1)

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
# normal = np.array([0.920,-0.302,0.251])
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
for i in range(fin_swia - inicio_swia):
    density_mean[i] = np.mean(
        density[inicio_swia + i - paso : inicio_swia + i]
    )  # toma desde atrás del ti así no se mete en la MPB nunca

ion_length = 2.28e07 / np.sqrt(np.mean(density_mean)) * 1e-5  # km
ion_min = 2.28e07 / np.sqrt(max(density_mean)) * 1e-5  # km
ion_max = 2.28e07 / np.sqrt(min(density_mean)) * 1e-5  # km
print(
    f"La longitud inercial de iones es {ion_length:1.3g} km (entre {ion_min:1.3g} y {ion_max:1.3g})"
)
