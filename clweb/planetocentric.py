import numpy as np
import sys

sys.path.append("..")
from funciones import find_nearest_final, fechas, tiempos

"""
Toma los datos en coordenadas pc y devuelve la latitud y longitud del cruce.
Considera el cruce como el punto medio entre t2 y t3.
"""

year, month, day, doy = 2016, "03", 16, 76
t2, t3 = 18.2204, 18.235

# year, month, day, doy = fechas()
# t2, t3 = tiempos()

path = f"../../../datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la laptop

mag = np.loadtxt(path + "MAG_pc.asc")
M = len(mag[:, 0])  # el numero de datos
B = mag[:, 6:9]
Bnorm = np.linalg.norm(B, axis=1)
Bxyz_paraperp = mag[:, 6:9]

hh = mag[:, 3]
mm = mag[:, 4]
ss = mag[:, 5]

t = hh + mm / 60 + ss / 3600  # hdec

M = np.size(t)  # el numero de datos

# la posición(x,y,z)
posicion = np.zeros((M, 3))
for i in range(9, 12):
    posicion[:, i - 9] = mag[:, i]

# la matriz diaria:
MD = np.zeros((M, 9))
MD[:, 0] = t
for i in range(1, 4):
    MD[:, i] = B[:, i - 1]
MD[:, 4] = Bnorm
for i in range(5, 8):
    MD[:, i] = posicion[:, i - 5] / 3390  # en radios marcianos
MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390  # altitud en km

cruce = np.where(t == find_nearest_final(t, (t2 + t3) / 2))[0][0]
x = posicion[cruce, 0]
y = posicion[cruce, 1]
z = posicion[cruce, 2]


r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
phi = np.arctan2(y, x)
theta = np.arccos(z / r)

longitude = phi / np.pi * 180
latitude = theta / np.pi * 180 - 90

if latitude < 0:
    print(f"latitude = {-latitude:.3g}ºN, longitude = {longitude:.3g}ºE")
else:
    print(f"latitude = {latitude:.3g}ºS, longitude = {longitude:.3g}ºE")
