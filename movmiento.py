import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag_1s, importar_swia
from funciones import fechas, tiempos, Mij

year, month, day, doy = fechas()
ti, tf = tiempos()

mag, t_mag, B, posicion = importar_mag_1s(year, month, day, ti, tf)
swia, t_swia, density, temperature, vel_mso = importar_swia(year, month, day, ti, tf)

v_maven = np.zeros((len(posicion), 3))
for i in range(len(posicion)-1):
    v_maven[i,:] = (posicion[i+1,:] - posicion[i,:]) / (t_mag[i+1] - t_mag[i]) /3600

#calculo la normal
M = len(t_mag)
Bnorm = np.linalg.norm(B, axis=1)

n_p = int(M / 2)

M_ij = Mij(B)

# ahora quiero los autovectores y autovalores
[lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

# Los ordeno de mayor a menor
idx = lamb.argsort()[::-1]
lamb = lamb[idx]
x = x[:, idx]
# ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
x1 = x[:, 0]
x2 = x[:, 1]
x3 = x[:, 2]

if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
    x3 = -x3

v_maven = np.zeros((len(posicion), 3))
vpara = np.zeros((len(posicion), 3))
vperp = np.zeros((len(posicion), 3))
for i in range(len(posicion)-1):
    v_maven[i,:] = (posicion[i+1,:] - posicion[i,:]) / (t_mag[i+1] - t_mag[i]) /3600
    vpara[i,:] = np.dot(v_maven[i,:], x3) * x3
    vperp[i,:] = v_maven[i,:] - vpara[i,:]
