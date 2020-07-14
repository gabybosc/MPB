import numpy as np
from importar_datos import importar_mag
import scipy.signal as signal
import sys
import matplotlib.pyplot as plt

sys.path.append("..")
from funciones import fechas, tiempos, donde, angulo, UTC_to_hdec, Mij

# year, month, day = 2016, "03", "16"  # fechas()
# year, month, day, doy = fechas()
# ti, tf = tiempos("Tiempo inicial y final de la magnetofunda")
# mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
#
#
# """ángulo entre el B de la magnetofunda y la corriente"""
# inicio_mf = donde(t, ti)
# fin_mf = donde(t, tf)
# j = [float(x) for x in input("Enter the normal in N.NN N.NN N.NN format\n").split()] #[-2.8, -19.2, -12.6]
# B_medio = np.mean(B, axis=0)
#
# ang = angulo(j, B_medio) * 180/np.pi
#
# print(f"El ángulo entre j y B magnetofunda es {ang}º")

year, month, day, doy = 2016, "03", 16, 76
ti_MVA, tf_MVA = UTC_to_hdec("18:13:33"), UTC_to_hdec("18:14:06")
ti, tf = UTC_to_hdec("17:55:00"), UTC_to_hdec("18:30:00")

date_entry = f"{year}-{month}-{day}"

mag, t, B, posicion = importar_mag(year, month, day, ti_MVA, tf_MVA)
# ya los importa cortados a los datos, entonces no hace falta que haga el cut yo

M = len(t)
Bnorm = np.linalg.norm(B, axis=1)

# la matriz diaria:
MD = np.zeros((M, 9))
MD[:, 0] = t
for i in range(1, 4):
    MD[:, i] = B[:, i - 1]
MD[:, 4] = Bnorm
for i in range(5, 8):
    MD[:, i] = posicion[:, i - 5] / 3390  # en radios marcianos
MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390  # altitud en km

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
if any(np.cross(x1, x2) - x3) > 0.01:
    print("Cambio el signo de x1 para que los av formen terna derecha")
    x1 = -x1

print("la normal del MVA es ", x3)

# las proyecciones
B1 = np.dot(B, x1)
B2 = np.dot(B, x2)
B3 = np.dot(B, x3)

# el B medio
B_medio_vectorial = np.mean(B, axis=0)

print("cociente de lambdas = ", lamb[1] / lamb[2])
B_norm_medio = np.linalg.norm(B_medio_vectorial)

# filtra
Tseg = 4
fp = 1 / Tseg / 16  # fp < fs para que sea pasabajos
fs = 1 / 16  # f normalizada, da lo mismo si es omega o frec
N, Wn = signal.buttord(fp, fs, 3, 50)
b, a = signal.butter(N, Wn, "low")
Bx_filtrado = signal.filtfilt(b, a, B[:, 0])
By_filtrado = signal.filtfilt(b, a, B[:, 1])
Bz_filtrado = signal.filtfilt(b, a, B[:, 2])

B_filtrado = np.transpose([Bx_filtrado, By_filtrado, Bz_filtrado])
M_ijfilt = Mij(B_filtrado)

# ahora quiero los autovectores y autovalores
[lambfilt, xfilt] = np.linalg.eigh(M_ijfilt)  # uso eigh porque es simetrica

# Los ordeno de mayor a menor
idx = lambfilt.argsort()[::-1]
lambfilt = lamb[idx]
xfilt = xfilt[:, idx]
# ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
x1_filt = xfilt[:, 0]
x2_filt = xfilt[:, 1]
x3_filt = xfilt[:, 2]

if x3_filt[0] < 0:  # si la normal aputna para adentro me la da vuelta
    x3_filt = -x3_filt
if any(np.cross(x1_filt, x2_filt) - x3_filt) > 0.01:
    print("Cambio el signo de x1 para que los av formen terna derecha")
    x1_filt = -x1_filt

print("la normal del MVA es ", x3_filt)

# las proyecciones
B1_filt = np.dot(B_filtrado, x1_filt)
B2_filt = np.dot(B_filtrado, x2_filt)
B3_filt = np.dot(B_filtrado, x3_filt)

# el B medio
B_medio_vectorial_filt = np.mean(B_filtrado, axis=0)

print("cociente de lambdas = ", lambfilt[1] / lambfilt[2])
B_norm_medio_filt = np.linalg.norm(B_medio_vectorial_filt)


f, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), sharex=True)

ax1.plot(B2, B1, zorder=1)
ax1.plot(B2_filt, B1_filt, zorder=1)
ax2.plot(B3, B1, zorder=1)
ax2.plot(B3_filt, B1_filt, zorder=1)
ax1.scatter(B2[0], B1[0], s=50, zorder=2, marker="o", color="r", label="start")
ax1.scatter(B2[-1], B1[-1], s=50, zorder=2, marker="x", color="r", label="end")
ax2.scatter(B3[0], B1[0], s=50, zorder=2, marker="o", color="r", label="start")
ax2.scatter(B3[-1], B1[-1], s=50, zorder=2, marker="x", color="r", label="end")
ax1.set_xlabel(f"B2 (nT)", fontsize=16)
ax2.set_xlabel(f"B3 (nT)", fontsize=16)
ax1.set_ylabel(f"B1 (nT)", fontsize=16)
ax2.set_ylabel(f"B1 (nT)", fontsize=16)
ax1.grid()
ax2.grid()
ax1.tick_params(axis="both", which="major", labelsize=14)
ax2.tick_params(axis="both", which="major", labelsize=14)
plt.suptitle("MAVEN MAG MVA", fontsize=18)
plt.legend(fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
