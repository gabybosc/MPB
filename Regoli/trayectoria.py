import numpy as np
import matplotlib.pyplot as plt
from time import time
import sys

sys.path.append("..")

from funciones import donde

path = "../../../datos/simulacion_leonardo/"
# path2 = "../../../datos/MAG_hires/"
path2 = "../../../datos/MAG_1s/2016/"
variables = np.load(path + "rho_v_B_J.npy")
posicion = np.load(path + "pos_mhd.npy")

# nave = np.loadtxt(
#     path2 + "mvn_mag_l2_2016076ss_20160316_v01_r01.sts",
#     skiprows=160,
#     usecols=(2, 3, 4, 11, 12, 13),
# )

nave = np.loadtxt(
    path2 + "mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts",
    skiprows=160,
    usecols=(2, 3, 4, 11, 12, 13),
)


"""
Quiero los datos en la trayectoria de la nave.
"""

"""
Por como están hechos los datos, cuando una de las posiciones es cero, todas las
variables se anulan. Voy a borrar esos datos entonces.
"""


ceros = [i for i in range(len(posicion)) if posicion[i, 1] == 0]

pos_cut = np.delete(posicion, ceros, axis=0)
var_cut = np.delete(variables, ceros, axis=0)

val = np.concatenate((pos_cut, var_cut), axis=1)

t = nave[:, 0] + nave[:, 1] / 60 + nave[:, 2] / 3600  # hdec
inicial = donde(t, 17)
final = donde(t, 19)

pos_nave = nave[inicial:final, 3:] / 3390  # en RM
xn = pos_nave[:, 0]
yn = pos_nave[:, 1]
zn = pos_nave[:, 2]

xmin = min(xn)
xmax = max(xn)
ymax = max(yn)
ymin = min(yn)
zmin = min(zn)
zmax = max(zn)

program_starts = time()
xcut = np.array([val[i, :] for i in range(len(val)) if xmin < val[i, 0] <= xmax])
ycut = np.array([xcut[i, :] for i in range(len(xcut)) if ymin < xcut[i, 1] <= ymax])
zcut = np.array([ycut[i, :] for i in range(len(ycut)) if zmin < ycut[i, 2] <= zmax])
program_ends = time()

print(f"El recorte tardó {program_ends-program_starts:.2f} s")


X = np.zeros(len(pos_nave))
Y = np.zeros(len(pos_nave))
Z = np.zeros(len(pos_nave))
xc = zcut[:, 0]
yc = zcut[:, 1]

program_starts = time()
for i in range(len(xn)):
    X[i] = xc[donde(xc, xn[i])]
    Y[i] = yc[donde(yc, yn[i])]
    # Y[i] = donde(zcut[:,1], yn[i])
    # Z[i] = donde(zcut[:,2], zn[i])
program_ends = time()

X = np.zeros(len(pos_nave))
program_starts1 = time()
for i in range(len(xn)):
    X[i] = xc[donde(xc, xn[i])]
program_ends1 = time()

X = np.zeros(len(pos_nave))
program_starts2 = time()
for i in range(len(xn)):
    X[i] = xc[donde(zcut[:, 0], xn[i])]
program_ends2 = time()

print(f"loop1 {program_ends1-program_starts1:.2f} s")
print(f"loop2 {program_ends2-program_starts2:.2f} s")


plt.scatter(X, xn)
plt.scatter(Y, yn)
plt.show()
Y[20] - yn[20]
