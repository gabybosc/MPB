import numpy as np
import glob as glob
import sys
import matplotlib.pyplot as plt
from fit_venus import fit_Xu

sys.path.append("..")
from funciones import angulo, donde


"""
Agarra los cruces de la MPB de VEX y va a devolverme una lista 
con el SZA de cada uno. Compara con el fit de Xu 2021.
"""

# path = glob.glob("../../../VEX.txt")

# Ajuste de Xu de la MPB:
x_Xu, yz_Xu = fit_Xu()
mag = np.loadtxt("../../../VEX.asc")

posicion = mag[:, 6:]  # en km

# orbita = posicion / 6050  # es la posicion en RV
orbita = posicion
XX = orbita[:, 0]
YZ = np.sqrt(orbita[:, 1] ** 2 + orbita[:, 2] ** 2)

plt.plot(XX, YZ)
plt.plot(x_Xu, yz_Xu)
plt.show()

"""
Vamos a tirar todos los puntos donde la órbita esté lejos del planeta
y me quedo solo con el dayside:
Me quedo con yz < 2, 0 < x < 2
"""

idx = [i for i in range(len(XX)) if 0 < XX[i] < 2]
idx2 = [i for i in range(len(YZ)) if YZ[i] < 2]

indice = list(set(idx) & set(idx2))  # los índices que estén en ambas listas

x_cut = XX[indice]
yz_cut = YZ[indice]

# plt.plot(x_cut, yz_cut)
# plt.plot(x_Xu, yz_Xu)
# plt.show()


"""
lista con los puntos donde los x, yz se parecen
"""

a = [donde(x_cut, x_Xu[i]) for i in range(len(x_Xu))]
b = [donde(yz_cut, yz_Xu[i]) for i in range(len(yz_Xu))]


ab = list(set(a) & set(b))

plt.plot(x_cut, yz_cut)
plt.plot(x_Xu, yz_Xu)
plt.scatter(x_cut[a], yz_cut[a])
plt.show()

"""
pos y pos_xu tienen la misma longitud. Como las órbitas son bien portadas
(valores siempre positivos), puedo simplemente hacer la resta entre las
normas y ver cuándo es mínima.
Me quedo con ese punto entonces.
"""

pos = np.transpose([x_cut[a], yz_cut[a]])
pos_xu = np.transpose([x_Xu, yz_Xu])

idx = np.linalg.norm(pos - pos_xu, axis=1).argmin()

pos[idx]

plt.plot(x_cut, yz_cut)
plt.plot(x_Xu, yz_Xu)
plt.scatter(x_cut[idx], yz_cut[idx])
plt.show()


def SZA(posicion, index):
    SZA = angulo(posicion[index, :], [1, 0]) * 180 / np.pi
    return SZA


sza_mpb = SZA(pos, idx)

print(sza_mpb)
