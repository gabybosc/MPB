import numpy as np

path = "../../../datos/simulacion_leonardo/"
# variables = np.loadtxt(path + "variables.txt")

posicion = np.load(path + "pos_cut.npy")  # cut es porque no tienen los ceros
variables = np.load(path + "var_cut.npy")

"""
Rho Ux Uy Uz Bx By Bz P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x b1y b1z e jx jy jz
"""

"""
Recorte de numpy
"""


val = np.concatenate((posicion, variables), axis=1)
limite = 1

# a = [np.abs(posicion[i, 2]) <= 0.05 for i in range(len(posicion))]
# idx = [i for i, x in enumerate(a) if x]
# zcero = posicion[idx, :]

zcero = np.array([val[i, :] for i in range(len(val)) if np.abs(val[i, 2]) <= 1.5])
ycero = np.array([zcero[i, :] for i in range(len(zcero)) if np.abs(zcero[i, 1]) <= 1.5])
zona_interes = np.array(
    [ycero[i, :] for i in range(len(ycero)) if 1 < ycero[i, 0] <= 2]
)

# quiero reordenar los x de forma creciente
reordenados = np.array(sorted(zona_interes, key=lambda f: f[0]))

np.save(path + "cubo_Daniel.npy", reordenados)
