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
limite = 0.05

zcero = np.array([val[i, :] for i in range(len(val)) if np.abs(val[i, 2]) <= limite])
ycero = np.array(
    [zcero[i, :] for i in range(len(zcero)) if np.abs(zcero[i, 1]) <= limite]
)
zona_interes = np.array(
    [ycero[i, :] for i in range(len(ycero)) if 0 < ycero[i, 0] <= 5]
)

# quiero reordenar los x de forma creciente
reordenados = np.array(sorted(zona_interes, key=lambda f: f[0]))

# np.save(path+'recorte.npy', reordenados)
