import numpy as np
from socket import gethostname


if gethostname() == "magneto2":
    path = f"../../../../../media/gabybosc/datos/Chuanfei/"
elif gethostname() == "gabybosc":
    path = "../../../datos/simulacion_chuanfei/"
# variables = np.loadtxt(path + "variables.txt")


"""
Rho Ux Uy Uz Bx By Bz P HpRho HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP
OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP b1x b1y b1z e jx jy jz
"""

"""
Antes de usar este tengo que usar convertir_en_un_txt.py
Recorte de numpy
"""

# posicion = np.load(path + "pos_cut.npy")  # cut es porque no tienen los ceros
# variables = np.load(path + "var_cut.npy")

# val = np.concatenate((posicion, variables), axis=1)
# limite = 1

# a = [np.abs(posicion[i, 2]) <= 0.05 for i in range(len(posicion))]
# idx = [i for i, x in enumerate(a) if x]
# zcero = posicion[idx, :]


# zcero = np.array([val[i, :] for i in range(len(val)) if np.abs(val[i, 2]) <= 2])
# ycero = np.array([zcero[i, :] for i in range(len(zcero)) if np.abs(zcero[i, 1]) <= 2])
# zona_interes = np.array(
#     [ycero[i, :] for i in range(len(ycero)) if 0 < ycero[i, 0] <= 2]
# )
#
# # quiero reordenar los x de forma creciente
# reordenados = np.array(sorted(zona_interes, key=lambda f: f[0]))
#
# np.savetxt(path + "cubo_2RM.gz", reordenados)


# si ya estoy en el eje x y sólo quiero seguir cortando
val = np.loadtxt(path + "nueva_simu/ejex.gz")

zcero = np.array([val[i, :] for i in range(len(val)) if 0 <= val[i, 1] <= 0.05])
reordenados = np.array(sorted(zcero, key=lambda f: f[0]))
np.savetxt(path + "nueva_simu/ejex+.gz", reordenados)
