import numpy as np

"""
Agarra todos los archivos que creé con IDL y los pasa a un único archivo.
Recorta si quiero.
"""

path = "../../../datos/simulacion_chuanfei/"
posicion = np.loadtxt(path + "salidas_idl/posicion.csv")  # r1, r2 (x, y o z)
densidades = np.loadtxt(path + "salidas_idl/rho.csv")  # rho, Hrho, Orho, O2rho, CO2rho
velocidad_h = np.loadtxt(path + "salidas_idl/velocidad_h.csv")  # Hvx Hvy Hvz
velocidad_O = np.loadtxt(path + "salidas_idl/velocidad_O.csv")  # Hvx Hvy Hvz
velocidad_O2 = np.loadtxt(path + "salidas_idl/velocidad_O2.csv")  # Hvx Hvy Hvz
velocidad_CO2 = np.loadtxt(path + "salidas_idl/velocidad_CO2.csv")  # Hvx Hvy Hvz
campo_B_b1 = np.loadtxt(path + "salidas_idl/campo.csv")  # Bx By Bz b1x b1y b1z
presion = np.loadtxt(path + "salidas_idl/presion.csv")  # Pe P HP OP O2P CO2P
corrientes = np.loadtxt(path + "salidas_idl/corrientes.csv")  # jx jy jz
grad_p = np.loadtxt(path + "salidas_idl/grad_p.csv")  # gradx grady gradz


"""
Unidades:
RM Mp/cc km/s nT nPa uA/m2 nPa/m
"""

arr = np.hstack(
    (
        posicion,
        densidades,
        velocidad_h,
        campo_B_b1,
        presion,
        corrientes,
        grad_p,
        velocidad_O,
        velocidad_O2,
        velocidad_CO2,
    )
)

# zona_interes = np.array(
#     [arr[i, :] for i in range(len(arr)) if -0.05 < arr[i, 1] <= 0.05]
# )

zona_interes = arr

reordenados = np.array(sorted(zona_interes, key=lambda f: f[0]))


np.savetxt(path + "filename.gz", reordenados)
