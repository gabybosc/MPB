import numpy as np

"""
Agarra todos los archivos que creé con IDL y los pasa a un único archivo.
Recorta si quiero.
Después uso recortar.py
"""

path = "../../../datos/simulacion_chuanfei/nueva_simu/"
posicion = np.loadtxt(path + "y=0_position.csv")  # r1, r2 (x, y o z)
densidades = np.loadtxt(path + "y=0_rho.csv")  # rho, Hrho, Orho, O2rho, CO2rho
velocity_h = np.loadtxt(path + "y=0_velocity_h.csv")  # Hvx Hvy Hvz
velocity_e = np.loadtxt(path + "y=0_velocity_e.csv")  # Hvx Hvy Hvz
velocity_O = np.loadtxt(path + "y=0_velocity_O.csv")  # Hvx Hvy Hvz
velocity_O2 = np.loadtxt(path + "y=0_velocity_O2.csv")  # Hvx Hvy Hvz
velocity_CO2 = np.loadtxt(path + "y=0_velocity_CO2.csv")  # Hvx Hvy Hvz
campo_B_b1 = np.loadtxt(path + "y=0_magnetic.csv")  # Bx By Bz b1x b1y b1z
presion = np.loadtxt(path + "y=0_pressure.csv")  # Pe P HP OP O2P CO2P
corrientes = np.loadtxt(path + "y=0_current.csv")  # jx jy jz
grad_p = np.loadtxt(path + "y=0_grad_p.csv")  # gradx grady gradz


"""
Unidades:
RM Mp/cc km/s nT nPa uA/m2 nPa/m
"""

arr = np.hstack(
    (
        posicion,
        densidades,
        velocity_h,
        campo_B_b1,
        presion,
        corrientes,
        grad_p,
        velocity_O,
        velocity_O2,
        velocity_CO2,
        velocity_e,
    )
)

# zona_interes = np.array(
#     [arr[i, :] for i in range(len(arr)) if -0.05 < arr[i, 1] <= 0.05]
# )
zona_interes = arr


reordenados = np.array(sorted(zona_interes, key=lambda f: f[0]))


np.savetxt(path + "filename.gz", reordenados)
