import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")

from funciones import donde

path = "../../../datos/simulacion_chuanfei/"
datos_hr = np.loadtxt(path + "ejex_new3_+.gz")  # high resolution


pi = 1000
pf = 1250

mu0 = 4e-7 * np.pi  # T m / A
mp = 1.67e-27  # proton mass, kg
e_SI = 1.6e-19  # C
g = 3.7  # Mars surface gravity, m/s²

"""
r1 r2 rho Hrho Orho O2rho CO2rho H_vx H_vy H_vz
Bx By Bz b1x b1y b1z Pe P HP OP
O2P CO2P jx jy jz gradx grady gradz O_vx O_vy
O_vz O2_vx O2_vy O2_vz CO2_vx CO2_vy CO2_vz e_vx e_vy e_vz
"""

x = datos_hr[:, 0]  # RM
z = datos_hr[:, 1]  # RM
B = datos_hr[:, 10:13]  # nT
b1 = datos_hr[:, 13:16]  # nT
J = datos_hr[:, 22:25]  # ua/m2
grad_p = datos_hr[:, 25:28]  # nPa/m
presion = {
    "e-": datos_hr[:, 16],
    "H+": datos_hr[:, 18],
    "O+": datos_hr[:, 19],
    "O2+": datos_hr[:, 20],
    "CO2+": datos_hr[:, 21],
}  # nPa
densities = {
    "e-": datos_hr[:, 2],
    "H+": datos_hr[:, 3],
    "O+": datos_hr[:, 4],
    "O2+": datos_hr[:, 5],
    "CO2+": datos_hr[:, 6],
}  # Mp/cc
velocities = {
    "H+": datos_hr[:, 7:10],
    "O+": datos_hr[:, 28:31],
    "O2+": datos_hr[:, 31:34],
    "CO2+": datos_hr[:, 34:37],
    "e-": datos_hr[:, 37:],
}  # km/s

xi = donde(x, 1.225)
xf = donde(x, 1.230)

B_upstream = np.mean(B[xi:xf, :], axis=0)  # nT

np_up = np.mean(densities["H+"][xi:xf])  # nT
velp_up = np.mean(velocities["H+"][xi:xf], axis=0)  # nT
vele_up = np.mean(velocities["e-"][xi:xf], axis=0)  # nT

print(f"B = {np.linalg.norm(B_upstream)}, n = {np_up}, v_p = {np.linalg.norm(velp_up)}, v_e = {np.linalg.norm(vele_up)}")
