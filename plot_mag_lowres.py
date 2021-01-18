import numpy as np
import matplotlib.pyplot as plt
from funciones import fechas
from funciones_plot import plot_datetime

"""
Plotea sólo los datos de MAG de baja resolución
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()

n = 150

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = "../../datos/"  # path a los datos desde la laptop

mag = np.loadtxt(
    path + f"MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts",
    skiprows=n,
)

dia = mag[:, 1]
t = mag[:, 6]  # el dia decimal
t = (t - dia) * 24  # hdec

M = np.size(t)  # el numero de datos

# el campo
# el campo
B = np.zeros((M, 3))
for i in range(7, 10):
    B[:, i - 7] = mag[:, i]

# la posición(x,y,z)
posicion = np.zeros((M, 3))
for i in range(11, 14):
    posicion[:, i - 11] = mag[:, i]


plt.figure()
plot_datetime(year, month, day, t, np.linalg.norm(B, axis=1))
plt.xlabel("t (UTC)")
plt.ylabel("|B|")
plt.title("MAG lowres hdec")


plt.show(block=False)
