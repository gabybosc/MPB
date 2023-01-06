import numpy as np
import matplotlib.pyplot as plt
from funciones import fechas
from funciones_plot import plot_datetime
from importar_datos import importar_mag_1s

"""
Plotea sólo los datos de MAG de baja resolución
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
mag, t, B, posicion = importar_mag_1s(year, month, day, 0.1, 24)


# plt.figure()
# # plot_datetime(year, month, day, t, np.linalg.norm(B, axis=1))
# plt.xlabel("t (UTC)")
# plt.ylabel("|B|")
# plt.ylim(ymin=0, ymax=70)
# plt.title("MAG lowres hdec")


# plt.show(block=False)

with open("outputs/hora_grupo1.txt", "a") as file:
    for l in np.unique(horas):
        file.write(l)
        file.write("\n")
