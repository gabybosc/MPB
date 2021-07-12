import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag
import sys

sys.path.append("..")

from funciones import fechas

"""
Plotea sólo B para todo el día
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = 0.01, 24

path = f"../../../datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la laptop

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

Bnorm = np.linalg.norm(B, axis=1)

plt.plot(t, Bnorm)
plt.ylabel("|B| (nT)")
plt.xlabel("Time (hdec)")
plt.show(block=False)
