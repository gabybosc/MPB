import numpy as np
import matplotlib.pyplot as plt
from funciones import donde, fechas, tiempos
from funciones_plot import plot_datetime
from importar_datos import importar_mag

"""
Plotea sólo los datos de MAG de alta resolución 
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = tiempos()

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

inicio = donde(t, ti)
fin = donde(t, tf)

if any(posicion[inicio:fin, 2]) > 0:
    print("Norte")

B_norm = np.linalg.norm(B, axis=1)
B_cut = B_norm[inicio:fin]
t_cut = t[inicio:fin]

plt.figure()
plt.plot(t_cut, B_cut)
plt.xlabel("t (hdec)")
plt.ylabel("|B|")
plt.title("MAG hires hdec")

plt.figure()
plot_datetime(year, month, day, t_cut, B_cut)
plt.xlabel("t (UTC)")
plt.ylabel("|B|")
plt.title("MAG hires UTC")

plt.show(block=False)
