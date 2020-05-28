import numpy as np
from importar_datos import importar_mag

import sys

sys.path.append("..")
from funciones import fechas, tiempos, donde, angulo, UTC_to_hdec

# year, month, day = 2016, "03", "16"  # fechas()
year, month, day, doy = fechas()
ti, tf = tiempos("Tiempo inicial y final de la magnetofunda")
mag, t, B, posicion = importar_mag(year, month, day, ti, tf)


"""ángulo entre el B de la magnetofunda y la corriente"""
inicio_mf = donde(t, ti)
fin_mf = donde(t, tf)
j = [float(x) for x in input("Enter the normal in N.NN N.NN N.NN format\n").split()] #[-2.8, -19.2, -12.6]
B_medio = np.mean(B, axis=0)

ang = angulo(j, B_medio) * 180/np.pi

print(f"El ángulo entre j y B magnetofunda es {ang}º")
