import numpy as np
import matplotlib.pyplot as plt
from funciones import SZA, donde
from importar_datos import importar_mag_1s


datos = np.loadtxt("outputs/para_sofi.txt")

dia = datos[:, 1]
hora = datos[:, 2]

i = donde(dia, 76)
j = donde(hora[i : i + 5], 18)

idx = i + j

mag, t, B, posicion = importar_mag_1s(
    2016, "03", int(dia[i] - 60), hora[idx] - 0.1, hora[idx] + 0.1
)  # tomo 6 min a cada lado de la MPB


"""
Estos parámetros son recontra aproximados, porque estoy considerando que el punto
que anoté en el txt es el punto medio de la MPB y que es igual de ancha para ambos
lados. Si realmente quisiera ver el B y eso debería hacer un mejor análisis.
"""
sza = SZA(posicion, int(len(posicion) / 2))
