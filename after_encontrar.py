import numpy as np
import csv as csv
import pandas as pd
from funciones import donde
from collections import Counter

catalogo = np.genfromtxt(f"outputs/grupo2.txt", dtype="str")
fecha = np.genfromtxt(f"outputs/grupo2.txt", dtype="str", usecols=0)
hora_bs = np.genfromtxt(f"outputs/grupo2.txt", dtype="str", usecols=1)
hora_mpb = np.genfromtxt(f"outputs/grupo2.txt", dtype="str", usecols=2)
theta = np.genfromtxt(f"outputs/grupo2.txt", dtype="float", usecols=3)
pdyn = np.genfromtxt(f"outputs/grupo2.txt", dtype="float", usecols=4)
total = np.genfromtxt(f"outputs/hoja_grupo2.txt", dtype="str", usecols=0)
"""
Si ya tengo más de 100, vamos a ver si parece ser un buen spread en tiempos y ángulos o si necesito seguir
"""

m = []
for i in range(len(total)):
    year = total[i].split("-")[0]
    month = total[i].split("-")[1]
    m.append(year + month)
print(Counter(m))  # la cantidad que tengo de cada mes, para elegir

m = []
for i in range(len(fecha)):
    year = fecha[i].split("-")[0]
    month = fecha[i].split("-")[1]
    m.append(year + month)
print(Counter(m))  # la cantidad que elegí de cada mes

# veo que para el grupo 2 le faltan al menos 2020-06, 2020-07, 2020-11
# voy a tener que mirar en el catálogo grupo 2 y elegir esos índices para mirar con encontrar_grupo2
