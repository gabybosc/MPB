import numpy as np
import matplotlib.pyplot as plt
from funciones import find_nearest, unix_to_decimal, unix_to_timestamp
import datetime as dt
import matplotlib.dates as md
plt.ion()
# np.set_printoptions(precision=4)

year = '2016'
month = '12'
day = '01'

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
swea = np.loadtxt(path + 'diff_en_flux.asc')

energy = swea[:, 7]
JE_total = swea[:, -1]

hh = swea[:,3]
mm = swea[:,4]
ss = swea[:,5]

tt = hh + mm/60 + ss/3600 #hdec

"""
Ejemplo funcional para E = 100 eV
index_100eV = np.where(energy == find_nearest(energy, 100))[0] #en todos estos índices se midió 100eV, es uno por segundo

t = tt[index_100eV] #los tiempos en los cuales se midió E = 100eV

JE_100eV = JE[index_100eV] #los valores de JE para 100eV

plt.semilogy(t, JE_100eV)
"""

energias = [100 + i*10 for i in range(11)]

plt.figure()
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[0]
    t = tt[index]
    JE = JE_total[index]
    plt.semilogy(t, JE, label = energia)
plt.xlabel('tiempo')
plt.ylabel('JE')
plt.legend()
