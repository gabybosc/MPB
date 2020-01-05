import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from funciones import find_nearest, fechas,tiempos
from funciones_plot import plot_datetime
from importar_datos import importar_mag

"""
Plotea los datos de B en coordenadas esféricas. Sirve para descartar campos corticales.
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = tiempos()
# year, month, day, doy = 2016, '03', 16, 76
# ti, tf = 17.85, 18.4

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

Bnorm = np.linalg.norm(B, axis=1)

Br = np.zeros(len(Bnorm))
posicion_normalizada = np.zeros((len(Bnorm), 3))

for i in range(len(Bnorm)):
    posicion_normalizada[i,:] = posicion[i,:] / np.linalg.norm(posicion[i,:])
    Br[i] = np.dot(B[i,:], posicion_normalizada[i,:])

Btg = Bnorm - Br

fig = plt.figure()
fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.12,right=0.95, hspace = 0.0, wspace=0.15)
plt.xticks( rotation=25 )

ax1 = plt.gca()
ax1 = plt.subplot2grid((2,1),(0,0))
ax1.plot(t, Bnorm)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('|B| (nT)')

ax2 = plt.subplot2grid((2,1),(1,0), sharex=ax1, sharey=ax1)
ax2.plot(t, Br, label = 'r')
ax2.plot(t, Btg, label = 'tg')
ax2.set_ylabel('B polares')
ax2.set_xlabel('Time (hdec)')
ax2.legend()
plt.show(block= False)
