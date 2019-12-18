import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from funciones import find_nearest, fechas,tiempos
from funciones_plot import plot_datetime
from importar_datos import importar_mag

"""
Plotea sólo los datos de MAG de alta resolución
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = tiempos()
# year, month, day, doy = 2016, '03', 16, 76
# ti, tf = 17.85, 18.4

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

Bnorm = np.linalg.norm(B, axis=1)
Bxyz_paraperp = mag[:,6:9]

# inicio = np.where(t == find_nearest(t, ti))[0][0]
# fin = np.where(t == find_nearest(t, tf))[0][0]
#
x = posicion[:, 0]
y = posicion[:, 1]
z = posicion[:, 2]

r = np.sqrt(x**2 + y**2 + z**2)
phi = np.arctan2(y,x)
theta = np.arccos(z/r)

polares = np.transpose(np.array([r,phi,theta]))
pnorm = polares/np.linalg.norm(polares, axis=0)

Bpolar = np.zeros((len(Bnorm),3))

for i in range(len(Bnorm)):
    Bpolar[i,:] = Bnorm[i] * pnorm[i,:]


# B_cut = B_norm[inicio:fin]
# t_cut = t[inicio:fin]
fig = plt.figure()
fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.12,right=0.95, hspace = 0.0, wspace=0.15)
plt.xticks( rotation=25 )

ax1 = plt.gca()
ax1 = plt.subplot2grid((2,1),(0,0))
ax1.plot(t, Bnorm)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('|B| (nT)')

ax2 = plt.subplot2grid((2,1),(1,0), sharex=ax1)
ax2.plot(t, Bpolar)
ax2.set_ylabel('B polar')
ax2.set_xlabel('Time (hdec)')
plt.show(block= False)
