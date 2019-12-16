import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cdflib as cdf
import datetime as dt
from funciones import find_nearest, deltaB, unix_to_decimal, unix_to_timestamp, Bpara_Bperp, fechas, tiempos
from importar_datos import importar_mag, importar_lpw, importar_swea, importar_swia
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import MultiCursor
import time
import matplotlib.dates as md
import matplotlib.cm as cm
import os as os
import pickle


"""
Con los tiempos t1 t2 t3 t4 ploteo los datos de MAG, SWEA, SWIA y LPW. Les dibujo encima líneas que correspondan a los tiempos.
Guarda las figuras como png y como pickle (para poder reabrirlas con python y que sean interactivas)
"""
#se fija que los datos no estén vacíos antes de cargarlos (para ahorrarme el error)


np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = tiempos()

path = f'../../../datos/clweb/{year}-{month}-{day}/'
datos = np.loadtxt('../outputs/t1t2t3t4.txt')
for j in range(len(datos)):
    if datos[j,0] == float(year) and datos[j,1] == float(doy) and int(datos[j,2]) == int(ti):
        i = j

t1 = datos[i,2]
t2 = datos[i,3]
t3 = datos[i,4]
t4 = datos[i,5]
tiempos = np.array([t1,t2,t3,t4])

mag, t, B, posicion = importar_mag(year, month, day)
Bnorm = np.linalg.norm(B, axis=1)

mag_low = np.loadtxt(path + 'mag_1s.sts', skiprows=160)
tlow = mag_low[:,6]  #el dia decimal
tlow = (tlow -int( doy)) * 24 #para que me de sobre la cantidad de horas

Mlow = np.size(tlow) #el numero de datos
#el campo
Blow = np.zeros((Mlow, 3))
for i in range(7,10):
    Blow[:,i-7] = mag_low[:, i]


ti = t1 - 0.15
tf = t4 + 0.15

inicio = np.where(t == find_nearest(t, ti))[0][0]
fin = np.where(t == find_nearest(t, tf))[0][0]
B_para, B_perp_norm, t_plot = Bpara_Bperp(Blow, tlow, ti, tf)
###############################################################################################SWEA

swea, t_swea, energias = importar_swea(year, month, day)

energy = swea[:,7]
JE_total = swea[:,-1]
inicio_swea = np.where(t_swea == find_nearest(t_swea, ti))[0][0]
fin_swea = np.where(t_swea == find_nearest(t_swea, tf))[0][0]

###############################################################################################SWIA

swia, t_swia, density = importar_swia(year, month, day)

inicio_swia = np.where(t_swia == find_nearest(t_swia, ti))[0][0]
fin_swia = np.where(t_swia == find_nearest(t_swia, tf))[0][0]
############################################################################################### LPW
lpw, t_lpw, e_density = importar_lpw(year, month, day)

inicio_lpw = np.where(t_lpw == find_nearest(t_lpw, ti))[0][0]
fin_lpw = np.where(t_lpw == find_nearest(t_lpw, tf))[0][0]



index = np.array((int(year), int(day)))


plt.clf()#clear figure
fig = plt.figure(1, constrained_layout=True)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(top = 0.93, bottom = 0.07, left = 0.05,right=0.95, hspace = 0.005, wspace=0.15)
fig.set_size_inches(15, 10)#con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((3,2),(0,0))
plt.plot(t_plot, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / B')
plt.plot(t_plot, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / B')
plt.setp(ax1.get_xticklabels(), visible=False)
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
ax1.set_ylabel(r'|$\Delta B$|/ B')
ax1.grid()
ax1.legend()

ax4 = plt.subplot2grid((3,2),(1,0), sharex=ax1)
ax4.plot(t[inicio:fin], B[inicio:fin,0], label='Bx')
ax4.plot(t[inicio:fin], B[inicio:fin,1], label='By')
ax4.plot(t[inicio:fin], B[inicio:fin,2], label='Bz')
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.set_ylabel('Bx, By, Bz (nT)')
ax4.legend()
ax4.grid()

ax3 = plt.subplot2grid((3,2),(2,0), sharex=ax1)
plt.plot(t[inicio:fin], Bnorm[inicio:fin])
ax3.grid()
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
ax3.set_ylabel('|B| (nT)')
ax3.set_xlabel('Tiempo (hdec)')

ax5 = plt.subplot2grid((3,2),(0,1), sharex=ax1)
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[0]
    JE = JE_total[index]
    plt.semilogy(t_swea[inicio_swea:fin_swea], JE[inicio_swea:fin_swea], label = f'{energia} eV')
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
ax5.grid()
ax5.set_xlabel('tiempo (hdec)')
ax5.set_ylabel('diff en flux')
ax5.legend()



ax7 = plt.subplot2grid((3,2),(1,1), sharex=ax1)
plt.setp(ax7.get_xticklabels(), visible=False)
ax7.set_ylabel('Densidad de p+ \n del SW (cm⁻³)')
plt.plot(t_swia[inicio_swia:fin_swia], density[inicio_swia:fin_swia])
ax7.grid()
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)


ax6 = plt.subplot2grid((3,2),(2,1), sharex=ax1)
ax6.set_ylabel('Densidad total \n de e- (cm⁻³)')
ax6.set_xlabel('Tiempo (hdec)')
plt.semilogy(t_lpw[inicio_lpw:fin_lpw], e_density[inicio_lpw:fin_lpw])
ax6.grid()
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)

plt.suptitle(f'MAVEN {year}-{month}-{day}')
multi = MultiCursor(fig.canvas, (ax1, ax7, ax3,ax4,ax5,ax6), color='r', lw=1)




plt.show(block=False)
plt.savefig(f'../outputs/figs_MPB/MPB_y{year}_d{doy}_t{int(t1)}.png', dpi=200)
pickle.dump(fig, open(f'../outputs/figs_MPB/MPB_y{year}_d{doy}_t{int(t1)}.pkl', 'wb'))
