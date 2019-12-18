import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cdflib as cdf
import datetime as dt
from funciones import find_nearest, deltaB, unix_to_decimal, unix_to_timestamp, Bpara_Bperp, fechas, tiempos,datenum
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

year, month, day, doy = 2016,'03', 16, 76#fechas()
ti, tf = 17.5, 18.5#tiempos()

path = f'../../../datos/clweb/{year}-{month}-{day}/'
datos = np.loadtxt('../outputs/t1t2t3t4.txt')
# for j in range(len(datos)):
#     if datos[j,0] == float(year) and datos[j,1] == float(doy) and int(datos[j,2]) == int(ti):
#         i = j

# t1 = datos[i,2]
# t2 = datos[i,3]
# t3 = datos[i,4]
# t4 = datos[i,5]
t1, t2,t3,t4 = 18.2167,18.2204,18.235,18.2476
tiempos = np.array([t1,t2,t3,t4])



mag, t, B, posicion = importar_mag(year, month, day, ti, tf)
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

swea, t_swea, energias = importar_swea(year, month, day, ti, tf)

energy = swea[:,7]
JE_total = swea[:,-1]
inicio_swea = np.where(t_swea == find_nearest(t_swea, ti))[0][0]
fin_swea = np.where(t_swea == find_nearest(t_swea, tf))[0][0]

###############################################################################################SWIA

swia, t_swia, density = importar_swia(year, month, day, ti, tf)

inicio_swia = np.where(t_swia == find_nearest(t_swia, ti))[0][0]
fin_swia = np.where(t_swia == find_nearest(t_swia, tf))[0][0]
############################################################################################### LPW
lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)

inicio_lpw = np.where(t_lpw == find_nearest(t_lpw, ti))[0][0]
fin_lpw = np.where(t_lpw == find_nearest(t_lpw, tf))[0][0]



index = np.array((int(year), int(day)))


swea, t_swea, energias = importar_swea(year, month, day, ti, tf)
energy = swea[:,7]
JE_total = swea[:,-1]

inicio_swea = np.where(t_swea == find_nearest(t_swea, ti))[0][0]
fin_swea = np.where(t_swea == find_nearest(t_swea, tf))[0][0]
###############################################################################################SWIA

swia, t_swia, density = importar_swia(year, month, day, ti,tf)

############################################################################################### LPW
lpw, t_lpw, e_density = importar_lpw(year, month, day,ti,tf)


year = int(year)
month = int(month)
day = int(day)
tiempo_mag = np.array([np.datetime64(datenum(year, month, day, x)) for x in t]) #datenum es una función mía
tiempo_swea = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swea])
tiempo_swia = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swia])
tiempo_lpw = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_lpw])
tiempo_low = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_plot]) #datenum es una función mía



tm1 = np.where(t == find_nearest(t, t1))
tm2 = np.where(t == find_nearest(t, t2))
tm3 = np.where(t == find_nearest(t, t3))
tm4 = np.where(t == find_nearest(t, t4))
tiempo_lim = [tiempo_mag[tm1],tiempo_mag[tm2],tiempo_mag[tm3],tiempo_mag[tm4]]

xfmt = md.DateFormatter('%H:%M:%S')

plt.clf()#clear figure
fig = plt.figure(1, constrained_layout=True)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(top = 0.93, bottom = 0.07, left = 0.05,right=0.95, hspace = 0.005, wspace=0.15)
fig.set_size_inches(15, 10)#con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((3,2),(0,0))
plt.plot(tiempo_low, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / B')
plt.plot(tiempo_low, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / B')
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel(r'|$\Delta B$|/ B')

ax2 = plt.subplot2grid((3,2),(1,0), sharex=ax1)
ax2.plot(tiempo_mag, B[:,0], label='Bx')
ax2.plot(tiempo_mag, B[:,1], label='By')
ax2.plot(tiempo_mag, B[:,2], label='Bz')
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_ylabel('Bx, By, Bz (nT)')

ax3 = plt.subplot2grid((3,2),(2,0), sharex=ax1)
plt.plot(tiempo_mag, Bnorm)
ax3.set_ylabel('|B| (nT)')
ax3.set_xlabel('Tiempo (UTC)')

ax4 = plt.subplot2grid((3,2),(0,1), sharex=ax1)
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[0]
    JE = JE_total[index]
    plt.semilogy(tiempo_swea[inicio_swea:fin_swea], JE[inicio_swea:fin_swea], label = f'{energia} eV')
ax4.set_ylabel('diff en flux')


ax5 = plt.subplot2grid((3,2),(1,1), sharex=ax1)
ax5.set_ylabel('Densidad de p+ \n del SW (cm⁻³)')
plt.plot(tiempo_swia, density)

ax6 = plt.subplot2grid((3,2),(2,1), sharex=ax1)
ax6.set_ylabel('Densidad total \n de e- (cm⁻³)')
ax6.set_xlabel('Tiempo (hdec)')
plt.semilogy(tiempo_lpw, e_density)


for ax in [ax1, ax2, ax4,ax5]:
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.legend()

for ax in [ax1, ax2, ax3,ax4,ax5,ax6]:
    ax.set_xlim(tiempo_mag[0],tiempo_mag[-1])
    ax.xaxis.set_major_formatter(xfmt)
    ax.grid()
    for xc in tiempo_lim:
        plt.axvline(x = xc, color = 'k', linewidth=1)

plt.suptitle(f'MAVEN {year}-{month}-{day}')
multi = MultiCursor(fig.canvas, (ax1, ax2, ax3,ax4,ax5,ax6), color='r', lw=1)




plt.show(block=False)
# plt.savefig(f'../outputs/figs_MPB/MPB_y{year}_d{doy}_t{int(t1)}.png', dpi=200)
# pickle.dump(fig, open(f'../outputs/figs_MPB/MPB_y{year}_d{doy}_t{int(t1)}.pkl', 'wb'))
