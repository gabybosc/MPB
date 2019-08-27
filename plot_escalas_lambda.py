"""
Plotea los archivos que devuelve escalas_lambda.py
"""

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from funciones import find_nearest, array_datenums
from funciones_plot import imshow_UTC, plot_datetime
import matplotlib.dates as md
plt.ion()


date_entry = input('Enter a date in YYYY-DDD format \n')
hora = input('Hora de la orbita \n')
# date_entry = '2016-066'
year, doy = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

tiempos_txt = np.loadtxt('outputs/t1t2t3t4.txt')
for i in range(len(tiempos_txt)):
    if int(year) == int(tiempos_txt[i,0]) and int(doy) == int(tiempos_txt[i,1]) and int(hora) == int(tiempos_txt[i,2]):
        tiempos = [tiempos_txt[i,2], tiempos_txt[i,3],tiempos_txt[i,4], tiempos_txt[i,5]]
timestamps = array_datenums(year, month, day, tiempos)#lo convierto a datenum

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop
datos = np.loadtxt(f'outputs/cociente_lambdas_d{doy}_t{hora}.txt', skiprows = 1)

periodo_ciclotron = datos[1:,0]
tiempo_central = datos[1:,1]
escalas = datos[0,2:] * 3600
cociente = np.transpose(datos[1:,2:])

ti = tiempo_central[0]-0.5
tf = tiempo_central[-1]+0.5
n = int(ti*32*3600)

mag = np.loadtxt(path + f'MAG_hires/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=n, usecols=(1,6,7,8,9))

dia = mag[:,0]
t = mag[:,1]  #el dia decimal
t = (t - dia) * 24 #hdec

M = np.size(t) #el numero de datos

#el campo
B = np.zeros((M, 3))
for i in range(2,5):
    B[:,i-2] = mag[:, i]

inicio = np.where(t == find_nearest(t, ti))[0][0]
fin = np.where(t == find_nearest(t, tf))[0][0]

B_norm = np.linalg.norm(B, axis = 1)
B_cut = B_norm[inicio:fin]
t_cut = t[inicio:fin]

inicio_MVA = np.where(t == find_nearest(t, tiempo_central[0]))[0][0]
fin_MVA = np.where(t == find_nearest(t, tiempo_central[-1]))[0][0]
B_MVA = B_norm[inicio_MVA:fin_MVA]
t_MVA = t[inicio_MVA:fin_MVA]


plt.figure()
imshow_UTC(year, month, day, tiempo_central, cociente, escalas, 'inferno')
# plot_datetime(year, month, day,t_MVA, B_MVA, 'cyan', '-', 1, 0.5) #no sé por qué se superpone mal, tiene mal los tiempos.
for tt in timestamps:
    plt.axvline(x = tt, color = 'g')
plt.title('Heatmap del cociente de lambdas en distintas escalas temporales \n y el campo magnético superpuesto')
plt.xlabel('Tiempo en el que está centrado (hh:mm:ss)')
plt.ylabel('Radio (s) \n |B| (nT)')

plt.figure()
plot_datetime(year, month, day,t_cut, B_cut, 'red', '-', 1, 1)
plt.ylabel('|B| (nT)')
plt.xlabel('Tiempo UTC (hh:mm:ss)')

"""
pensar criterio para la escala maxima
"""
