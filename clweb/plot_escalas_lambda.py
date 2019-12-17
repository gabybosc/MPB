"""
Plotea los archivos que devuelve escalas_lambda.py
"""

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from funciones import find_nearest, array_datenums, fechas
from funciones_plot import imshow_UTC, plot_datetime
import matplotlib.dates as md
import os
plt.ion()

year, month, day, doy = fechas()
hora = input('Hora (HH) \n')

tiempos_txt = np.loadtxt('../outputs/t1t2t3t4.txt')
for i in range(len(tiempos_txt)):
    if int(year) == int(tiempos_txt[i,0]) and int(doy) == int(tiempos_txt[i,1]) and int(hora) == int(tiempos_txt[i,2]):
        tiempos = [tiempos_txt[i,2], tiempos_txt[i,3],tiempos_txt[i,4], tiempos_txt[i,5]]
timestamps = array_datenums(year, month, day, tiempos)#lo convierto a datenum

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
if os.path.isfile(path + 'mag_filtrado.txt'):
    mag = np.loadtxt(path + 'mag_filtrado.txt', skiprows=2)
    M = len(mag[:,0]) #el numero de datos
    B = mag[:, :3]

    Bnorm = mag[:,-1]
    mag = np.loadtxt(path + 'MAG.asc')
    Bxyz_paraperp = mag[:,6:9]
else:
    mag = np.loadtxt(path + 'MAG.asc')
    M = len(mag[:,0]) #el numero de datos
    B = mag[:, 6:9]
    Bnorm = np.linalg.norm(B, axis=1)


datos = np.loadtxt(f'../outputs/cociente_lambdas_d{doy}_t{hora}.txt', skiprows = 1)

periodo_ciclotron = datos[1:,0]
tiempo_central = datos[1:,1]
escalas = datos[0,2:] * 3600
cociente = np.transpose(datos[1:,2:])

ti = tiempo_central[0]-0.5
tf = tiempo_central[-1]+0.5
n = int(ti*32*3600)


hh = mag[:,3]
mm = mag[:,4]
ss = mag[:,5]

t = hh + mm/60 + ss/3600 #hdec

M = np.size(t) #el numero de datos

inicio = np.where(t == find_nearest(t, ti))[0][0]
fin = np.where(t == find_nearest(t, tf))[0][0]

B_cut = Bnorm[inicio:fin]
t_cut = t[inicio:fin]

inicio_MVA = np.where(t == find_nearest(t, tiempo_central[0]))[0][0]
fin_MVA = np.where(t == find_nearest(t, tiempo_central[-1]))[0][0]
B_MVA = Bnorm[inicio_MVA:fin_MVA]
t_MVA = t[inicio_MVA:fin_MVA]



plt.figure()
imshow_UTC(year, month, day, tiempo_central, cociente, escalas, 'inferno', 3)
# plot_datetime(year, month, day,t_MVA, B_MVA, 'cyan', '-', 1, 0.5) #no sé por qué se superpone mal, tiene mal los tiempos.
for tt in timestamps:
    plt.axvline(x = tt, color = 'g') #plotea los tiempos t1t2t3t4
plt.title(f'Heatmap del cociente de lambdas en distintas escalas temporales \n  para el día {day}-{month}-{year}')
plt.xlabel('Tiempo en el que está centrado (hh:mm:ss)')
plt.ylabel('Radio (s) \n |B| (nT)')

plt.figure()
plot_datetime(year, month, day,t_cut, B_cut, 'red', '-', 1, 1)
plt.ylabel('|B| (nT)')
plt.xlabel('Tiempo UTC (hh:mm:ss)')

"""
pensar criterio para la escala maxima
"""
