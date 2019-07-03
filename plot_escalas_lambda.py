"""
Plotea los archivos que devuelve escalas_lambda.py
"""

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from funciones import find_nearest

date_entry = input('Enter a date in YYYY-DDD format \n')\
# date_entry = '2016-066'
year, doy = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop
datos = np.loadtxt(f'outputs/cociente_lambdas_{doy}.txt', skiprows = 1)

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

plt.figure()
plt.imshow(cociente,aspect = 'auto',origin = 'lower', extent=(tiempo_central[0], tiempo_central[-1], escalas[0], escalas[-1]), cmap='inferno', vmax=30)
plt.plot(tiempo_central, periodo_ciclotron / 3600)
plt.colorbar()
plt.xlabel('Tiempo en el que está centrado (hdec)')
plt.ylabel('Diámetro (hdec)')
plt.title('Heatmap del cociente de lambdas en distintas escalas temporales')

plt.figure()
plt.plot(t_cut, B_cut)
plt.xlabel('t (hdec)')
plt.ylabel('|B|')
plt.grid()
plt.title('MAG hires')
plt.show(block= False)


"""
pasar todoa  segundos
ver si es radio o diametro
poner el modulo de B arriba o superpuesto
pensar criterio para la escala maxima
"""
