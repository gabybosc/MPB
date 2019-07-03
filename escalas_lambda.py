"""
Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
En el eje x están los tiempos en los cuales se centran las ventanas en las que barre. En el eje y están las escalas, que van a ir desde 1s hasta 40s de diámetro. El mapa de colores me va a dar el valor del cociente.
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import time as time
from funciones import find_nearest, Mij

date_entry = input('Enter a date in YYYY-DDD format \n')
# date_entry = '2016-076'
year, doy = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

ti = float(input('Tiempo inicial del barrido\n'))
tf = float(input('Tiempo final del barrido\n'))
# ti = 18.215
# tf = 18.245
n = int(ti*32*3600)
nn = int(tf*32 * 3600)
# path = '../../../MAVEN/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop
mag = np.loadtxt(path + f'MAG_hires/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=n, usecols=(0,1,2,3,4,5,6,7,8,9,10))
# mag = np.genfromtxt(path + f'MAG_hires/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skip_header=n, skip_footer=nn)


tiempo_central = np.zeros(int((tf-ti)*3600)) #la cantidad de segundos entre tf y ti
tiempo_central[0] = ti
for i in range(len(tiempo_central)-1):
    tiempo_central[i+1] = tiempo_central[i] + 1/3600 #el tiempo central se va barriendo cada 5 segundos

print(f'tiempo final efectivo = {tiempo_central[-1]}')

escalas = np.zeros(40)
escalas[0] = 1/3600 #la escala más chica es de 1s
for i in range(len(escalas)-1):
    escalas[i+1] = escalas[i] + 1/3600

# print(f'escala mayor = {escalas[-1]*3600}s')

dia = mag[:,1]
t = mag[:,6]  #el dia decimal
t = (t - dia) * 24 #para que me de sobre la cantidad de horas

M = np.size(t) #el numero de datos

#el campo
B = np.zeros((M, 3))
for i in range(7,10):
    B[:,i-7] = mag[:, i]


program_starts = time.time()
cociente = np.zeros((len(tiempo_central), len(escalas)))
for i in range(len(tiempo_central)):
    if i*100%len(tiempo_central) ==0:
        print(f'{i*100/len(tiempo_central)}%')
    for j in range(len(escalas)):
        inicio = np.where(t == find_nearest(t, tiempo_central[i] - escalas[j]))[0][0]
        fin = np.where(t == find_nearest(t, tiempo_central[i] + escalas[j]))[0][0]

        #ahora empieza el MVA con los datos que elegí
        B_cut = B[inicio:fin+1,:]

        M_ij = Mij(B_cut)

        #ahora quiero los autovectores y autovalores
        [lamb, x] = np.linalg.eigh(M_ij) #uso eigh porque es simetrica

        #Los ordeno de mayor a menor
        idx = lamb.argsort()[::-1]
        lamb = lamb[idx]
        x = x[:,idx]
        #ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
        x1 = x[:,0]
        x2 = x[:,1]
        x3 = x[:,2]
        if x3[0] < 0: #si la normal aputna para adentro me la da vuelta
            x3 = - x3

        cociente[i,j] = lamb[1]/lamb[2]

program_ends = time.time()

print(f'El loop tardó {program_ends-program_starts} s')

m = 1.67E-27 #kg
q = 1.6E-19 #C
periodo_ciclotron = 2*np.pi*m / (q * np.linalg.norm(B_cut, axis=1)) * 1E9 #en s
periodo_diezmado = np.zeros(len(tiempo_central))
k = len(periodo_ciclotron) / len(tiempo_central)
for i in range(len(periodo_diezmado)):
    periodo_diezmado[i] = periodo_ciclotron[int(i*k)]

matriz = np.zeros((len(tiempo_central)+1, len(escalas)+2))
matriz[0,2:] = escalas
matriz[1:,0] = periodo_diezmado
matriz[1:,1] = tiempo_central
matriz[1:,2:] = cociente

with open(f'outputs/cociente_lambdas_{doy}.txt','w') as file:
    file.write('La primera columna es el período de ciclotron, la primera fila son las escalas, la segunda columna es el tiempo central.\n')
    for i in range(len(matriz[:,0])):
        for j in range(len(matriz[0,:])):
            file.write(f'{matriz[i,j]}\t')
        file.write('\n')
