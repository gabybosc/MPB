import numpy as np
import matplotlib.pyplot as plt
from funciones import fechas, tiempos, find_nearest_final, find_nearest_inicial

year, month, day, doy = 2017, 11, 24, 328#fechas()
ti_MVA, tf_MVA = 12.24, 12.28#tiempos()

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop

mag = np.loadtxt(path + 'MAG_pc.asc')
M = len(mag[:,0]) #el numero de datos
B = mag[:, 6:9]
Bnorm = np.linalg.norm(B, axis=1)
Bxyz_paraperp = mag[:,6:9]

hh = mag[:,3]
mm = mag[:,4]
ss = mag[:,5]

t = hh + mm/60 + ss/3600 #hdec

M = np.size(t) #el numero de datos

#la posición(x,y,z)
posicion = np.zeros((M, 3))
for i in range(9,12):
    posicion[:,i-9] = mag[:, i]

#la matriz diaria:
MD = np.zeros((M, 9))
MD[:, 0] = t
for i in range(1,4):
    MD[:, i] = B[:,i-1]
MD[:,4] = Bnorm
for i in range(5,8):
    MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390 #altitud en km

x = posicion[:,0]/3390
y = posicion[:,1]/3390
z = posicion[:,2]/3390

r = np.sqrt(x**2+y**2)

inicio = np.where(t == find_nearest_inicial(t, ti_MVA))[0][0]
fin = np.where(t == find_nearest_final(t, tf_MVA))[0][0]

MD_cut = MD[inicio:fin]
t_cut = t[inicio:fin]
posicion_cut = posicion[inicio:fin]
n_p = int(len(posicion_cut)/2)
