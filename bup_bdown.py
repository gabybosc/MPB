from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from scipy.stats import norm
from datetime import datetime
from funciones import hodograma, error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, set_axes_equal,datenum

np.set_printoptions(precision=4)

path = '../datos/marzo 2016/16/' #path a los datos
datos = np.loadtxt(path + 'mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148) #lee todo y me da todo
n =2
datos = datos[:-n, :] #borra las ultimas 2 filas, que es ya el dia siguiente (no sé si siempre)
lpw = np.loadtxt(path + 'mvn_kp_insitu_20160316_v14_r03_orbita18h.csv') #son los datos entre las 18 y las 19h
t_lpw = lpw[:,0] + lpw[:,1]/60 + lpw[:,2]/3600

ti = 18.227
tf = 18.235
dia = datos[:,1]
t = datos[:,6]  #el dia decimal
t = (t - dia) * 24 #para que me de sobre la cantidad de horas

M = np.size(t) #el numero de datos

#tengo que asegurarme de que no haya agujeros en mis datos
for i in range(M-1):
    if t[i+1] - t[i] > 24 * 1.5e-5: #1.1e-5 es 1s, le doy un poco más
        print('salto en la linea {} de {} segundos'.format(i+144, (t[i+1] - t[i]) / (24*1.1e-5)))

#el campo
B = np.zeros((M, 3))
for i in range(7,10):
    B[:,i-7] = datos[:, i]

#la posición(x,y,z)
posicion = np.zeros((M, 3))
for i in range(11,14):
    posicion[:,i-11] = datos[:, i]

#la matriz diaria:
MD = np.zeros((M, 9))
MD[:, 0] = t
for i in range(1,4):
    MD[:, i] = B[:,i-1]
MD[:,4] = np.linalg.norm(B, axis=1)#la norma de B
for i in range(5,8):
    MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390 #altitud en km

#si quiero elegir entre ciertas horas:
t1 = find_nearest_inicial(t, ti)
t2 = find_nearest_final(t, tf)
inicio = np.where(t == t1)[0][0]
fin = np.where(t == t2)[0][0]

t_1 = np.where(t == find_nearest(t,18.2167))[0][0]
t_2 = np.where(t == find_nearest(t,18.2204))[0][0]
t_3 = np.where(t == find_nearest(t,18.235))[0][0]
t_4 = np.where(t == find_nearest(t,18.2476))[0][0]

tiempo_mag = np.array([np.datetime64(datenum(2016, 3, 16, x)) for x in t]) #datenum es una función mía

#########
Bu = np.zeros((180, 4))
for i in range(180):
    paso = 18.2 + 0.0001 * i
    inicio_up = np.where(t == find_nearest_inicial(t, paso))[0][0]
    fin_up = t_1
    B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0) #nT
    # B_upstream_max = max(np.linalg.norm(B[inicio_up:fin_up, :], axis=1)) #nT
    Bu[i, 0:3] = B_upstream
    Bu[i, 3] = np.linalg.norm(B_upstream)
plt.figure(1)
plt.plot(tiempo_mag[t_1-180:t_1], Bu[:,0:3])
plt.figure(2)
plt.plot(tiempo_mag[t_1-180:t_1],Bu[:,3])

Bd = np.zeros((200, 4))
for i in range(200):
    paso = 18.248 + 0.0001 * i
    inicio_down = t_4
    fin_down = np.where(t == find_nearest_final(t, paso))[0][0]
    B_downstream = np.mean(B[inicio_down:fin_down, :], axis=0) #nT
    Bd[i, 0:3] = B_downstream
    Bd[i, 3] = np.linalg.norm(B_downstream)
plt.figure(3)
plt.plot(tiempo_mag[t_4:t_4+200],Bd[:,0:3])
plt.figure(4)
plt.plot(tiempo_mag[t_4:t_4+200],Bd[:,3])


plt.show(block=False)
