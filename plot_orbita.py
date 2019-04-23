from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from scipy.stats import norm
from datetime import datetime
from funciones import hodograma, error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, set_axes_equal
# from funcion_flujo_energia_cdf import flujo_energia

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

#################
#ahora empieza el MVA con los datos que elegí
MD_cut = MD[inicio : fin+1, :]
M_cut = np.size(MD_cut[:,0])
B_cut = B[inicio:fin+1, :]

M_ij = Mij(B_cut)
# np.savetxt('outs/matriz_%d'%dia[0], M_ij) #guarda la matriz de cada dia

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

#las proyecciones
B1 = np.dot(B_cut, x1)
B2 = np.dot(B_cut, x2)
B3 = np.dot(B_cut, x3)

#el B medio
B_medio_vectorial = np.array([np.mean(B_cut[0]), np.mean(B_cut[1]), np.mean(B_cut[2])])


###############
#elijo la órbita a partir de 17.85
x_nave = posicion[:,0] / 3390
y_nave = posicion[:,1] / 3390
z_nave = posicion[:,2] / 3390

t_orb = np.where(t == find_nearest(t, 17.85))[0][0]
x_final = np.where(x_nave == find_nearest(x_nave[fin:-1], x_nave[t_orb]))[0][0] #me da la posición más cercana a la inicial

#posicion cada diez minutos
arr = np.zeros(34)
arr[0] = 17.85
for i in range(33):
    arr[i+1] = arr[i] + 0.01667

idx = [np.where(t == find_nearest(t, xx))[0][0] for xx in arr]

marte = plt.Circle((0, 0), 1, color='r', alpha = 0.2)
marte1 = plt.Circle((0, 0), 1, color='r', alpha = 0.2)
marte2 = plt.Circle((0, 0), 1, color='r', alpha = 0.2)


fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
ax.add_artist(marte)

ax.plot(x_nave[t_orb:x_final], y_nave[t_orb:x_final])
ax.scatter(x_nave[idx], y_nave[idx])
ax.set_aspect('equal', adjustable='datalim')
ax.set_xlabel(r'$X_{MSO}$')
ax.set_ylabel(r'$Y_{MSO}$')

fig1, ax1 = plt.subplots() # note we must use plt.subplots, not plt.subplot
ax1.add_artist(marte1)

ax1.plot(x_nave[t_orb:x_final], z_nave[t_orb:x_final])
ax1.scatter(x_nave[idx], z_nave[idx])
ax1.set_aspect('equal', adjustable='datalim')
ax1.set_xlabel(r'$X_{MSO}$')
ax1.set_ylabel(r'$Z_{MSO}$')

fig2, ax2 = plt.subplots() # note we must use plt.subplots, not plt.subplot
ax2.add_artist(marte2)

ax2.plot(y_nave[t_orb:x_final], z_nave[t_orb:x_final])
ax2.scatter(y_nave[idx], z_nave[idx])
ax2.set_aspect('equal', adjustable='datalim')
ax2.set_xlabel(r'$Y_{MSO}$')
ax2.set_ylabel(r'$Z_{MSO}$')

# plt.plot(x_nave[t_orb:x_final], np.sqrt(y_nave[t_orb:x_final]**2 + z_nave[t_orb:x_final]**2))

plt.show()
