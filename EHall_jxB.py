import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
import matplotlib.dates as md
import datetime as dt
from funciones import find_nearest, unix_to_decimal, plot_select, set_axes_equal,find_nearest_final, find_nearest_inicial, Mij, datenum

np.set_printoptions(precision=4)
'''
Este código calcula EHall a partir del rotor del campo B: E ~ rotB x B, pero sin calcular numéricamente el rotor, sino que usando la ecuación con B_para y B_perp
'''

###########DATOS
path = '../datos/marzo 2016/16/'
cdf_swia = cdf.CDF(path + 'mvn_swi_l2_onboardsvymom_20160316_v01_r01.cdf')
lpw = np.loadtxt(path + 'mvn_kp_insitu_20160316_v14_r03_orbita18h.csv') #son los datos entre las 18 y las 19h
datos = np.loadtxt(path + 'mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148)
n =2
datos = datos[:-n, :]

t_lpw = lpw[:,0] + lpw[:,1]/60 + lpw[:,2]/3600

t_unix = cdf_swia.varget('time_unix')
density = cdf_swia.varget('density') #cgs
v_mso_imported = cdf_swia.varget('velocity_mso') #SI

t_swia = unix_to_decimal(t_unix)
inicio_swia = np.where(t_swia == find_nearest(t_swia, 17.85))[0][0]
fin_swia = np.where(t_swia == find_nearest(t_swia, 18.4))[0][0]

t_swia_cut = t_swia[inicio_swia:fin_swia]
density_cut = density[inicio_swia:fin_swia]

#datos de B
dia = datos[:,1]
t_mag = datos[:,6]
t_mag = (t_mag - dia) * 24
B = np.empty((len(datos), 3))
for i in range(7,10):
    B[:,i-7] = datos[:, i] #nT

posicion = np.empty((len(datos), 3))
for i in range(11,14):
    posicion[:,i-11] = datos[:, i] #km

inicio = np.where(t_mag == find_nearest_inicial(t_mag, 17.85))[0][0]
fin = np.where(t_mag == find_nearest_final(t_mag, 18.4))[0][0]

B_cut = B[inicio:fin+1, :]

M_ij = Mij(B_cut)
[lamb, x] = np.linalg.eigh(M_ij) #uso eigh porque es simetrica
#Los ordeno de mayor a menor
idx = lamb.argsort()[::-1]
lamb = lamb[idx]
x = x[:,idx]
#ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
x1 = x[:,0]
x2 = x[:,1]
x3 = x[:,2]

inicio_mpb = np.where(t_mag == find_nearest(t_mag, 18.227))[0][0]
fin_mpb = np.where(t_mag == find_nearest(t_mag, 18.235))[0][0]
ancho_mpb = 63 #km

inicio_up = np.where(t_mag == find_nearest_inicial(t_mag, 18.20))[0][0]
fin_up = np.where(t_mag == find_nearest_final(t_mag, 18.21))[0][0]
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0) #nT

inicio_down = np.where(t_mag == find_nearest_inicial(t_mag, 18.24))[0][0]
fin_down = np.where(t_mag == find_nearest_final(t_mag, 18.26))[0][0]
B_downstream = np.mean(B[inicio_down:fin_down,:], axis=0) #nT

mu = 4* np.pi * 1E-7 #Henry/m

J_s = np.cross(x3, (B_upstream-B_downstream)) / mu #nA/m

J_v = J_s / (1000*ancho_mpb) #nA/m²

print('La corriente superficial es Js = {0:1.3g} nA/m, la corriente en volumen es Jv = {1:1.3g} nA/m²'.format(np.linalg.norm(J_s), np.linalg.norm(J_v)))


#en realidad, quiero la densidad de electrones antes de la mpb, elijo el intervalo en el cual voy a graficar
e_density = lpw[:,3]
ti_lpw = np.where(t_lpw == find_nearest(t_lpw, 18.20))[0][0]
tf_lpw = np.where(t_lpw == find_nearest(t_lpw, 18.227))[0][0]
n_e = np.nanmean(e_density[ti_lpw:tf_lpw]) #hace el mean ignorando los nans #cm⁻³
n_e = n_e * 1E6 #m⁻³
# n_e = 1E7
q_e = 1.6E-19 #carga electron #C

E_Hall = np.cross(J_v * 1E-9, B[inicio_up:fin_down, :] * 1E-9) / (q_e * n_e) #V/m
E_Hall_max = np.cross(J_v * 1E-9, B[inicio_down+38, :] * 1E-9) / (q_e * n_e) #V/m
print('El campo de Hall max es {} V/m'.format(E_Hall_max))

tiempo_mag = np.array([np.datetime64(datenum(2016, 3, 16, x)) for x in t_mag[inicio_up:fin_down]]) #datenum es una función mía
t1 = np.where(t_mag[inicio_up:fin_down] == find_nearest(t_mag, 18.2193))[0][0]
t2 = np.where(t_mag[inicio_up:fin_down] == find_nearest(t_mag, 18.227))[0][0]
t3 = np.where(t_mag[inicio_up:fin_down] == find_nearest(t_mag, 18.235))[0][0]
t4 = np.where(t_mag[inicio_up:fin_down] == find_nearest(t_mag, 18.2476))[0][0]


plt.figure()
plt.subplots_adjust(bottom=0.2)
plt.xticks( rotation=25 )
ax = plt.gca()
xfmt = md.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(xfmt)
ax.plot(tiempo_mag, np.linalg.norm(E_Hall, axis=1))
for xc in [tiempo_mag[t1],tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
plt.ylabel('E Hall (V/m)')
plt.xlabel('Tiempo')
plt.grid()


plt.show(block = False)
