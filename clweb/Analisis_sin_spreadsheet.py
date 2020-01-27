import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
import scipy.signal as signal
import datetime as dt
import pdb
import os
from matplotlib.mlab import normpdf
from scipy.stats import norm
from funciones import error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, next_available_row, datenum, unix_to_decimal, UTC_to_hdec, fechas, tiempos
from funciones_MVA import ajuste_conico, plot_velocidades, plot_FLorentz, plot_bootstrap, bootstrap
from funciones_plot import hodograma, set_axes_equal
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from MVA_sin_spreadsheet import MVA

"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
"""


# year, month, day, doy = fechas()
# ti_MVA, tf_MVA = tiempos()
year, month, day, doy = 2016, '03', 16, 76
ti_MVA, tf_MVA = UTC_to_hdec('18:13:33'), UTC_to_hdec('18:14:06')

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
mag = np.loadtxt(path + 'MAG.asc')
swia = np.loadtxt(path + 'SWIA.asc')
lpw = np.loadtxt(path + 'LPW.asc')

x3, normal_boot, normal_fit, t, B, posicion, inicio, fin, B_cut, t1, t2, t3, t4,B_medio_vectorial = MVA(year, month, day, doy, ti_MVA, tf_MVA, mag)

#########
#buscamos el ángulo entre las normales
angulo_mva = np.arccos(np.clip(np.dot(normal_fit, x3), -1.0, 1.0)) #el clip hace que si por algun motivo el dot me da >1 (i.e. 1,00002), me lo convierte en 1
angulo_boot = np.arccos(np.clip(np.dot(normal_boot, x3), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!


##############
#Calculo la velocidad de la nave
v_punto = np.zeros((fin-inicio, 3))
norma_v = np.zeros(fin-inicio)
posicion_cut = posicion[inicio : fin+1, :]
t_cut = t[inicio : fin+1] * 3600 #en segundos
for i in range(fin-inicio):
    v_punto[i,:] = (posicion[inicio+1,:] - posicion[inicio]) / (t_cut[i+1]-t_cut[i]) #en km/s
    norma_v[i] = np.linalg.norm(v_punto[i,:])
#veamos que no cambia mucho punto a punto, usemos la norma
diff = max(norma_v)- min(norma_v)
#la velocidad promedio
v_media = np.array([np.mean(v_punto[:,0]), np.mean(v_punto[:,1]), np.mean(v_punto[:,2])])

#ahora quiero ver si la nave atraviesa perpendicularmente a la MPB
v_media_norm = v_media/np.linalg.norm(v_media) #la normalizo por cuestion de  cuentas nomas
angulo_v_fit = np.arccos(np.clip(np.dot(normal_fit, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_v_mva = np.arccos(np.clip(np.dot(x3, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_v_boot = np.arccos(np.clip(np.dot(normal_boot, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!

B_intermedio = B_medio_vectorial/np.linalg.norm(B_medio_vectorial) #la normalizo por cuestion de  cuentas nomas
angulo_B_fit = np.arccos(np.clip(np.dot(normal_fit, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_B_mva = np.arccos(np.clip(np.dot(x3, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_B_boot = np.arccos(np.clip(np.dot(normal_boot, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!

######
##Espesor de la MPB
#ahora veamos v_para
deltat_14 = (t4 - t1) * 3600
deltat_23 = (t3 - t2) * 3600

v_para = np.dot(v_media, normal_fit) * normal_fit
x_14_fit = v_para * deltat_14 #en km# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x_23_fit = v_para * deltat_23

#si ahora proyecto sobre la normal de la MVA
v_para_MVA = np.dot(v_media, x3) * x3
x_14_MVA = v_para_MVA * deltat_14 #en km
x_23_MVA = v_para_MVA * deltat_23

#si ahora proyecto sobre la normal del bootstrap
v_para_boot = np.dot(np.array([-2.487,  0.479,  2.836]), normal_boot) * normal_boot
x_14_boot = v_para_boot * deltat_14 #en km
x_23_boot = v_para_boot * deltat_23

# plot_velocidades(X1, Y1, Z1, R, normal_fit, x3, v_media, v_para, v_para_MVA)

###########
#giroradio


#########
###análisis de corrientes

inicio_up = np.where(t == find_nearest_inicial(t, t1-0.015))[0][0] #las 18:12:00
fin_up = np.where(t == find_nearest_final(t, t1))[0][0] #las 18:13:00
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0) #nT

inicio_down = np.where(t == find_nearest_inicial(t, t4))[0][0] #las 18:14:51
fin_down = np.where(t == find_nearest_final(t, t4+0.015))[0][0] #las 18:15:52
B_downstream = np.mean(B[inicio_down:fin_down,:], axis=0) #nT



omega = np.arccos(np.dot(B_upstream,B_downstream)/(np.linalg.norm(B_upstream)*np.linalg.norm(B_downstream)))

mu = 4* np.pi * 1E-7 #Henry/m

J_s_MVA = np.cross(x3, (B_upstream-B_downstream)) / mu #nA/m
ancho_mpb = np.linalg.norm(x_23_MVA) #considero que es el tiempo corto a lo largo de la normal del MVA
J_v_MVA = J_s_MVA / (1000*ancho_mpb) #nA/m²

J_s_fit = np.cross(normal_fit, (B_upstream-B_downstream)) / mu #nA/m
ancho_mpb_fit = np.linalg.norm(x_23_fit) #considero que es el tiempo corto a lo largo de la normal del ajuste
J_v_fit = J_s_fit / (1000*ancho_mpb) #nA/m²

J_s_boot = np.cross(normal_boot, (B_upstream-B_downstream)) / mu #nA/m
ancho_mpb_boot = np.linalg.norm(x_23_boot) #km
J_v_boot = J_s_boot / (1000*ancho_mpb_boot) #nA/m²

fuerza_mva = np.cross(J_v_MVA * 1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
fuerza_fit = np.cross(J_v_fit *1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
fuerza_boot = np.cross(J_v_boot * 1E-9, B[inicio_down,:]*1E-9) #N/m^3


e_density = lpw[:,-1]
t_lpw = lpw[:,3] + lpw[:,4]/60 + lpw[:,5]/3600

ti_lpw = np.where(t_lpw == find_nearest(t_lpw, t2))[0][0]
tf_lpw = np.where(t_lpw == find_nearest(t_lpw, t3))[0][0]
n_e = np.nanmean(e_density[ti_lpw:tf_lpw]) #hace el mean ignorando los nans #cm⁻³
n_e = n_e * 1E6 #m⁻³
if np.isnan(n_e):
    n_e = 1E7
    print('LPW no tiene datos de densidad, asumí n_e = 1E7')
q_e = 1.6E-19 #carga electron #C

E_Hall = np.cross(J_v_MVA * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
E_Hall_fit = np.cross(J_v_fit * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
E_Hall_boot = np.cross(J_v_boot * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m


#
# plot_FLorentz(X1, Y1, Z1, R, J_v, B_upstream, B_downstream, fuerza_mva, x3)
#
# plt.figure()
# plt.plot(t[inicio_up:fin_down], fuerza_mva)
# plt.plot(t[inicio_up:fin_down], fuerza_ajuste)
# plt.legend(['Fx', 'Fy', 'Fz'])
# plt.xlabel('Tiempo')
# plt.ylabel('Fuerza (N/m^3)')
#
plt.show(block=False)
