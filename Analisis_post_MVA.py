import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
from matplotlib.mlab import normpdf
from scipy.stats import norm
import datetime as dt
from funciones import error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, datenum, unix_to_decimal
from funciones_MVA import ajuste_conico, plot_velocidades, plot_FLorentz
from funciones_plot import hodograma, set_axes_equal

"""
Toma las normales y lo que devuelva MVA.py (tanto hires como lowres) y ahce un análisis con esto.
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
"""
###corremos el MVA e importamos las variables que quiero (las normales)
import MVA_lowres
from MVA_lowres import *


#########
#buscamos el ángulo entre las normales
angulo_mva = np.arccos(np.clip(np.dot(normal_fit, x3), -1.0, 1.0)) #el clip hace que si por algun motivo el dot me da >1 (i.e. 1,00002), me lo convierte en 1
angulo_boot = np.arccos(np.clip(np.dot(normal_boot, x3), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!

print('El ángulo entre la normal del ajuste y la del MVA = {0:1.3g}º'.format(angulo_mva * 180/np.pi))
print('El ángulo entre la normal del bootstrap y la del MVA = {0:1.3g}º'.format(angulo_boot * 180/np.pi))


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
print('La velocidad minima es = {0:1.3g} km/s, la diferencia entre la minima y la maxima es = {1:1.3g} km/s'.format(min(norma_v), diff))
#la velocidad promedio
v_media = np.array([np.mean(v_punto[:,0]), np.mean(v_punto[:,1]), np.mean(v_punto[:,2])])
print('La velocidad media de la nave a través de la MPB es = {0:1.3g} km/s'.format(np.linalg.norm(v_media)))

#ahora quiero ver si la nave atraviesa perpendicularmente a la MPB
v_media_norm = v_media/np.linalg.norm(v_media) #la normalizo por cuestion de  cuentas nomas
angulo_v_ajuste = np.arccos(np.clip(np.dot(normal_fit, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_v_mva = np.arccos(np.clip(np.dot(x3, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_v_boot = np.arccos(np.clip(np.dot(normal_boot, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre la velocidad y la normal del ajuste es = {0:1.3g}º'.format(angulo_v_ajuste * 180/np.pi))
print('El ángulo entre la velocidad y la normal del MVA es = {0:1.3g}º'.format(angulo_v_mva * 180/np.pi))
print('El ángulo entre la velocidad y la normal del bootstrap es = {0:1.3g}º'.format(angulo_v_boot * 180/np.pi))

B_intermedio = B_medio_vectorial/np.linalg.norm(B_medio_vectorial) #la normalizo por cuestion de  cuentas nomas
angulo_B_ajuste = np.arccos(np.clip(np.dot(normal_fit, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_B_mva = np.arccos(np.clip(np.dot(x3, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_B_boot = np.arccos(np.clip(np.dot(normal_boot, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre el campo B y la normal del ajuste es = {0:1.3g}º'.format(angulo_B_ajuste * 180/np.pi))
print('El ángulo entre el campo B y la normal del MVA es = {0:1.3g}º'.format(angulo_B_mva * 180/np.pi))
print('El ángulo entre el campo B y la normal del bootstrap es = {0:1.3g}º'.format(angulo_B_boot * 180/np.pi))

# ######
# ##Espesor de la MPB
# #ahora veamos v_para
# v_para = np.dot(v_media, normal_fit) * normal_fit
# deltat_14 = (t4 - t1) * 3600
# deltat_23 = (t3 - t2) * 3600
# x_14_fit = v_para * deltat_14 #en km# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# x_23_fit = v_para * deltat_23
# print('El ancho de la MPB proyectando sobre la normal de Vignes para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14_fit), np.linalg.norm(x_23_fit)))
#
# #si ahora proyecto sobre la normal de la MVA
# v_para_MVA = np.dot(v_media, x3) * x3
# x_14_MVA = v_para_MVA * deltat_14 #en km
# x_23_MVA = v_para_MVA * deltat_23
# print('El ancho de la MPB proyectando sobre la normal del MVA para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14_MVA), np.linalg.norm(x_23_MVA)))
#
# #si ahora proyecto sobre la normal del bootstrap
# v_para_boot = np.dot(np.array([-2.487,  0.479,  2.836]), normal_boot) * normal_boot
# x_14_boot = v_para_boot * deltat_14 #en km
# x_23_boot = v_para_boot * deltat_23
# print('El ancho de la MPB proyectando sobre la normal del bootstrap para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14_boot), np.linalg.norm(x_23_boot))
#
# for i in [t_1,t_2,t_3,t_4]:
#     sza = np.arccos(np.clip(np.dot(posicion[i]/np.linalg.norm(posicion[i]), [1,0,0]), -1.0, 1.0))
#     altitud = MD[i, 8]
#     print('El tiempo = {2:1.6g} tiene SZA = {0:1.3g}º y altitud = {1:1.3g} km '.format(sza * 180/np.pi, altitud, t[i]))
#
#
#
# plot_velocidades(X1, Y1, Z1, R, normal_fit, x3, v_media, v_para, v_para_MVA)
#
#
# #########
# ###análisis de corrientes
#
# inicio_up = np.where(t == find_nearest_inicial(t, t1-0.015))[0][0] #las 18:12:00
# fin_up = np.where(t == find_nearest_final(t, t1))[0][0] #las 18:13:00
# B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0) #nT
#
# inicio_down = np.where(t == find_nearest_inicial(t, t4))[0][0] #las 18:14:51
# fin_down = np.where(t == find_nearest_final(t, t4+0.015))[0][0] #las 18:15:52
# B_downstream = np.mean(B[inicio_down:fin_down,:], axis=0) #nT
#
# print('B upstream es {0:1.3g} nT y su módulo es {0:1.3g}'.format(B_upstream, np.linalg.norm(B_upstream)))
# print('B downstream es {0:1.3g} nT y su módulo es {0:1.3g}'.format(B_downstream, np.linalg.norm(B_downstream)))
#
# omega = np.arccos(np.dot(B_upstream,B_downstream)/(np.linalg.norm(B_upstream)*np.linalg.norm(B_downstream)))
# print('El ángulo omega es {0:1.3g}º'.format(omega * 180/np.pi))
#
# mu = 4* np.pi * 1E-7 #Henry/m
#
# J_s = np.cross(x3, (B_upstream-B_downstream)) / mu #nA/m
#
# ancho_mpb = np.linalg.norm(x_23_MVA) #considero que es el tiempo corto a lo largo de la normal del MVA
# J_v = J_s / (1000*ancho_mpb) #nA/m²
#
# print('La corriente superficial con la normal del MVA es Js = {0:1.3g} mA/m, |Js| = {0:1.3g} mA/m'.format(J_s *1E-6, np.linalg.norm(J_s)*1E-6))
# print('La corriente en volumen con la normal del MVA es Jv = {0:1.3g} nA/m², |Jv| = {0:1.3g} nA/m²'.format(J_v, np.linalg.norm(J_v)))
#
#
# J_s_ajuste = np.cross(normal_fit, (B_upstream-B_downstream)) / mu #nA/m
# ancho_mpb_ajuste = np.linalg.norm(x_23_fit) #considero que es el tiempo corto a lo largo de la normal del ajuste
# J_v_ajuste = J_s_ajuste / (1000*ancho_mpb) #nA/m²
# print('La corriente superficial con la normal del ajuste es Js = {0:1.3g} mA/m, |Js| = {0:1.3g} mA/m'.format(J_s_ajuste *1E-6, np.linalg.norm(J_s_ajuste)*1E-6))
# print('La corriente en volumen con la normal del ajuste es Jv = {0:1.3g} nA/m², |Jv| = {0:1.3g} nA/m²'.format(J_v_ajuste, np.linalg.norm(J_v_ajuste)))
#
#
# # fuerza_mva = np.cross(J_v * 1E-9, B[inicio_up:fin_down, :] * 1E-9) #A/m² * T = N/m³
# # fuerza_ajuste = np.cross(J_v_ajuste * 1E-9, B[inicio_up:fin_down, :] * 1E-9) #A/m² * T = N/m³
# fuerza_mva = np.cross(J_v * 1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
# fuerza_ajuste = np.cross(J_v_ajuste * 1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
# print('La fuerza de lorentz del MVA es {0:1.3g} V/m, su magnitud es {0:1.3g}'.format(fuerza_mva, np.linalg.norm(fuerza_mva)))
# print('La fuerza de lorentz del ajuste es {0:1.3g} V/m, su magnitud es {0:1.3g}'.format(fuerza_ajuste, np.linalg.norm(fuerza_ajuste)))
#
#
#
# e_density = lpw.varget('data')[:,3]
# t_unix = lpw.varget('time_unix')
# t_lpw = unix_to_decimal(t_unix)
# ti_lpw = np.where(t_lpw == find_nearest(t_lpw, t2))[0][0]
# tf_lpw = np.where(t_lpw == find_nearest(t_lpw, t3))[0][0]
# n_e = np.nanmean(e_density[ti_lpw:tf_lpw]) #hace el mean ignorando los nans #cm⁻³
# n_e = n_e * 1E6 #m⁻³
# # n_e = 1E7
# q_e = 1.6E-19 #carga electron #C
#
# E_Hall = np.cross(J_v * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
# print('El campo de Hall del MVA es {} mV/m, su magnitud es {}'.format(E_Hall*1E3, np.linalg.norm(E_Hall)*1E3))
#
# E_Hall_ajuste = np.cross(J_v_ajuste * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
# print('El campo de Hall del ajuste es {} mV/m, su magnitud es {}'.format(E_Hall_ajuste*1E3, np.linalg.norm(E_Hall_ajuste)*1E3))
#
# J_s_boot = np.cross(normal_boot, (B_upstream-B_downstream)) / mu #nA/m
#
# ancho_mpb_boot = 77 #km, esta copiado de lo que da en el otro script
# J_v_boot = J_s_boot / (1000*ancho_mpb_boot) #nA/m²
# print('La corriente superficial con la normal del bootstrap es Js = {} mA/m, |Js| = {} mA/m'.format(J_s_boot *1E-6, np.linalg.norm(J_s_boot)*1E-6))
# print('La corriente en volumen con la normal del bootstrap es Jv = {} nA/m², |Jv| = {} nA/m²'.format(J_v_boot, np.linalg.norm(J_v_boot)))
#
# fuerza_mva_boot = np.cross(J_v_boot * 1E-9, B[inicio_down,:]*1E-9) #N/m^3
# print('La fuerza de lorentz es {} V/m, su magnitud es {}'.format(fuerza_mva_boot, np.linalg.norm(fuerza_mva_boot)))
# E_Hall_boot = np.cross(J_v_boot * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
# print('El campo de Hall es {} mV/m, su magnitud es {}'.format(E_Hall_boot*1E3, np.linalg.norm(E_Hall_boot)*1E3))
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
# plt.show(block=False) #para que siga andando aunque no lo haya cerrado
