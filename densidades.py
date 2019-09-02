import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
from funciones import find_nearest, unix_to_decimal

np.set_printoptions(precision=4)

"""
debuggear
cambiar para que use cdflib
"""

###########DATOS
path = '../../datos/'
swia = cdf.CDF(path + 'SWIA/mvn_swi_l2_onboardsvymom_20160316_v01_r01.cdf')
mag = np.loadtxt(path + 'MAG_1s/2016/mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148)
n =2
mag = mag[:-n, :]

t1 = 18.2167
t2 = 18.2204
t3 = 18.235
t4 = 18.2476

t_unix = swia.varget('time_unix')
density = swia.varget('density')
temperature = swia.varget('temperature_mso')
vel_mso_xyz = swia.varget('velocity_mso')

t_swia = unix_to_decimal(t_unix)
inicio = np.where(t_swia == find_nearest(t_swia, 17.9))[0][0]
fin = np.where(t_swia == find_nearest(t_swia, 18.4))[0][0]

t_swia_cut = t_swia[inicio:fin]
density_cut = density[inicio:fin]
temperature_cut = temperature[inicio:fin]
vel_mso_xyz = vel_mso_xyz[inicio:fin] #km/s

#los datos del campo B
dia = mag[:,1]
t_mag = mag[:,6]
t_mag = (t_mag - dia) * 24
B = np.empty((len(mag), 3))
for i in range(7,10):
    B[:,i-7] = mag[:, i]

posicion = np.empty((len(mag), 3))
for i in range(11,14):
    posicion[:,i-11] = mag[:, i]


#quiero diezmar el tiempo y el campo para que sean los mismos que tiene swia
idx = np.zeros(len(t_swia_cut))
for i in range(len(idx)):
    idx[i] = np.where(t_mag == find_nearest(t_mag, t_swia_cut[i]))[0][0]
idx = idx.astype(int)

t_diezmado = t_mag[idx] #lo diezmó

B_cut = B[idx]
posicion_cut = posicion[idx]

####Qué son esos tiempos??
ti_swia = np.where(t_swia == find_nearest(t_swia, 18.05))[0][0]
tf_swia = np.where(t_swia == find_nearest(t_swia, 18.18))[0][0]

ti_mag = np.where(t_diezmado == find_nearest(t_diezmado, 18.1))[0][0]
tf_mag = np.where(t_diezmado == find_nearest(t_diezmado, 18.2))[0][0]

inicio_mpb = np.where(t_diezmado == find_nearest(t_diezmado, 18.2156))[0][0]
fin_mpb = np.where(t_diezmado == find_nearest(t_diezmado, 18.2506))[0][0]
ancho_mpb = np.linalg.norm(posicion_cut[inicio_mpb] - posicion_cut[fin_mpb]) #km

####################

density_mean = np.empty(tf_swia-ti_swia) #upstream
for i in range(tf_swia-ti_swia):
    density_mean[i] = np.mean(density[ti_swia+i:ti_swia+i+30])

#Elijo la densidad promedio media
idx_d = int(np.abs(np.argmin(density_mean) + np.argmax(density_mean))/2)
density_avg = density_mean[idx_d] #cm

ion_length = 2.28E07 / np.sqrt(density_avg) * 1E-5 #km
print('Ion inertial length = {0:1.3g} km'.format(ion_length))


B_avg = np.empty((len(idx), 3))
B_avg_normalized = np.empty((len(idx),3))
temp_para_xyz = np.empty((len(idx), 3))
v_maven = np.empty((len(idx),3))
v_planet = np.empty((len(idx),3))


for i in range(len(idx)-1):
    B_avg[i,:] = np.mean(B_cut[i:i+30,:], axis=0) * 1E-5 #lo paso a gauss
    B_avg_normalized[i,:] = B_avg[i,:] / np.linalg.norm(B_avg[i,:]) #adimensional
    temp_para_xyz[i,:] = np.dot(B_avg_normalized[i,:], temperature_cut[i,:]) * B_avg_normalized[i,:] #eV
    v_maven[i,:] = (posicion_cut[i+1,:] - posicion_cut[i,:]) / (t_diezmado[i+1]-t_diezmado[i]) / 3600 #en km/s
    v_maven[-1,:] = v_maven[-2,:] #porque si no, me queda vacío
    v_planet[i,:] = np.nanmean(vel_mso_xyz[i:i+30,:] + v_maven[i:i+30,:], axis=0) #km/s
    v_planet[-1,:] = v_planet[-2,:] #porque si no, me queda vacío

temp_para = np.linalg.norm(temp_para_xyz, axis = 1) #eV
temp_perp = np.linalg.norm(temperature_cut - temp_para_xyz, axis=1) #eV
v_alfven = np.linalg.norm(2.18E11 * B_avg / np.sqrt(density_avg), axis=1) * 1E-5 #km/s
vel_mso = np.linalg.norm(vel_mso_xyz, axis = 1) #km/s


ion_gyroradius = np.empty(tf_mag-ti_mag)
for i in range(tf_mag-ti_mag):
    ion_gyroradius[i] = 1.02E02 * np.sqrt(temp_perp[ti_mag+i])/np.linalg.norm(B_avg[ti_mag+i,:]) * 1E-5#km

print('ion gyroradius mean = {0:1.3g} km'.format(np.nanmean(ion_gyroradius, axis=0))) #nanmean ignora los nans

print('v_alfven / v_mso = {0:1.3g} '.format(np.mean(v_alfven[ti_mag:tf_mag] / vel_mso[ti_mag:tf_mag])))

E_convective = np.cross(-v_planet, B_avg)
E_convective_normalized = np.linalg.norm(E_convective, axis=1) #nV/km

print('El E convectivo medio es {0:1.3g} nV/km'.format(np.mean(E_convective_normalized[ti_mag:tf_mag]))) #en SI nV/km


#########presión térmica
ti_up = np.where(t_swia == find_nearest(t_swia, t1-0.15))[0][0]
tf_up = np.where(t_swia == find_nearest(t_swia, t1))[0][0]
ti_down = np.where(t_swia == find_nearest(t_swia, t4))[0][0]
tf_down = np.where(t_swia == find_nearest(t_swia, t4+0.15))[0][0]

density_up = np.empty(tf_up-ti_up) #upstream
for i in range(tf_up-ti_up):
    density_up[i] = np.mean(density[ti_up+i:tf_up+i+30]) * 1E6

density_down = np.empty(tf_down-ti_down) #downstream
for i in range(tf_down-ti_down):
    density_down[i] = np.mean(density[ti_down+i:tf_down+i+30]) * 1E6


#densidad en kg/m^3
mp = 1.67e-27 #masa del proton en kg
#mp = 1.5e-10 #masa del proton en joules/c^2
rho_u = mp*density_up
rho_d = mp*density_down

#presion suponiendo gas ideal (en Pa=J/m^3)
kB = 1.38e-23 #cte de Boltzmann en J/K
#por ahora supongo T = 2*Ti, tendria que ser T=Ti+Te
temperatura_swia_norm = np.linalg.norm(temperature, axis=1)
T_up = 2*np.mean(temperatura_swia_norm[ti_up:tf_up])*(11604.5) #en K
# T_down = 2*np.mean(temperatura_static_norm[ti_down:tf_down])*(11604.5) #en K ###en realidad como no tengo casi prptones del sw de este lado debería usar los datos de STATIC en lugar de swia. O buscar valor.
T_down = 0.5 * 11604.5 #K
##con Te estimadas de las distribuciones
P_up = density_up*kB*T_up
P_down = density_down*kB*T_down

#######Presión magnética:
inicio_up = np.where(t_diezmado == find_nearest_inicial(t_diezmado, t1-0.015))[0][0] #las 18:12:00
fin_up = np.where(t_diezmado == find_nearest_final(t_diezmado, t1))[0][0] #las 18:13:00
B_upstream = np.mean(B_cut[inicio_up:fin_up, :], axis=0) #nT

inicio_down = np.where(t_diezmado == find_nearest_inicial(t_diezmado, t4))[0][0] #las 18:14:51
fin_down = np.where(t_diezmado == find_nearest_final(t_diezmado, t4+0.015))[0][0] #las 18:15:52
B_downstream = np.mean(B_cut[inicio_down:fin_down,:], axis=0) #nT

#######cociente de presiones:
mu = (np.pi*4)*(1e-7) #permeabilidad mag del vacio en N/A² =  T m/A

P_Bup = np.linalg.norm(B_upstream*1E-9)**2/(2*mu) # T A/ m = Pa
P_Bdown = np.linalg.norm(B_downstream*1E-9)**2/(2*mu)

#beta del plasma
beta_up = P_up/P_Bup
beta_down = P_down / P_Bdown


plt.figure()
for xc in [18.2193,18.2272,18.2331,18.2476]:
    plt.axvline(x = xc, color = 'black', linewidth=0.5)
plt.plot(t_swia_cut, density_cut, label='density')
plt.axvline(x = t_diezmado[ti_mag], color ='red')
plt.axvline(x = t_diezmado[tf_mag], color ='red')
plt.axvline(x = t_swia[ti_swia], color ='blue')
plt.axvline(x = t_swia[tf_swia], color ='blue')
plt.scatter(t_swia[ti_swia+15:tf_swia+15], density_mean, color = 'r', label='density avg')
plt.ylabel('ion density (cm⁻³)')
plt.xlabel('time (hdec)')
plt.legend()

plt.figure()
plt.plot(t_diezmado, temp_para, label=r'T$\parallel$')
plt.plot(t_diezmado, temp_perp, label=r'T $\perp$')
plt.axvline(x = t_diezmado[ti_mag], color ='red')
plt.axvline(x = t_diezmado[tf_mag], color ='red')
plt.axvline(x = t_swia[ti_swia], color ='blue')
plt.axvline(x = t_swia[tf_swia], color ='blue')
for xc in [18.2193,18.2272,18.2331,18.2476]:
    plt.axvline(x = xc, color = 'black', linewidth=0.5)
plt.ylabel('temperature (eV)')
plt.legend()

plt.figure()
plt.plot(t_diezmado, v_alfven, label='vel Alfven')
plt.plot(t_diezmado, vel_mso, label='vel mso')
plt.axvline(x = t_diezmado[ti_mag], color ='red')
plt.axvline(x = t_diezmado[tf_mag], color ='red')
plt.axvline(x = t_swia[ti_swia], color ='blue')
plt.axvline(x = t_swia[tf_swia], color ='blue')
for xc in [18.2193,18.2272,18.2331,18.2476]:
    plt.axvline(x = xc, color = 'black', linewidth=0.5)
plt.ylabel('Velocity (km/s)')
plt.legend()

plt.figure()
plt.plot(t_diezmado[ti_mag:tf_mag], ion_gyroradius)
plt.axvline(x = t_diezmado[ti_mag], color ='red')
plt.axvline(x = t_diezmado[tf_mag], color ='red')
plt.axvline(x = t_swia[ti_swia], color ='blue')
plt.axvline(x = t_swia[tf_swia], color ='blue')
plt.axhline(y = np.mean(ion_gyroradius, axis=0), color = 'r', label='mean')
plt.ylabel('Ion gyroradius (km)')
plt.legend()

plt.figure()
plt.axvline(x = 18.2193, color = 'black', linewidth=0.5, label ='MPB')
plt.axvline(x = t_diezmado[ti_mag], color ='red')
plt.axvline(x = t_diezmado[tf_mag], color ='red')
plt.axvline(x = t_swia[ti_swia], color ='blue')
plt.axvline(x = t_swia[tf_swia], color ='blue')
plt.plot(t_diezmado[0:inicio_mpb], E_convective_normalized[0:inicio_mpb], label='E convectivo')
plt.ylabel('E (nV/km)')
plt.legend()


plt.show(block = False)
