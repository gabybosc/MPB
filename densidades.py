import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.mlab import normpdf
from scipy.stats import norm
import spacepy.pycdf as cdf
from funciones import find_nearest, unix_to_decimal, plot_select, set_axes_equal

np.set_printoptions(precision=4)


###########DATOS
path = 'datos/marzo 2016/16/'
cdf_swia = cdf.CDF(path + 'mvn_swi_l2_onboardsvymom_20160316_v01_r01.cdf')
datos = np.loadtxt(path + 'mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148)
n =2
datos = datos[:-n, :]

t_unix = np.asarray(cdf_swia['time_unix'][...]) #no olvidarse los ...!
density = np.asarray(cdf_swia['density'][...])
temperature = np.asarray(cdf_swia['temperature_mso'][...])
vel_mso_xyz = np.asarray(cdf_swia['velocity_mso'][...])

t_swia = unix_to_decimal(t_unix) 
inicio = np.where(t_swia == find_nearest(t_swia, 17.9))[0][0]
fin = np.where(t_swia == find_nearest(t_swia, 18.4))[0][0]

t_swia_cut = t_swia[inicio:fin]
density_cut = density[inicio:fin]
temperature_cut = temperature[inicio:fin]
vel_mso_xyz = vel_mso_xyz[inicio:fin] #km/s

#los datos del campo B
dia = datos[:,1]
t_sts_raw = datos[:,6]
t_sts_raw = (t_sts_raw - dia) * 24
B = np.empty((len(datos), 3))
for i in range(7,10):
    B[:,i-7] = datos[:, i]

posicion = np.empty((len(datos), 3))
for i in range(11,14):
    posicion[:,i-11] = datos[:, i]


#quiero diezmar el tiempo y el campo para que sean los mismos que tiene swia
idx = np.empty(449)
for i in range(449):
    idx[i] = np.where(t_sts_raw == find_nearest(t_sts_raw, t_swia_cut[i]))[0][0]
idx = idx.astype(int)

t_sts = t_sts_raw[idx] #lo diezmó

B_cut = B[idx]
posicion_cut = posicion[idx]


ti_swia = np.where(t_swia == find_nearest(t_swia, 18.05))[0][0]
tf_swia = np.where(t_swia == find_nearest(t_swia, 18.18))[0][0]

ti_mag = np.where(t_sts == find_nearest(t_sts, 18.1))[0][0]
tf_mag = np.where(t_sts == find_nearest(t_sts, 18.2))[0][0]

inicio_mpb = np.where(t_sts == find_nearest(t_sts, 18.2156))[0][0]
fin_mpb = np.where(t_sts == find_nearest(t_sts, 18.2506))[0][0]
ancho_mpb = np.linalg.norm(posicion_cut[inicio_mpb] - posicion_cut[fin_mpb]) #km

####################

density_mean = np.empty(tf_swia-ti_swia)
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
    v_maven[i,:] = (posicion_cut[i+1,:] - posicion_cut[i,:]) / (t_sts[i+1]-t_sts[i]) / 3600 #en km/s
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

E_Hall = np.empty((len(idx),3))
for i in range(len(idx)):
    E_Hall[i,:] = E_convective[i,:] * (v_alfven/vel_mso)[i] * ion_length / ancho_mpb #nV/km


plt.figure()
for xc in [18.2193,18.2272,18.2331,18.2476]:
    plt.axvline(x = xc, color = 'black', linewidth=0.5)
plt.plot(t_swia_cut, density_cut, label='density')
plt.scatter(t_swia[ti_swia+15:tf_swia+15], density_mean, color = 'r', label='density avg')
plt.ylabel('ion density (cm⁻³)')
plt.xlabel('time (hdec)')
plt.legend()

plt.figure()
plt.plot(t_sts, temp_para, label=r'T$\parallel$')
plt.plot(t_sts, temp_perp, label=r'T $\perp$')
for xc in [18.2193,18.2272,18.2331,18.2476]:
    plt.axvline(x = xc, color = 'black', linewidth=0.5)
plt.ylabel('temperature (eV)')
plt.legend()

plt.figure()
plt.plot(t_sts, v_alfven, label='vel Alfven')
plt.plot(t_sts, vel_mso, label='vel mso')
for xc in [18.2193,18.2272,18.2331,18.2476]:
    plt.axvline(x = xc, color = 'black', linewidth=0.5)
plt.ylabel('Velocity (km/s)')
plt.legend()

plt.figure()
plt.plot(t_sts[ti_mag:tf_mag], ion_gyroradius)
plt.axhline(y = np.mean(ion_gyroradius, axis=0), color = 'r', label='mean')
plt.ylabel('Ion gyroradius (km)')
plt.legend()

plt.figure()
plt.axvline(x = 18.2193, color = 'black', linewidth=0.5, label ='MPB')
plt.plot(t_sts[0:inicio_mpb], E_convective_normalized[0:inicio_mpb], label='E convectivo')
plt.plot(t_sts[0:inicio_mpb], np.linalg.norm(E_Hall, axis=1)[0:inicio_mpb], label='E Hall')
plt.ylabel('E (nV/km)')
plt.legend()


plt.show(block = False)
