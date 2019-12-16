import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
from funciones import find_nearest, unix_to_decimal, find_nearest_inicial, find_nearest_final, fechas

np.set_printoptions(precision=4)

"""
NO ANDA, NECESITA TEMPERATURA Y MÁS COSAS
En el cálculo de beta = P_termica / P_mag no puedo usar los datos de SWIA para calcular P_t downstream, ya que no hay datos de viento solar en esa zona.
"""

##########CONSTANTES
mp = 1.67e-27 #masa del proton en kg
kB = 1.38e-23 #cte de Boltzmann en J/K
#mp = 1.5e-10 #masa del proton en joules/c^2
q_e = 1.602e-19 #carga del electrón en C

###########DATOS
year, month, day, doy = fechas()
ti = input('Hora del cruce (HH)\n')

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
swia = np.loadtxt(path+'SWIA.asc')
mag = np.loadtxt(path + 'mag_filtrado.txt', skiprows=2)

datos = np.loadtxt('../outputs/t1t2t3t4.txt')
for j in range(len(datos)):
    if datos[j,0] == float(year) and datos[j,1] == float(doy) and int(datos[j,2]) == int(ti):
        i = j

t1 = datos[i,2]
t2 = datos[i,3]
t3 = datos[i,4]
t4 = datos[i,5]

# temperature = swia.varget('temperature_mso') #eV
# vel_mso_xyz = swia.varget('velocity_mso') #km/s

t_swia = swia[:,3] + swia[:,4]/60 + swia[:,5]/3600 #hdec
density = swia[:,-1]
inicio = np.where(t_swia == find_nearest(t_swia, 17.9))[0][0]
fin = np.where(t_swia == find_nearest(t_swia, 18.4))[0][0]

t_swia_cut = t_swia[inicio:fin]
density_cut = density[inicio:fin]

#los datos del campo B
M = len(mag[:,0]) #el numero de datos
B = mag[:, :3]
Bnorm = mag[:,-1]

mag = np.loadtxt(path + 'MAG.asc')
hh = mag[:,3]
mm = mag[:,4]
ss = mag[:,5]

t_mag = hh + mm/60 + ss/3600 #hdec

posicion = np.zeros((M, 3))
for i in range(9,12):
    posicion[:,i-9] = mag[:, i]

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

inicio_mpb = np.where(t_diezmado == find_nearest(t_diezmado, t1))[0][0]
fin_mpb = np.where(t_diezmado == find_nearest(t_diezmado, t4))[0][0]
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


corriente_comp = np.array(np.transpose([q_e * density_cut * vel_mso_xyz[:,columnas] * 1E9 for columnas in [0,1,2]]))#A/m²
corriente_norma = q_e * vel_mso * density_cut * 1E9 #A/m²

plt.figure()
plt.plot(t_diezmado,corriente_norma*1E6)
plt.xlabel('Tiempo')
plt.ylabel('|j| (mA/m²)')

plt.figure()
plt.plot(t_diezmado,corriente_comp*1E6)
plt.xlabel('Tiempo')
plt.ylabel('Corriente (mA/m²)')
plt.legend(['jx', 'jy', 'jz'])


ion_gyroradius = np.empty(tf_mag-ti_mag)
for i in range(tf_mag-ti_mag):
    ion_gyroradius[i] = 1.02E02 * np.sqrt(temp_perp[ti_mag+i])/np.linalg.norm(B_avg[ti_mag+i,:]) * 1E-5#km

print('ion gyroradius mean = {0:1.3g} km'.format(np.nanmean(ion_gyroradius, axis=0))) #nanmean ignora los nans

print('v_alfven / v_mso = {0:1.3g} '.format(np.mean(v_alfven[ti_mag:tf_mag] / vel_mso[ti_mag:tf_mag])))

E_convective = np.cross(-v_planet, B_avg)
E_convective_normalized = np.linalg.norm(E_convective, axis=1) #nV/km

print('El E convectivo medio es {0:1.3g} nV/km'.format(np.mean(E_convective_normalized[ti_mag:tf_mag]))) #en SI nV/km




plt.figure()
for xc in [t1,t2,t3,t4]:
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
for xc in [t1,t2,t3,t4]:
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
for xc in [t1,t2,t3,t4]:
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
