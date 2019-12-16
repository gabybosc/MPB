import numpy as np
import cdflib as cdf
from funciones import fechas, find_nearest_final, find_nearest_inicial, unix_to_decimal


qp = 1.6e-19 #Coulombs (q electrica del proton)
mp = 1.67e-27 #kg (masa del proton)

year, month, day, doy = fechas()
ti = input('Hora del cruce (HH)\n')
in_out = input('Inbound? (y/n)\n')

path = '../../datos/'
swia = cdf.CDF(path + f'SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf')
mag = np.loadtxt(path + f'MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=160)

datos = np.loadtxt('outputs/t1t2t3t4.txt')
for j in range(len(datos)):
    if datos[j,0] == float(year) and datos[j,1] == float(doy) and int(datos[j,2]) == int(ti):
        i = j

if in_out=='y':
    t1 = datos[i,2]
    t2 = datos[i,3]
    t3 = datos[i,4]
    t4 = datos[i,5]
else:
    t4 = datos[i,2]
    t3 = datos[i,3]
    t2 = datos[i,4]
    t1 = datos[i,5]


t_unix = swia.varget('time_unix')
t_swia = unix_to_decimal(t_unix)
density = swia.varget('density') #cm⁻³
temperature = swia.varget('temperature_mso') #eV
vel_mso_xyz = swia.varget('velocity_mso') #km/s
vel_mso = np.linalg.norm(vel_mso_xyz, axis = 1) #km/s

dia = mag[:,1]
t_mag = mag[:,6]
t_mag = (t_mag - dia) * 24
B = np.empty((len(mag), 3))
for i in range(7,10):
    B[:,i-7] = mag[:, i]



inicio_up = np.where(t_swia == find_nearest_inicial(t_swia, t1-0.015))[0][0]
fin_up = np.where(t_swia == find_nearest_final(t_swia, t1))[0][0]
V_upstream = np.mean(vel_mso[inicio_up:fin_up])

bi = np.where(t_mag == find_nearest_inicial(t_mag, t1-0.015))[0][0]
bf = np.where(t_mag == find_nearest_inicial(t_mag, t1))[0][0]
B_upstream = np.mean(B[bi:bf, :], axis=0) #nT
Bu = np.linalg.norm(B_upstream)

density_mean = np.zeros(fin_up-inicio_up)
for i in range(fin_up-inicio_up):
    density_mean[i] = np.mean(density[inicio_up+i:inicio_up+i+30])

#Elijo la densidad promedio media
idx_d = int(np.abs(np.argmin(density_mean) + np.argmax(density_mean))/2)
density_avg = density_mean[idx_d] #cm

# density_avg es valor medio de la densidad numerica upstream (solar wind) en unidades de 1/cm^3
di = 2.28e2 / np.sqrt(density_avg)

print(f'La longitud inercial de protones del SW es {di:.3g} km')


# W_ci la frecuencia de giroradio de los iones upstream en Hz
#Bu es el modulo del campo medio upstream
W_ci = qp*Bu*(1e-9)/mp #rad/s (pase B de nT a T)

ri = V_upstream / W_ci # giroradio convectivo iones en km
# V_upstream es la velocidad del viento solar upstream en la componente normal en km/s

print(f'El giroradio es {ri:.3g} km')
