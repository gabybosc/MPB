import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.mlab import normpdf
from scipy.stats import norm
import spacepy.pycdf as cdf
import matplotlib.dates as md
import datetime as dt
from funciones import find_nearest, unix_to_decimal, plot_select, set_axes_equal
from datetime import datetime

np.set_printoptions(precision=4)

'''
Este código calcula EHall a partir del E convectivo
'''

###########DATOS
path = 'datos/marzo 2016/16/'
cdf_swia = cdf.CDF(path + 'mvn_swi_l2_onboardsvymom_20160316_v01_r01.cdf')
datos = np.loadtxt(path + 'mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148)
n =2
datos = datos[:-n, :]

density = cdf_swia['density'] #cgs -  1/cm3
time_unix = cdf_swia['time_unix'] #s
v_mso_xyz = cdf_swia['velocity_mso'] #SI - km/s

t_unix = np.asarray(time_unix[...]) #no olvidarse los ...!
density = np.asarray(density[...]) #cgs
v_mso_imported = np.asarray(v_mso_xyz[...]) #SI

t_swia = unix_to_decimal(t_unix)
inicio = np.where(t_swia == find_nearest(t_swia, 17.9))[0][0]
fin = np.where(t_swia == find_nearest(t_swia, 18.4))[0][0]

t_swia_cut = t_swia[inicio:fin]
density_cut = density[inicio:fin]
v_mso_xyz = v_mso_imported[inicio:fin]

#los datos del campo B
dia = datos[:,1]
t_sts = datos[:,6]
t_sts = (t_sts - dia) * 24
B = np.empty((len(datos), 3))
for i in range(7,10):
    B[:,i-7] = datos[:, i] #nT

posicion = np.empty((len(datos), 3))
for i in range(11,14):
    posicion[:,i-11] = datos[:, i] #km


#quiero diezmar el tiempo y el campo para que sean los mismos que tiene swia
idx = np.empty(449)
for i in range(449):
    idx[i] = np.where(t_sts == find_nearest(t_sts, t_swia_cut[i]))[0][0]
idx = idx.astype(int)

t_sts = t_sts[idx] #lo diezmó

B_cut = B[idx]
posicion_cut = posicion[idx]


ti_cdf = np.where(t_swia == find_nearest(t_swia, 18.1))[0][0]-15
tf_cdf = np.where(t_swia == find_nearest(t_swia, 18.2))[0][0]

ti_sts = np.where(t_sts == find_nearest(t_sts, 18.1))[0][0]
tf_sts = np.where(t_sts == find_nearest(t_sts, 18.2))[0][0]

inicio_mpb = np.where(t_sts == find_nearest(t_sts, 18.2193))[0][0]
fin_mpb = np.where(t_sts == find_nearest(t_sts, 18.2476))[0][0]
ancho_mpb = np.linalg.norm(posicion_cut[inicio_mpb] - posicion_cut[fin_mpb])

####################
density_mean = np.empty(tf_cdf-ti_cdf)
for i in range(tf_cdf-ti_cdf):
    density_mean[i] = np.mean(density[ti_cdf+i:ti_cdf+i+30]) #cgs - cm⁻³

#Elijo la densidad promedio media
idx_d = int(np.abs(np.argmin(density_mean) + np.argmax(density_mean))/2)
density_avg = density_mean[idx_d] #cgs - cm⁻³

ion_length = 2.28E07 / np.sqrt(density_avg) #cgs -  cm
print('Ion inertial length = {0:1.3g} km'.format(ion_length *1E-5))

B_avg = np.empty((len(idx), 3))
v_maven = np.empty((len(idx),3))
v_planet = np.empty((len(idx),3))
v_avg = np.empty(len(idx))

for i in range(len(idx)):
    B_avg[i,:] = np.mean(B_cut[i:i+30,:], axis=0)
    v_maven[i,:] = (posicion_cut[ti_sts+1,:] - posicion_cut[ti_sts,:]) / (t_sts[ti_sts+1]-t_sts[ti_sts]) / 3600 #en km/s
for i in range(len(idx)):
    v_planet[i,:] = np.nanmean(v_mso_xyz[i:i+30,:] + v_maven[i:i+30,:], axis=0) #km/s


v_mso = np.linalg.norm(v_mso_xyz, axis = 1)
for i in range(len(idx)):
    v_avg[i] = np.mean(v_mso[i:i+30], axis=0)

E_convective = np.cross(-v_planet * 1E3, B_avg * 1E-9)
E_convective_normalized = np.linalg.norm(E_convective, axis=1) #V/m

ti_funda = np.where(t_sts == find_nearest(t_sts, 18.06))[0][0] #los limites de la magnetofunda
tf_funda = np.where(t_sts == find_nearest(t_sts, 18.2193))[0][0]
print('El E convectivo medio en la magnetofunda es {0:1.3g} mV/m'.format(np.mean(E_convective_normalized[ti_funda:tf_funda]*1E3))) #en V/m

t1 = t_unix[np.where(t_swia == find_nearest(t_swia, 18.2193))[0][0]]
t2 = t_unix[np.where(t_swia == find_nearest(t_swia, 18.2272))[0][0]]
t3 = t_unix[np.where(t_swia == find_nearest(t_swia, 18.235))[0][0]]
t4 = t_unix[np.where(t_swia == find_nearest(t_swia, 18.2476))[0][0]]
t_bs = t_unix[np.where(t_swia == find_nearest(t_swia, 18.06))[0][0]]


n = len(E_convective_normalized)
t_inicial = t_unix[inicio]
t_final = t_unix[fin]
timestamps = np.linspace(t_inicial, t_final, n)
dates = [dt.datetime.utcfromtimestamp(ts) for ts in timestamps] #me lo da en UTC
datenums = md.date2num(dates)

# plt.figure()
# plt.subplots_adjust(bottom=0.2)
# plt.xticks( rotation=25 )
# ax = plt.gca()
xfmt = md.DateFormatter('%H:%M')
# ax.xaxis.set_major_formatter(xfmt)
# plt.plot(datenums[0:inicio_mpb], E_convective_normalized[0:inicio_mpb]*3E4)
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'black', linewidth=1)
# # plt.axhline(y = 0.0298, color = 'red', linewidth=1, label='E Hall (V/m)')
# plt.ylabel('E convectivo (V/m)')
# plt.xlabel('Tiempo')
#
#
# plt.figure()
# ax1 = plt.gca()
# ax1.xaxis.set_major_formatter(xfmt)
# # plt.plot(datenums, v_alfven*1E-5, label='Velocidad de Alfven')
# plt.plot(datenums, v_mso*1E-5, label='Velocidad de los protones') #onboard vel moment
# for xc in [md.date2num(dt.datetime.utcfromtimestamp(t1)),md.date2num(dt.datetime.utcfromtimestamp(t2)),md.date2num(dt.datetime.utcfromtimestamp(t3)),md.date2num(dt.datetime.utcfromtimestamp(t4))]:
#     plt.axvline(x = xc, color = 'black', linewidth=0.5)
# plt.xlabel('Tiempo')
# plt.ylabel('Velocidad (km/s)')
# plt.legend()

#
# fig, ax1 = plt.subplots()
#
# color = 'tab:blue'
# ax1.set_xlabel('Tiempo')
# ax1.set_ylabel('E convectivo (mV/m)', color=color)
# plt.plot(datenums[0:inicio_mpb], E_convective_normalized[0:inicio_mpb]*1E3, color = color)
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'black', linewidth=1)
# ax1.tick_params(axis='y', labelcolor=color)
# ax1.xaxis.set_major_formatter(xfmt)
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
# color = 'tab:red'
# ax2.set_ylabel('Velocidad (km/s)', color=color)  # we already handled the x-label with ax1
# plt.plot(datenums[0:inicio_mpb], v_avg[0:inicio_mpb], label='Velocidad de los protones', color=color)
# ax2.tick_params(axis='y', labelcolor=color)
#
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
#
# plt.show(block=False)


#
# fig = plt.figure(1)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
# fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.10,right=0.95, hspace = 0.005, wspace=0.15)
# plt.xticks( rotation=25 )
# xfmt = md.DateFormatter('%H:%M')
#
# ax1 = plt.gca()
# ax1.xaxis.set_major_formatter(xfmt)
# ax1 = plt.subplot2grid((4,1),(0,0))
# ax1.set_ylabel(r'$|E_{cv}|$ (mV/m)')
# plt.setp(ax1.get_xticklabels(), visible=False)
# plt.plot(datenums[0:inicio_mpb], E_convective_normalized[0:inicio_mpb]*1E3)
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'black', linewidth=1)
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t_bs)), color = 'm', linewidth=1)
# plt.grid()
#
# ax2 = plt.subplot2grid((4,1),(2,0), sharex=ax1)
# ax2.xaxis.set_major_formatter(xfmt)
# plt.setp(ax2.get_xticklabels(), visible=False)
# ax2.set_ylabel(r'$|v_{ion}^{SW}|$ (km/s)')  # we already handled the x-label with ax1
# plt.plot(datenums[0:inicio_mpb], v_avg[0:inicio_mpb], color='orange')
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'black', linewidth=1)
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t_bs)), color = 'm', linewidth=1)
# plt.setp(ax2.get_xticklabels(), visible=False)
# plt.grid()
#
# ax3 = plt.subplot2grid((4,1),(1,0), sharex=ax1)
# ax3.xaxis.set_major_formatter(xfmt)
# ax3.set_ylabel('|B| (nT)')  # we already handled the x-label with ax1
# plt.plot(datenums[0:inicio_mpb], np.linalg.norm(B_cut[0:inicio_mpb], axis=1), color='green')
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'black', linewidth=1)
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t_bs)), color = 'm', linewidth=1)
# plt.setp(ax3.get_xticklabels(), visible=False)
# plt.grid()
#
# ax4 = plt.subplot2grid((4,1),(3,0), sharex=ax1)
# ax4.xaxis.set_major_formatter(xfmt)
# ax4.set_ylabel(r'$v_{ion}^{SW}$ (km/s)')
# plt.plot(datenums[0:inicio_mpb], v_mso_xyz[0:inicio_mpb])
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'black', linewidth=1)
# plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t_bs)), color = 'm', linewidth=1)
# plt.legend([r'$v_x$ MSO', r'$v_y$ MSO', r'$v_z$ MSO'])
# ax4.set_xlabel('Tiempo')
# plt.grid()

dstart = datetime(2016,3,16,17,52)
dend = datetime(2016,3,16,18,14)

fig1 = plt.figure(2)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig1.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.10,right=0.95, hspace = 0.005, wspace=0.15)
plt.xticks( rotation=25 )
xfmt = md.DateFormatter('%H:%M')

ax1 = plt.gca()
ax1.xaxis.set_major_formatter(xfmt)
ax1 = plt.subplot2grid((2,1),(1,0))
ax1.set_ylabel(r'$|E_{cv}|$ (mV/m)')
plt.plot(datenums[0:inicio_mpb], E_convective_normalized[0:inicio_mpb]*1E3)
plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'k', linewidth=1)
ax1.axvspan(xmin = md.date2num(dt.datetime.utcfromtimestamp(t_bs)), xmax = md.date2num(dt.datetime.utcfromtimestamp(t1)), facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
ax1.legend()
plt.grid()
# ax1.set_title('MAVEN MAG SWEA LPW SWIA 16 de marzo de 2016')

ax2 = plt.subplot2grid((2,1),(0,0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
ax2.set_ylabel(r'$E_{cv}$ (mV/m)')
plt.plot(datenums[0:inicio_mpb], E_convective[0:inicio_mpb]*1E3)
plt.axvline(x = md.date2num(dt.datetime.utcfromtimestamp(t1)), color = 'k', linewidth=1)
plt.grid()
ax2.axvspan(xmin = md.date2num(dt.datetime.utcfromtimestamp(t_bs)), xmax = md.date2num(dt.datetime.utcfromtimestamp(t1)), facecolor = 'r', alpha = 0.2)
ax2.legend([r'$E_{xMSO}$',r'$E_{yMSO}$',r'$E_{zMSO}$'], loc='lower left')
ax1.set_xlabel('Tiempo')
plt.xlim(dstart, dend)
plt.setp(ax2.get_xticklabels(), visible=False)

fig1.tight_layout()

plt.show(block=False)
