import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
import scipy.signal as signal
import datetime as dt
import pdb
from scipy.stats import norm
from funciones import error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, next_available_row, datenum, unix_to_decimal, UTC_to_hdec, Bpara_Bperp
from funciones_MVA import ajuste_conico, plot_velocidades, plot_FLorentz, plot_bootstrap, bootstrap
from funciones_plot import hodograma, set_axes_equal
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import matplotlib.dates as md
import matplotlib.cm as cm
import os as os


"""
Hace el MVA y nada más. Es para correr rápido buscando la hoja de corriente.
"""


np.set_printoptions(precision=4)

date_entry = input('Enter a date in YYYY-DDD format \n')
year, doy = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

# ti_MVA = float(input("t_incial = "))
# tf_MVA = float(input("t_final = "))
ti_MVA = UTC_to_hdec(input('Tiempo inicial hh:mm:ss\n'))
tf_MVA = UTC_to_hdec(input('Tiempo final hh:mm:ss\n'))
while tf_MVA < ti_MVA:
    print('t final no puede ser menor a t inicial. \n')
    # ti_MVA = float(input('t_inicial = '))
    # tf_MVA = float(input("t_final = "))
    ti_MVA = UTC_to_hdec(input('Tiempo inicial hh:mm:ss\n'))
    tf_MVA = UTC_to_hdec(input('Tiempo final hh:mm:ss\n'))


datos_tiempo = np.loadtxt('outputs/t1t2t3t4.txt')
idx_d = np.where(int(doy) == datos_tiempo[:,1].astype(int))[0]
idx_h = np.where(int(ti_MVA) == datos_tiempo[:,2].astype(int))[0]
idx = np.intersect1d(idx_d, idx_h)[0]
t1 = datos_tiempo[idx,2]
t2 = datos_tiempo[idx,3]
t3 = datos_tiempo[idx,4]
t4 = datos_tiempo[idx,5]
tiempos = np.array([t1,t2,t3,t4])

path = '../../datos/' #path a los datos desde la laptop
ni = int(ti*32*3600 - 10000)
nf = int((24-tf)*32*3600 - 10000)
mag = np.genfromtxt(path + f'MAG_hires/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skip_header=ni, skip_footer=nf)


dia = mag[:,1]
t = mag[:,6]  #el dia decimal
t = (t - dia) * 24 #para que me de sobre la cantidad de horas

M = np.size(t) #el numero de datos

#tengo que asegurarme de que no haya agujeros en mis datos
for i in range(M-1):
    if t[i+1] - t[i] > 24 * 1.5e-5: #1.1e-5 es 1s, le doy un poco más
        print('salto en la linea {} de {} segundos'.format(i+144, (t[i+1] - t[i]) / (24*1.1e-5)))

#el campo
B = np.zeros((M, 3))
for i in range(7,10):
    B[:,i-7] = mag[:, i]

#la posición(x,y,z)
posicion = np.zeros((M, 3))
for i in range(11,14):
    posicion[:,i-11] = mag[:, i]

#la matriz diaria:
MD = np.zeros((M, 9))
MD[:, 0] = t
for i in range(1,4):
    MD[:, i] = B[:,i-1]
MD[:,4] = np.linalg.norm(B, axis=1)#la norma de B
for i in range(5,8):
    MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390 #altitud en km


inicio = np.where(t == find_nearest_inicial(t, ti_MVA))[0][0]
fin = np.where(t == find_nearest_final(t, tf_MVA))[0][0]

#################

#ahora empieza el MVA con los datos que elegí
MD_cut = MD[inicio : fin+1, :]
M_cut = np.size(MD_cut[:,0])
posicion_cut = posicion[inicio:fin+1,:]
B_cut = B[inicio:fin,:]

M_ij = Mij(B_cut)

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

av = np.concatenate([x1,x2,x3])

if x3[0] < 0: #si la normal aputna para adentro me la da vuelta
    x3 = - x3
if any(np.cross(x1,x2) - x3) > 0.01:
    print('Cambio el signo de x1 para que los av formen terna derecha')
    x1 = -x1

#las proyecciones
B1 = np.dot(B_cut, x1)
B2 = np.dot(B_cut, x2)
B3 = np.dot(B_cut, x3)


#el B medio
B_medio_vectorial = np.mean(B_cut, axis=0)
altitud = np.mean(MD_cut[:,8])
SZA = np.arccos(np.clip(np.dot(posicion_cut[0,:]/np.linalg.norm(posicion_cut[0,:]), [1,0,0]), -1.0, 1.0))* 180/np.pi

B_norm_medio = np.linalg.norm(B_medio_vectorial)

hodograma(B1, B2, B3, 'nT', 'MAVEN MAG MVA')

#el error
phi, delta_B3 = error(lamb, B_cut, M_cut, x)

###############
####fit
orbita = posicion[np.where(t == find_nearest_inicial(t, t1-1))[0][0] : np.where(t == find_nearest_final(t, t4+1))[0][0], :] / 3390 #radios marcianos

t_nave = find_nearest(t,(t2+t3)/2) #el tiempo en el medio de la hoja de corriente
index = np.where(t == t_nave)[0][0]
x0 = 0.78
e = 0.9
normal_fit, X1, Y1, Z1, R, L0 = ajuste_conico(posicion, index, orbita, x3)

B3_fit = np.dot(B_cut, normal_fit)

###############
##Bootstrap

N_boot = 1000
normal_boot, phi_boot, delta_B3_boot, out, out_phi = bootstrap(N_boot, B_cut, M_cut)

muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

B3_boot = np.dot(B_cut, normal_boot)

#######
#Errores
if phi[2,1] > phi[2,0]:
    error_normal = phi[2,1]*57.2958
else:
    error_normal = phi[2,0]*57.2958
    #quiero ver si el error más grande es phi31 o phi32

if sigma31 > sigma32:
    error_boot = sigma31
else:
    error_boot = sigma32

angulo_mva = np.arccos(np.clip(np.dot(normal_fit, x3), -1.0, 1.0))


print(f'SZA = {SZA:.3g}º y altitud = {int(altitud)}km')
print(f'MVA entre los tiempos {ti_MVA} y {tf_MVA}')
print(f'Cociente de lambdas = {lamb[1]/lamb[2]:.4g}')
print(f'El ángulo entre las normales de MVA y del fit es {angulo_mva * 180/np.pi}º')

ti = t1 - 0.15
tf = t4 + 0.15

B_para, B_perp_norm, j_inicial, j_final = Bpara_Bperp(B, t, ti, tf)
t_plot = t[j_inicial+12:j_final+12]


###############################################################################################SWEA
file_size_swea = os.path.getsize(path +  f'SWEA/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf')
if file_size_swea > 10000000:
    swea = cdf.CDF(path + f'SWEA/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf')

    flux_all = swea.varget('diff_en_fluxes')
    energia = swea.varget('energy')
    t_unix = swea.varget('time_unix')

    tu = unix_to_decimal(t_unix)
    ti_swea = np.where(tu == find_nearest(tu, ti))[0][0]
    tf_swea = np.where(tu == find_nearest(tu, tf))[0][0]
    t_swea = tu[ti_swea:tf_swea]
    flux = flux_all[ti_swea:tf_swea]
    flux_plot = np.transpose(flux)[::-1]


else:
    print('no hay datos de SWEA')

###############################################################################################SWIA
file_size_swia = os.path.getsize(path +f'SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf')
if file_size_swia > 2300000:
    swia = cdf.CDF(path + f'SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf')

    t_unix = swia.varget('time_unix')
    density = swia.varget('density')

    t_swia = unix_to_decimal(t_unix)
    inicio_swia = np.where(t_swia == find_nearest(t_swia, ti))[0][0]
    fin_swia = np.where(t_swia == find_nearest(t_swia, tf))[0][0]

    density_cut = density[inicio_swia:fin_swia]
else:
    print('no hay datos de SWIA')



index = np.array((int(year), dia[0]))


plt.clf()#clear figure
fig = plt.figure(1, constrained_layout=True)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(top = 0.93, bottom = 0.07, left = 0.05,right=0.95, hspace = 0.005, wspace=0.15)
fig.set_size_inches(15, 10)#con este tamaño ocupa toda la pantalla de la laptop

ax1 = plt.subplot2grid((3,2),(0,0))
plt.plot(t_plot, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / B')
plt.plot(t_plot, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / B')
plt.setp(ax1.get_xticklabels(), visible=False)
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
plt.axvline(x = ti_MVA, color = 'r', linewidth=1)
plt.axvline(x = tf_MVA, color = 'r', linewidth=1)
ax1.set_ylabel(r'|$\Delta B$|/ B')
ax1.grid()
ax1.legend()

ax4 = plt.subplot2grid((3,2),(1,0), sharex=ax1)
ax4.plot(t_plot, B[j_inicial:j_final,1], label='By')
ax4.plot(t_plot, B[j_inicial:j_final,0], label='Bx')
ax4.plot(t_plot, B[j_inicial:j_final,2], label='Bz')
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
plt.axvline(x = ti_MVA, color = 'r', linewidth=1)
plt.axvline(x = tf_MVA, color = 'r', linewidth=1)
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.set_ylabel('Bx, By, Bz (nT)')
ax4.legend()
ax4.grid()

ax3 = plt.subplot2grid((3,2),(2,0), sharex=ax1)
plt.plot(t_plot, MD[j_inicial:j_final,4])
ax3.grid()
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
plt.axvline(x = ti_MVA, color = 'r', linewidth=1)
plt.axvline(x = tf_MVA, color = 'r', linewidth=1)
ax3.set_ylabel('|B| (nT)')
ax3.set_xlabel('Tiempo (hdec)')

ax5 = plt.subplot2grid((3,2),(0,1), sharex=ax1)
ax5.set_ylabel('Energia', picker=True)#, bbox=dict(facecolor='red'))
plt.setp(ax5.get_xticklabels(), visible=False)
for xc in tiempos:
    plt.axvline(x = xc, color = 'k', linewidth=1)
plt.axvline(x = ti_MVA, color = 'r', linewidth=1)
plt.axvline(x = tf_MVA, color = 'r', linewidth=1)
if file_size_swea > 10000000:
    im = plt.imshow(flux_plot, aspect = 'auto',origin = 'lower', extent=(t_swea[0], t_swea[-1],  energia[-1], energia[0]), cmap='inferno', norm=LogNorm(vmin=1E4, vmax=1E9))
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes("top", size="7%", pad="1%")
    cb = plt.colorbar(im, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")



ax7 = plt.subplot2grid((3,2),(1,1), sharex=ax1)
plt.setp(ax7.get_xticklabels(), visible=False)
ax7.set_ylabel('Densidad de p+ \n del SW (cm⁻³)')
if file_size_swia > 2300000:
    plt.plot(t_swia[inicio_swia:fin_swia], density_cut)
    plt.axvline(x = ti_MVA, color = 'r', linewidth=1)
    plt.axvline(x = tf_MVA, color = 'r', linewidth=1)
    ax7.grid()
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)

plt.suptitle(f'MAVEN {year}-{doy}')




plt.show(block=False)
#############
