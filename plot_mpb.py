import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, MultiCursor
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
from datetime import datetime
from funciones import find_nearest, deltaB, onpick1, unix_to_decimal, plot_select, set_axes_equal, datenum
from funcion_flujo_energia_cdf import flujo_energia
import matplotlib.dates as md
import datetime as dt
import time
from datetime import datetime
import matplotlib.cm as cm

np.set_printoptions(precision=4)

#DATOS DE PDS
# dia = input("dia del mes = ")
# diaa = input('dia del año = ')
dia = 16
diaa = 76
path = '../datos/marzo 2016/{}/'.format(dia) #path a los datos
datos = np.loadtxt(path + 'mvn_mag_l2_20160{}ss1s_201603{}_v01_r01.sts'.format(diaa, dia), skiprows=148) #lee todo y me da todo
n =2
datos = datos[:-n, :] #borra las ultimas 2 filas, que es ya el dia siguiente (no sé si siempre)
cdf_swea = cdf.CDF(path + 'mvn_swe_l2_svyspec_201603{}_v04_r01.cdf'.format(dia))
cdf_swia = cdf.CDF(path + 'mvn_swi_l2_onboardsvymom_201603{}_v01_r01.cdf'.format(dia))
lpw = np.loadtxt(path + 'mvn_kp_insitu_20160316_v14_r03.tab') #son los datos entre las 18 y las 19h

t_lpw = lpw[:,0] + lpw[:,1]/60 + lpw[:,2]/3600

dia = datos[:,1]
t = datos[:,6]  #el dia decimal
t = (t - dia) * 24 #hdec
year = datos[0,0]

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
MD[:, 8] = np.sqrt(posicion[:,0]**2 + posicion[:,1]**2 + posicion[:,2]**2) - 3390 #altitud en km


t1 = find_nearest(t, 17.85)
t2 = find_nearest(t, 18.4)
j_inicial = np.where(t == t1)[0][0]
j_final =  np.where(t == t2)[0][0]

#el deltaB es una norma. Corresponde al t del medio del intervalo.
#Lo hago en ventanas de 60s, moviendose de a 1s.
B_para = np.empty(j_final - j_inicial)
B_perp = np.empty((j_final - j_inicial, 3))
B_perp_norm = np.empty(j_final - j_inicial)
for j in range(j_inicial, j_final):
    Mi = j
    Mf = j + 25
    M_delta = 25
    B_delta = B[Mi:Mf]
    t_delta = t[Mi:Mf]
    #hasta aca importa los datos
    #ahora quiero que haga el delta para cada j
    deltaB_para, deltaB_perp = deltaB(B_delta)
    # pero quiero que lo guarde en otro B. Va a corresponder al t de la mtiad del intervalo
    B_para[j-j_inicial] = deltaB_para[12]
    B_perp[j-j_inicial, :] = deltaB_perp[12, :]
    B_perp_norm[j-j_inicial] = np.linalg.norm(deltaB_perp[12,:])


########## SWEA
t_flux, flux_cut, E_flux = flujo_energia(17.85, 18.4, cdf_swea) #t_flux es diferente, pues viene del cdf

t_unix_e = cdf_swea.varget('time_unix')

t_swea = unix_to_decimal(t_unix_e)
inicio_e = np.where(t_swea == find_nearest(t_swea, 17.85))[0][0]
fin_e = np.where(t_swea == find_nearest(t_swea, 18.4))[0][0]

ne = len(flux_cut)
ti_swea = t_unix_e[inicio_e]
tf_swea = t_unix_e[fin_e]
timestamps_swea = np.linspace(ti_swea, tf_swea, ne)
dates_swea = [dt.datetime.utcfromtimestamp(ts) for ts in timestamps_swea] #me lo da en UTC
datenums_swea = md.date2num(dates_swea)


######## densidad SWIA
t_unix = cdf_swia.varget('time_unix')
density = cdf_swia.varget('density')

t_swia = unix_to_decimal(t_unix)
inicio = np.where(t_swia == find_nearest(t_swia, 17.85))[0][0]
fin = np.where(t_swia == find_nearest(t_swia, 18.4))[0][0]

density_cut = density[inicio:fin]

n = len(density_cut)
ti_swia = t_unix[inicio]
tf_swia = t_unix[fin]
timestamps_swia = np.linspace(ti_swia, tf_swia, n)
dates_swia = [dt.datetime.utcfromtimestamp(ts) for ts in timestamps_swia] #me lo da en UTC
datenums_swia = md.date2num(dates_swia)

####### densidad electrones
e_density = lpw[:,3]
ti_lpw = np.where(t_lpw == find_nearest(t_lpw, 17.85))[0][0]
tf_lpw = np.where(t_lpw == find_nearest(t_lpw, 18.4))[0][0]

tiempo_lpw = np.array([np.datetime64(datetime(2016,3,16,int(x[0]), int(x[1]), int(x[2]))) for x in lpw[ti_lpw:tf_lpw,:]]) #un array de datetimes con las horas entre 17.85 y 18.4

############# MAG
tiempo_mag = np.array([np.datetime64(datenum(2016, 3, 16, x)) for x in t[j_inicial:j_final]]) #datenum es una función mía

###########el plot

t1 = datetime(2016,3,16,18,13,00)
t2 = np.where(t[j_inicial:j_final] == find_nearest(t, 18.2201))[0][0]
t3 = np.where(t[j_inicial:j_final] == find_nearest(t, 18.235))[0][0]
t4 = np.where(t[j_inicial:j_final] == find_nearest(t, 18.2476))[0][0]
t_bs = np.where(t[j_inicial:j_final] == find_nearest(t, 18.05))[0][0]

dstart = datetime(2016,3,16,17,50)
dend = datetime(2016,3,16,18,25)

plt.clf()#clear figure
fig = plt.figure(1)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.10,right=0.95, hspace = 0.005, wspace=0.15)
plt.xticks( rotation=25 )
xfmt = md.DateFormatter('%H:%M:%S')

ax1 = plt.gca()
ax1.xaxis.set_major_formatter(xfmt)
ax1 = plt.subplot2grid((5,1),(0,0))
ax1.plot(tiempo_mag, B[j_inicial:j_final,0], label='Bx MSO')
ax1.plot(tiempo_mag, B[j_inicial:j_final,1], label='By MSO')
ax1.plot(tiempo_mag, B[j_inicial:j_final,2], label='Bz MSO')
# for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
#     plt.axvline(x = xc, color = 'k', linewidth=1)
ax1.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[t_bs], facecolor = 'b', alpha = 0.2)
ax1.axvspan(xmin = tiempo_mag[t_bs], xmax = t1, facecolor = 'r', alpha = 0.2)
ax1.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,19,13), facecolor = 'yellow', alpha = 0.2)
# plt.axvline(x = tiempo_mag[t_bs], color = 'm', linewidth=1)
ax1.set_ylabel('Componentes de B (nT)')
ax1.legend(loc='lower left',fontsize='small')
ax1.grid()
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_title('MAVEN MAG SWEA LPW SWIA 16 de marzo de 2016')

ax2 = plt.subplot2grid((5,1),(1,0), sharex=ax1)
ax2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, MD[j_inicial:j_final,4])
# for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
#     plt.axvline(x = xc, color = 'k', linewidth=1)
# plt.axvline(x = tiempo_mag[t_bs], color = 'm', linewidth=1)
ax2.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[t_bs], facecolor = 'b', alpha = 0.2)
ax2.axvspan(xmin = tiempo_mag[t_bs], xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
ax2.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,19,13), facecolor = 'yellow', alpha = 0.2, label = 'MPR')
plt.grid()
plt.ylabel('|B| (nT)')
plt.setp(ax2.get_xticklabels(), visible=False)

ax3 = plt.subplot2grid((5,1),(2,0), sharex=ax1)
ax3.xaxis.set_major_formatter(xfmt)
ax3.set_ylabel('Flujo de e- del SW \n (cm⁻² sr⁻¹ s⁻¹)')
line, = ax3.semilogy(datenums_swea, flux_cut[:,0], linewidth=1, label = '{} eV'.format(int(E_flux[0])), picker=5)
for j in range(1, len(flux_cut[0,:])):
    ax3.semilogy(datenums_swea, flux_cut[:,j], linewidth=1, label = '{} eV'.format(int(E_flux[j])))
    ax3.legend(loc='lower right', bbox_to_anchor=(1.5, 0))
# for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
#     plt.axvline(x = xc, color = 'k', linewidth=1)
# plt.axvline(x = tiempo_mag[t_bs], color = 'm', linewidth=1)
ax3.axvspan(xmin = tiempo_mag[t_bs], xmax = t1, facecolor = 'r', alpha = 0.2)
ax3.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,19,13), facecolor = 'yellow', alpha = 0.2)
# ax3.legend(fontsize='small', loc='left')
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.grid()

ax4 = plt.subplot2grid((5,1),(3,0), sharex=ax1)
plt.semilogy(tiempo_lpw, e_density[ti_lpw:tf_lpw])
# for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
#     plt.axvline(x = xc, color = 'k', linewidth=1)
# plt.axvline(x = tiempo_mag[t_bs], color = 'm', linewidth=1)
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.axvspan(xmin = tiempo_mag[t_bs], xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
ax4.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,19,13), facecolor = 'yellow', alpha = 0.2, label = 'MPR')
ax4.grid()
ax4.set_ylabel('Densidad total \n de e- (cm⁻³)')

ax5 = plt.subplot2grid((5,1),(4,0), sharex=ax1)
ax5.xaxis.set_major_formatter(xfmt)
plt.plot(datenums_swia, density_cut)
# for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
#     plt.axvline(x = xc, color = 'k', linewidth=1)
# plt.axvline(x = tiempo_mag[t_bs], color = 'm', linewidth=1)
ax5.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[t_bs], facecolor = 'b', alpha = 0.2, label = 'Viento solar no perturbado')
ax5.axvspan(xmin = tiempo_mag[t_bs], xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
ax5.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,19,13), facecolor = 'yellow', alpha = 0.2, label = 'MPR')
ax5.set_ylabel('Densidad de p+ \n del SW (cm⁻³)')
ax5.set_xlabel('Tiempo (UTC)')
ax5.grid()
ax5.legend(loc='upper left',fontsize='small')
ax5.set_xlim(dstart, dend)

ax3.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[t_bs], facecolor = 'b', alpha = 0.2)
ax4.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[t_bs], facecolor = 'b', alpha = 0.2)
plt.tight_layout()

###### B para y B perp
fig2 = plt.figure(2)
fig2.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.10,right=0.95, hspace = 0.005, wspace=0.15)
ax = plt.gca()
ax.xaxis.set_major_formatter(xfmt)
plt.semilogy(tiempo_mag, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / |B|')
plt.semilogy(tiempo_mag, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / |B|')
for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
ax.axvspan(xmin = datetime(2016,3,16,18,5), xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
ax.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,19,13), facecolor = 'yellow', alpha = 0.2, label = 'MPR')
plt.grid()
plt.legend()
plt.xlabel('Tiempo (UTC)')
plt.ylabel('Variación de relativa de B')
plt.title('MAVEN MAG 16 de marzo de 2016')
ax.set_xlim(datetime(2016,3,16,18,5),datetime(2016,3,16,18,20))

##########3 mpb zoom

for i in range(3,8):

    fig = plt.figure(i, figsize=(8,30))#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
    fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.12,right=0.95, hspace = 0.0, wspace=0.15)
    plt.xticks( rotation=25 )
    xfmt = md.DateFormatter('%H:%M')

    axz1 = plt.gca()
    axz1.xaxis.set_major_formatter(xfmt)
    axz1 = plt.subplot2grid((6,1),(0,0))
    axz1.plot(tiempo_mag, B[j_inicial:j_final,0], label='Bx MSO')
    axz1.plot(tiempo_mag, B[j_inicial:j_final,1], label='By MSO')
    axz1.plot(tiempo_mag, B[j_inicial:j_final,2], label='Bz MSO')
    for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
        plt.axvline(x = xc, color = 'k', linewidth=1.5)
    axz1.set_ylabel('Componentes de B (nT)')
    axz1.legend(loc='lower right')
    axz1.grid()
    axz1.set_xlim(datetime(2016,3,16,18,10),datetime(2016,3,16,18,20))
    axz1.set_title('MAVEN MAG SWEA LPW SWIA 16 de marzo de 2016')
    plt.setp(axz1.get_xticklabels(), visible=False)

    axz2 = plt.subplot2grid((6,1),(1,0), sharex=axz1)
    axz2.xaxis.set_major_formatter(xfmt)
    plt.plot(tiempo_mag, MD[j_inicial:j_final,4])
    for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
        plt.axvline(x = xc, color = 'k', linewidth=1.5)
    plt.grid()
    plt.ylabel('|B| (nT)')
    axz2.set_xlim(datetime(2016,3,16,18,10),datetime(2016,3,16,18,20))
    axz1.axvspan(xmin = datetime(2016,3,16,18,10), xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
    axz1.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,20), facecolor = 'yellow', alpha = 0.2, label = 'MPR')
    axz2.axvspan(xmin = datetime(2016,3,16,18,10), xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
    axz2.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,20), facecolor = 'yellow', alpha = 0.2, label = 'MPR')

    if i >3:

        axz3 = plt.subplot2grid((6,1),(2,0), sharex=axz1)
        axz3.xaxis.set_major_formatter(xfmt)
        axz3.set_ylabel('Flujo dif de energia \n de e- del SW \n (cm⁻² sr⁻¹ s⁻¹)')
        line, = axz3.semilogy(datenums_swea, flux_cut[:,0], linewidth=1, label = '{} eV'.format(int(E_flux[0])), picker=5)
        for j in range(1, len(flux_cut[0,:])):
            axz3.semilogy(datenums_swea, flux_cut[:,j], linewidth=1, label = '{} eV'.format(int(E_flux[j])))
            axz3.legend(loc='lower right', bbox_to_anchor=(1.5, 0))
        for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
            plt.axvline(x = xc, color = 'k', linewidth=1.5)
        axz3.legend(fontsize='small', loc='right')
        axz3.grid()
        axz3.set_xlim(datetime(2016,3,16,18,10),datetime(2016,3,16,18,20))

        plt.setp(axz2.get_xticklabels(), visible=False)
        axz3.axvspan(xmin = datetime(2016,3,16,18,10), xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
        axz3.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,20), facecolor = 'yellow', alpha = 0.2, label = 'MPR')

    if i>4:
        axz4 = plt.subplot2grid((6,1),(3,0), sharex=axz1)
        plt.semilogy(tiempo_lpw, e_density[ti_lpw:tf_lpw])
        for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
            plt.axvline(x = xc, color = 'k', linewidth=1.5)
        axz4.grid()
        axz4.set_ylabel('Densidad total \n de e- (cm⁻³)')
        axz4.set_xlim(datetime(2016,3,16,18,10),datetime(2016,3,16,18,20))
        plt.setp(axz3.get_xticklabels(), visible=False)
        axz4.axvspan(xmin = datetime(2016,3,16,18,10), xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
        axz4.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,20), facecolor = 'yellow', alpha = 0.2, label = 'MPR')

    if i >5:
        axz5 = plt.subplot2grid((6,1),(4,0), sharex=axz1)
        axz5.xaxis.set_major_formatter(xfmt)
        plt.plot(datenums_swia, density_cut)
        for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
            plt.axvline(x = xc, color = 'k', linewidth=1.5)
        axz5.set_ylabel('Densidad de p+ \n del SW (cm⁻³)')
        axz5.grid()
        axz5.set_xlim(datetime(2016,3,16,18,10),datetime(2016,3,16,18,20))
        plt.setp(axz4.get_xticklabels(), visible=False)
        axz5.axvspan(xmin = datetime(2016,3,16,18,10), xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
        axz5.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,20), facecolor = 'yellow', alpha = 0.2, label = 'MPR')

    if i>6:
        axz6 = plt.subplot2grid((6,1),(5,0), sharex=axz1)
        axz6.semilogy(tiempo_mag, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / |B|')
        axz6.semilogy(tiempo_mag, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / |B|')
        for xc in [t1,tiempo_mag[t2],tiempo_mag[t3],tiempo_mag[t4]]:
            plt.axvline(x = xc, color = 'k', linewidth=1.5)
        axz6.grid()
        axz6.legend()
        axz6.set_xlabel('Tiempo (UTC)')
        axz6.set_ylabel('Variación relativa \n de B')
        axz6.set_xlim(datetime(2016,3,16,18,10),datetime(2016,3,16,18,20))

        plt.setp(axz5.get_xticklabels(), visible=False)
        axz6.axvspan(xmin = datetime(2016,3,16,18,10), xmax = t1, facecolor = 'r', alpha = 0.2, label = 'Magnetofunda')
        axz6.axvspan(xmin = tiempo_mag[t4], xmax = datetime(2016,3,16,18,20), facecolor = 'yellow', alpha = 0.2, label = 'MPR')



    # plt.tight_layout()
plt.show(block=False)
