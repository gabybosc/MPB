from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, MultiCursor
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
from datetime import datetime
from funciones import find_nearest, deltaB, onpick1
from funcion_flujo_energia_cdf import flujo_energia
import time
import matplotlib.cm as cm

np.set_printoptions(precision=4)

#DATOS DE PDS
dia = 16#input("dia del mes = ")
diaa = 76#input('dia del año = ')
path = '../datos/marzo 2016/{}/'.format(dia) #path a los datos
datos = np.loadtxt(path + 'mvn_mag_l2_20160{}ss1s_201603{}_v01_r01.sts'.format(diaa, dia), skiprows=148) #lee todo y me da todo
n =2
datos = datos[:-n, :] #borra las ultimas 2 filas, que es ya el dia siguiente (no sé si siempre)
cdf_swea = cdf.CDF(path + 'mvn_swe_l2_svyspec_201603{}_v04_r01.cdf'.format(dia))


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
MD[:,4] = np.sqrt(B[:,0]**2 + B[:,1]**2 + B[:,2]**2)#el modulo de B
for i in range(5,8):
    MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
MD[:, 8] = np.sqrt(posicion[:,0]**2 + posicion[:,1]**2 + posicion[:,2]**2) - 3390 #altitud en km

#Si quiero elegir manualmente la orbita:
# plt.plot(t, MD[:,4])
# plt.xlabel('t (hdec)')
# plt.ylabel('|B|')
# plt.title('Orbitas')
# plt.show()

ti = 18.16#float(input("Tiempo inicial = "))
tf = 18.28#float(input("Tiempo final = "))
orbit_number = 5#float(input('Número de órbita = '))
# if tf < ti:
#     print('El tiempo inicial no puede ser mayor al final')
#     ti = float(input("Tiempo inicial = "))
#     tf = float(input("Tiempo final = "))

t1 = find_nearest(t, ti)
t2 = find_nearest(t, tf)
j_inicial = np.where(t == t1)[0][0]
j_final =  np.where(t == t2)[0][0]

#el deltaB es una norma. Corresponde al t del medio del intervalo.
#Lo hago en ventanas de 60s, moviendose de a 1s.
B_para = np.zeros(j_final - j_inicial)
B_perp = np.zeros((j_final - j_inicial, 3))
B_perp_norm = np.zeros(j_final - j_inicial)
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

theta = np.arccos(B[:,2]/MD[:,4]) * 57.2958#cos(theta) = Bz/|B|
phi = np.arctan2(B[:,0], B[:,1])* 57.2958#tg(phi) = By/Bx
theta_cut = theta[j_inicial+12:j_final+12]
phi_cut = phi[j_inicial+12:j_final+12]

t_plot = t[j_inicial+12:j_final+12]
t_flux, flux_cut, E_flux = flujo_energia(t1, t2, cdf_swea) #t_flux es diferente, pues viene del cdf

index = np.array((year, dia[0], orbit_number))

happy = False
while not happy:
    val = []
    while len(val) < 4:
        plt.clf()#clear figure
        fig = plt.figure(1)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
        fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.05,right=0.95, hspace = 0.005, wspace=0.15)

        ax1 = plt.subplot2grid((3,3),(0,0))#421
        plt.plot(t_plot, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / B')
        plt.plot(t_plot, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / B')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r'|$\Delta B$|/ B')
        ax1.grid()
        ax1.legend()
        ax2 = plt.subplot2grid((3,3),(1,0), sharex=ax1,sharey = ax1)
        plt.plot(t_plot, B_perp[:,0], label='x')
        plt.plot(t_plot, B_perp[:,1], label='y')
        plt.plot(t_plot, B_perp[:,2], label='z')
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.set_ylabel(r'$\Delta B \perp$/ B')
        ax2.grid()
        ax2.legend()
        ax3 = plt.subplot2grid((3,3),(0,1), sharex=ax1)
        plt.plot(t[j_inicial + 12: j_final +12], MD[j_inicial:j_final,4])
        ax3.grid()
        ax3.set_ylabel('|B| (nT)')
        plt.setp(ax3.get_xticklabels(), visible=False)

        ax4 = plt.subplot2grid((3,3),(2,0), sharex=ax1)
        ax4.plot(t[j_inicial+12:j_final+12], B[j_inicial:j_final,1], label='By')
        ax4.plot(t[j_inicial+12:j_final+12], B[j_inicial:j_final,0], label='Bx')
        ax4.plot(t[j_inicial + 12: j_final +12], B[j_inicial:j_final,2], label='Bz')
        ax4.set_ylabel('Bx, By, Bz (nT)')
        ax4.legend()
        ax4.set_xlabel('Tiempo (hdec)')
        ax4.grid()

        ax6 = plt.subplot2grid((3,3),(0,2), sharex=ax1)
        ax6.set_ylabel(r'$\Theta (º)$')#, bbox=dict(facecolor='red'))
        ax6.plot(t_plot, theta_cut, linewidth=0.5)
        plt.setp(ax6.get_xticklabels(), visible=False)
        ax6.grid()
        ax7 = plt.subplot2grid((3,3),(1,2), sharex=ax1)
        plt.plot(t_plot, phi_cut, linewidth=0.5)
        # plt.setp(ax7.get_xticklabels(), visible=False)
        ax7.set_ylabel(r'$\phi (º)$')
        ax7.grid()
        ax7.set_xlabel('Tiempo (h)')

        ax5 = plt.subplot2grid((3,3),(1,1), rowspan=2, sharex=ax1)
        ax5.set_ylabel('Flujo (cm⁻² sr⁻¹ s⁻¹)', picker=True)#, bbox=dict(facecolor='red'))
        line, = ax5.semilogy(t_flux, flux_cut[:,0], linewidth=1, label = '{0:1.4g} eV'.format(E_flux[0]), picker=5)
        fig.canvas.mpl_connect('pick_event', onpick1)
        for j in range(1, len(flux_cut[0,:])):
            ax5.semilogy(t_flux, flux_cut[:,j], linewidth=1, label = '{0:1.4g} eV'.format(E_flux[j]))
        ax5.legend(loc='lower right', bbox_to_anchor=(1.5, 0))
        ax5.set_xlabel('Tiempo (hdec)')
        ax5.grid()

        multi = MultiCursor(fig.canvas, (ax1, ax2,ax3,ax4,ax5,ax6,ax7), color='black', lw=0.5)

        zoom_ok = False
        print('\nSpacebar when ready to click:\n')
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress()
        print('Click to select MPB: ')
        val = np.asarray(plt.ginput(4))[:,0]
        print('Selected values: ', val)
        outs = np.concatenate((index, val))
        outs = sorted(outs[3:7])

    print('Happy? Keyboard click for yes, mouse click for no.')
    happy = plt.waitforbuttonpress()
with open('tiempos.txt','ab') as f:
     np.savetxt(f, [outs], fmt='%5g', delimiter=',')


plt.show(block=False)
