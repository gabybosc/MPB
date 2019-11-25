import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cdflib as cdf
import datetime as dt
from funciones import find_nearest, deltaB, unix_to_decimal, unix_to_timestamp, Bpara_Bperp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import matplotlib.dates as md
import matplotlib.cm as cm
import os as os
import pickle


"""
Con los tiempos t1 t2 t3 t4 ploteo los datos de MAG, SWEA, SWIA y LPW. Les dibujo encima líneas que correspondan a los tiempos.
Guarda las figuras como png y como pickle (para poder reabrirlas con python y que sean interactivas)
"""
#se fija que los datos no estén vacíos antes de cargarlos (para ahorrarme el error)


np.set_printoptions(precision=4)

fechas = np.loadtxt('../outputs/t1t2t3t4.txt')
# for i in range(len(fechas)):
i = 16
for j in range(1):
    year = int(fechas[i,0])
    doy = int(fechas[i, 1])
    t1 = fechas[i,2]
    t2 = fechas[i,3]
    t3 = fechas[i,4]
    t4 = fechas[i,5]
    tiempos = np.array([t1,t2,t3,t4])

    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date
    year = date_orbit.strftime("%Y")
    doy = date_orbit.strftime("%j")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    # path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
    path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
    if 1==2:#os.path.isfile(path + 'mag_filtrado.txt'):
        mag = np.loadtxt(path + 'mag_filtrado.txt', skiprows=2)
        M = len(mag[:,0]) #el numero de datos
        B = mag[:, :3]

        Bnorm = mag[:,-1]
        mag = np.loadtxt(path + 'MAG.asc')
        Bxyz_paraperp = mag[:,6:9]
    else:
        mag = np.loadtxt(path + 'MAG.asc')
        M = len(mag[:,0]) #el numero de datos
        B = mag[:, 6:9]
        Bnorm = np.linalg.norm(B, axis=1)
        Bxyz_paraperp = mag[:,6:9]

    hh = mag[:,3]
    mm = mag[:,4]
    ss = mag[:,5]

    t = hh + mm/60 + ss/3600 #hdec

    M = np.size(t) #el numero de datos

    #el campo
    B = np.zeros((M, 3))
    for i in range(6,9):
        B[:,i-6] = mag[:, i]

    Bnorm = mag[:,-1]

    #la posición(x,y,z)
    posicion = np.zeros((M, 3))
    for i in range(9,12):
        posicion[:,i-9] = mag[:, i]

    #la matriz diaria:
    MD = np.zeros((M, 9))
    MD[:, 0] = t
    for i in range(1,4):
        MD[:, i] = B[:,i-1]
    MD[:,4] = Bnorm
    for i in range(5,8):
        MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
    MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390 #altitud en km

    ti = t1 - 0.15
    tf = t4 + 0.15

    B_para, B_perp_norm, j_inicial, j_final = Bpara_Bperp(Bxyz_paraperp, t, ti, tf)
    t_plot = t[j_inicial+12:j_final+12]

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]
    ###############################################################################################SWEA

    swea = np.loadtxt(path + 'SWEA.asc')

    energy = swea[:, 7]
    JE_total = swea[:, -1]

    t_swea = np.unique(swea[:,3] + swea[:,4]/60 + swea[:,5]/3600) #hdec
    # ti = 10.5
    # tf = 11
    inicio_swea = np.where(t_swea == find_nearest(t_swea, ti))[0][0]
    fin_swea = np.where(t_swea == find_nearest(t_swea, tf))[0][0]

    energias = [100 + i*20 for i in range(6)]


    ###############################################################################################SWIA
    swia = np.loadtxt(path + 'SWIA.asc')

    density = swia[:,-1]

    t_swia = swia[:,3] + swia[:,4]/60 + swia[:,5]/3600 #hdec
    inicio_swia = np.where(t_swia == find_nearest(t_swia, ti))[0][0]
    fin_swia = np.where(t_swia == find_nearest(t_swia, tf))[0][0]


    ############################################################################################### LPW
    lpw = np.loadtxt(path + 'LPW.asc')

    e_density = lpw[:,-1]

    t_lpw = lpw[:,3] + lpw[:,4]/60 + lpw[:,5]/3600
    inicio_lpw = np.where(t_lpw == find_nearest(t_lpw, ti))[0][0]
    fin_lpw = np.where(t_lpw == find_nearest(t_lpw, tf))[0][0]





    index = np.array((int(year), int(day)))


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
    ax1.set_ylabel(r'|$\Delta B$|/ B')
    ax1.grid()
    ax1.legend()

    ax4 = plt.subplot2grid((3,2),(1,0), sharex=ax1)
    ax4.plot(t[inicio:fin], B[inicio:fin,0], label='Bx')
    ax4.plot(t[inicio:fin], B[inicio:fin,1], label='By')
    ax4.plot(t[inicio:fin], B[inicio:fin,2], label='Bz')
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.set_ylabel('Bx, By, Bz (nT)')
    ax4.legend()
    ax4.grid()

    ax3 = plt.subplot2grid((3,2),(2,0), sharex=ax1)
    plt.plot(t[inicio:fin], Bnorm[inicio:fin])
    ax3.grid()
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)
    ax3.set_ylabel('|B| (nT)')
    ax3.set_xlabel('Tiempo (hdec)')

    ax5 = plt.subplot2grid((3,2),(0,1), sharex=ax1)
    for energia in energias:
        index = np.where(energy == find_nearest(energy, energia))[0]
        JE = JE_total[index]
        plt.semilogy(t_swea[inicio_swea:fin_swea], JE[inicio_swea:fin_swea], label = f'{energia} eV')
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)
    ax5.set_xlabel('tiempo (hdec)')
    ax5.set_ylabel('diff en flux')
    ax5.legend()



    ax7 = plt.subplot2grid((3,2),(1,1), sharex=ax1)
    plt.setp(ax7.get_xticklabels(), visible=False)
    ax7.set_ylabel('Densidad de p+ \n del SW (cm⁻³)')
    plt.plot(t_swia[inicio_swia:fin_swia], density[inicio_swia:fin_swia])
    ax7.grid()
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)


    ax6 = plt.subplot2grid((3,2),(2,1), sharex=ax1)
    ax6.set_ylabel('Densidad total \n de e- (cm⁻³)')
    ax6.set_xlabel('Tiempo (hdec)')
    plt.semilogy(t_lpw[inicio_lpw:fin_lpw], e_density[inicio_lpw:fin_lpw])
    ax6.grid()
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)

    plt.suptitle(f'MAVEN {year}-{month}-{day}')




    plt.show(block=False)
    plt.savefig(f'../outputs/figs_MPB/MPB_y{year}_d{doy}_t{int(t1)}.png', dpi=200)
    pickle.dump(fig, open(f'../outputs/figs_MPB/MPB_y{year}_d{doy}_t{int(t1)}.pkl', 'wb'))
