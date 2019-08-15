import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cdflib as cdf
import datetime as dt
from funciones import find_nearest, deltaB, unix_to_decimal, unix_to_timestamp
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

fechas = np.loadtxt('outputs/t1t2t3t4.txt')
for i in range(len(fechas)):
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
    path = '../../datos/' #path a los datos desde la laptop
    mag = np.loadtxt(path + f'MAG_1s/2016/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=148) #datos MAG 1s (para plotear no quiero los datos pesados)
    n =2
    mag = mag[:-n, :] #borra las ultimas 2 filas, que es ya el dia siguiente (no sé si siempre)

    dia = mag[:,1]
    t = mag[:,6]  #el dia decimal
    t = (t - dia) * 24 #hdec

    M = np.size(t) #el numero de datos

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
    MD[:,4] = np.linalg.norm(B, axis = 1)
    for i in range(5,8):
        MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
    MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390 #altitud en km

    ti = t1 - 0.15
    tf = t4 + 0.15

    j_inicial = np.where(t == find_nearest(t, ti))[0][0]
    j_final =  np.where(t == find_nearest(t, tf))[0][0]

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

    ############################################################################################### LPW
    file_size_lpw = os.path.getsize(path +  f'LPW/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf')
    if file_size_lpw > 3000000:
        lpw = cdf.CDF(path + f'LPW/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf')

        t_unix = lpw.varget('time_unix')
        e_density = lpw.varget('data')[:,3]

        t_lpw = unix_to_decimal(t_unix)
        inicio_lpw = np.where(t_lpw == find_nearest(t_lpw, ti))[0][0]
        fin_lpw = np.where(t_lpw == find_nearest(t_lpw, tf))[0][0]

    else:
        print('no hay datos de LPW')



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
    ax1.set_ylabel(r'|$\Delta B$|/ B')
    ax1.grid()
    ax1.legend()

    ax4 = plt.subplot2grid((3,2),(1,0), sharex=ax1)
    ax4.plot(t[j_inicial+12:j_final+12], B[j_inicial:j_final,1], label='By')
    ax4.plot(t[j_inicial+12:j_final+12], B[j_inicial:j_final,0], label='Bx')
    ax4.plot(t[j_inicial + 12: j_final +12], B[j_inicial:j_final,2], label='Bz')
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.set_ylabel('Bx, By, Bz (nT)')
    ax4.legend()
    ax4.grid()

    ax3 = plt.subplot2grid((3,2),(2,0), sharex=ax1)
    plt.plot(t[j_inicial + 12: j_final +12], MD[j_inicial:j_final,4])
    ax3.grid()
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)
    ax3.set_ylabel('|B| (nT)')
    ax3.set_xlabel('Tiempo (hdec)')

    ax5 = plt.subplot2grid((3,2),(0,1), sharex=ax1)
    ax5.set_ylabel('Energia', picker=True)#, bbox=dict(facecolor='red'))
    plt.setp(ax5.get_xticklabels(), visible=False)
    for xc in tiempos:
        plt.axvline(x = xc, color = 'k', linewidth=1)
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
        ax7.grid()
        for xc in tiempos:
            plt.axvline(x = xc, color = 'k', linewidth=1)


    ax6 = plt.subplot2grid((3,2),(2,1), sharex=ax1)
    ax6.set_ylabel('Densidad total \n de e- (cm⁻³)')
    ax6.set_xlabel('Tiempo (hdec)')
    if file_size_lpw > 3000000:
        plt.semilogy(t_lpw[inicio_lpw:fin_lpw], e_density[inicio_lpw:fin_lpw])
        ax6.grid()
        for xc in tiempos:
            plt.axvline(x = xc, color = 'k', linewidth=1)

    plt.suptitle(f'MAVEN {year}-{doy}')




    plt.show(block=False)
    plt.savefig(f'outputs/MPB_y{year}_d{doy}_t{int(t1)}.png', dpi=200)
    pickle.dump(fig, open(f'outputs/MPB_y{year}_d{doy}_t{int(t1)}.pkl', 'wb'))
