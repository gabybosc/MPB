import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
from funciones import find_nearest, unix_to_decimal, plot_select

# np.set_printoptions(precision=4)
def flujo_energia(t1, t2, cdf_file):
    flux = cdf_file.varget('diff_en_fluxes')
    energia = cdf_file.varget('energy')
    t_unix = cdf_file.varget('time_unix')

    t = unix_to_decimal(t_unix)
    ti = np.where(t == find_nearest(t, t1))[0][0]
    tf = np.where(t == find_nearest(t, t2))[0][0]
    t = t[ti:tf]
    flux = flux[ti:tf]

    #elijo la energia mas baja
    Ei = find_nearest(energia, 20)
    inicio = np.where(energia == Ei)[0][0]

    #empiezo un array al que le voy a meter las energias
    E = Ei
    i = inicio
    E = np.append(E, find_nearest(energia, E+20))
    i = np.append(i, np.where(energia == E[1])[0][0])

    for j in range(len(energia)):
        if E[j+1] > 100:
            break
        else:
            E = np.append(E, find_nearest(energia, E[j+1]+20))
            i = np.append(i, np.where(energia == E[j+2])[0][0])

    E = np.append(E, find_nearest(energia, 220))
    E = np.append(E,find_nearest(energia, 300))
    i = np.append(i, np.where(energia == E[-2])[0][0])
    i = np.append(i, np.where(energia == E[-1])[0][0])

    flux_cut = np.zeros((len(flux), len(i)))

    for j in range(len(i)):
        flux_cut[:, j] = flux[:, int(i[j])]

    #borramos el punto donde el log(0) = inf
    for j in range(len(flux_cut[:,0])):
        for i in range(len(flux_cut[0,:])):
            if np.log(flux_cut[j,i]) < 1:
                flux_cut[j, i] = None

    # plot_select(t, flux_cut, E)
    # for j in range(len(flux_cut[0,:])):
    #     plt.semilogy(t, flux_cut[:,j], label = E[j])
    # plt.legend()
    # plt.xlabel('tiempo')
    # plt.ylabel('flujo')
    # plt.show()

    # ti = float(input("Tiempo inicial = "))
    # tf = float(input("Tiempo final = "))
    #
    # t1 = find_nearest(t, ti)
    # t2 = find_nearest(t, tf)
    # inicio = np.where(t == t1)[0][0]
    # fin = np.where(t == t2)[0][0]
    #
    return(t, flux_cut, E)
