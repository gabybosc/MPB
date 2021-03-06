import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
from funciones import find_nearest, unix_to_decimal
from funciones_plot import plot_select
import datetime as dt
import matplotlib.cm as cm
import matplotlib.dates as md

# np.set_printoptions(precision=4)
def flujo_energia(t1, t2, cdf_file):

    # np.set_printoptions(precision=4)

    path = '../../datos/SWEA/' #path a los datos
    cdf_file = cdf.CDF(path + 'mvn_swe_l2_svyspec_20160402_v04_r01.cdf')

    flux = cdf_file.varget('diff_en_fluxes')
    energia = cdf_file.varget('energy')
    t_unix = cdf_file.varget('time_unix')

    tu = unix_to_decimal(t_unix)
    ti = np.where(tu == find_nearest(tu, t1))[0][0]
    tf = np.where(tu == find_nearest(tu, t2))[0][0]
    t = tu[ti:tf]
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

    t1 = t_unix[np.where(tu == find_nearest(tu, 18.2193))[0][0]]
    t2 = t_unix[np.where(tu == find_nearest(tu, 18.2272))[0][0]]
    t3 = t_unix[np.where(tu == find_nearest(tu, 18.2331))[0][0]]
    t4 = t_unix[np.where(tu == find_nearest(tu, 18.2476))[0][0]]

    n = len(flux_cut)
    t_inicial = t_unix[ti]
    t_final = t_unix[tf]
    timestamps = np.linspace(t_inicial, t_final, n)
    dates = [dt.datetime.utcfromtimestamp(ts) for ts in timestamps] #me lo da en UTC
    datenums = md.date2num(dates)

    # flux = cdf_file.varget('diff_en_fluxes')
    # energia = cdf_file.varget('energy')
    # t_unix = cdf_file.varget('time_unix')
    #
    # tu = unix_to_decimal(t_unix)
    # ti = np.where(tu == find_nearest(tu, t1))[0][0]
    # tf = np.where(tu == find_nearest(tu, t2))[0][0]
    # t = tu[ti:tf]
    # flux = flux[ti:tf]
    #
    #
    # #elijo la energia mas baja
    # Ei = find_nearest(energia, 20)
    # inicio = np.where(energia == Ei)[0][0]
    #
    # #empiezo un array al que le voy a meter las energias
    # E = Ei
    # i = inicio
    # E = np.append(E, find_nearest(energia, E+20))
    # i = np.append(i, np.where(energia == E[1])[0][0])
    #
    # for j in range(len(energia)):
    #     if E[j+1] > 100:
    #         break
    #     else:
    #         E = np.append(E, find_nearest(energia, E[j+1]+20))
    #         i = np.append(i, np.where(energia == E[j+2])[0][0])
    #
    # E = np.append(E, find_nearest(energia, 220))
    # E = np.append(E,find_nearest(energia, 300))
    # i = np.append(i, np.where(energia == E[-2])[0][0])
    # i = np.append(i, np.where(energia == E[-1])[0][0])
    #
    # flux_cut = np.zeros((len(flux), len(i)))
    #
    # for j in range(len(i)):
    #     flux_cut[:, j] = flux[:, int(i[j])]
    #
    # #borramos el punto donde el log(0) = inf
    # for j in range(len(flux_cut[:,0])):
    #     for i in range(len(flux_cut[0,:])):
    #         if np.log(flux_cut[j,i]) < 1:
    #             flux_cut[j, i] = None


    return(t, flux_cut, E)
