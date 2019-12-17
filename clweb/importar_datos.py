import numpy as np
import os
from funciones import find_nearest

def importar_mag(year, month, day, ti, tf):
    path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
    if os.path.isfile(path + 'mag_filtrado.txt'):
        mag = np.loadtxt(path + 'mag_filtrado.txt', skiprows=2)
        M = len(mag[:,0]) #el numero de datos
        B = mag[:, :3]
        Bnorm = mag[:,-1]


    else:
        mag = np.loadtxt(path + 'MAG.asc')
        M = len(mag[:,0]) #el numero de datos
        B = mag[:, 6:9]
        Bnorm = np.linalg.norm(B, axis=1)

    mag = np.loadtxt(path + 'MAG.asc')
    hh = mag[:,3]
    mm = mag[:,4]
    ss = mag[:,5]

    t = hh + mm/60 + ss/3600 #hdec

    posicion = np.zeros((M, 3))
    for i in range(9,12):
        posicion[:,i-9] = mag[:, i]

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    return(mag, t_cut, B, posicion)

    ###############################################################################################SWEA
def importar_swea(year, month, day, ti, tf):
    path = f'../../../datos/clweb/{year}-{month}-{day}/'
    swea = np.loadtxt(path + 'SWEA.asc')

    energy = swea[:, 7]
    JE_total = swea[:, -1]

    t_swea = np.unique(swea[:,3] + swea[:,4]/60 + swea[:,5]/3600) #hdec

    energias = [50 + i*50 for i in range(3)]

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]

    return(swea, t_cut, energias)

    ###############################################################################################SWIA
def importar_swia(year, month, day, ti, tf):
    path = f'../../../datos/clweb/{year}-{month}-{day}/'
    swia = np.loadtxt(path + 'SWIA.asc')

    density = swia[:,-1]

    t_swia = swia[:,3] + swia[:,4]/60 + swia[:,5]/3600 #hdec

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]

    return(swia, t_cut, density)

    ############################################################################################### LPW
def importar_lpw(year, month, day, ti, tf):
    path = f'../../../datos/clweb/{year}-{month}-{day}/'
    lpw = np.loadtxt(path + 'LPW.asc')

    e_density = lpw[:,-1]

    t_lpw = lpw[:,3] + lpw[:,4]/60 + lpw[:,5]/3600

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]

    return(lpw, t_cut, e_density)