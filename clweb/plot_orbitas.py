import numpy as np
import matplotlib.pyplot as plt
from funciones import fechas, find_nearest_inicial,find_nearest_final, tiempos

"""
Hace la fig2 del poster de la AGU2019
"""

def datos():
    year, month, day, doy = fechas()
    ti, tf = tiempos()
    path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
    mag = np.loadtxt(path + 'MAG.asc')

    M = len(mag[:,0]) #el numero de datos

    hh = mag[:,3]
    mm = mag[:,4]
    ss = mag[:,5]

    t = hh + mm/60 + ss/3600 #hdec
    t_utc = [f'{int(hh[j])}:{int(mm[j])}' for j in range( len(t))]

    posicion = np.zeros((M, 3))
    for i in range(9,12):
        posicion[:,i-9] = mag[:, i]/3390

    inicio = np.where(t == find_nearest_inicial(t, ti))[0][0]
    fin = np.where(t == find_nearest_final(t, tf))[0][0]

    posicion_cut = posicion[inicio:fin,:]
    t_cut = t[inicio:fin]

    return(t_cut, posicion_cut, year, month, day)

def datos_fijos(year, month, day, ti, tf):
    path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
    mag = np.loadtxt(path + 'MAG.asc')

    M = len(mag[:,0]) #el numero de datos

    hh = mag[:,3]
    mm = mag[:,4]
    ss = mag[:,5]

    t = hh + mm/60 + ss/3600 #hdec
    t_utc = [f'{int(hh[j])}:{int(mm[j])}' for j in range( len(t))]

    posicion = np.zeros((M, 3))
    for i in range(9,12):
        posicion[:,i-9] = mag[:, i]/3390

    inicio = np.where(t == find_nearest_inicial(t, ti))[0][0]
    fin = np.where(t == find_nearest_final(t, tf))[0][0]

    posicion_cut = posicion[inicio:fin,:]
    t_cut = t[inicio:fin]

    return(t_cut, posicion_cut, year, month, day)

def BS_MPB(L, e, x0):
    THETA = np.linspace(0, np.pi * 3/4, 100)
    PHI = np.linspace(0, 2 * np.pi, 100)

    r1 = L / (1 + e * np.cos(THETA))

    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)
    yz = np.sqrt(Y1**2+Z1**2)
    return(X1,yz)

def marte(x_bs, yz_bs, x_mpb, y_mpb):
    fig, ax = plt.subplots()
    ax.plot()
    ax.plot(x_bs, yz_bs, color='#07aec7', linestyle='-.')
    ax.plot(x_mpb, yz_mpb, color='#FF1493', linestyle='-.')
    ax.axis('equal')
    ax.set_xlim(0,3)
    ax.set_ylim(0,2.5)
    circle = plt.Circle((0, 0), 1, color='#c1440e', clip_on=True)
    ax.add_artist(circle)
    ax.set_title('MAVEN MSO coordinates', fontsize=16)
    ax.set_xlabel(r'$X_{MSO}$ ($R_M$)', fontsize=14)
    ax.set_ylabel(r'$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)', fontsize=14)

def orbitas(posicion, year, month, day):
    x = posicion[:,0]
    y = posicion[:,1]
    z = posicion[:,2]
    proyeccion = np.sqrt(y**2+z**2)
    fig = plt.gcf()
    ax = fig.gca()
    ax.plot(x, proyeccion, label=f'{day}-{month}-{year}')


x_bs, yz_bs = BS_MPB(2.04, 1.03, 0.64)
x_mpb, yz_mpb = BS_MPB(0.96, 0.9, 0.78)
marte(x_bs, yz_bs, x_mpb, yz_mpb)


# t, posicion_cut, year, month, day = datos()
t, posicion_cut, year, month, day = datos_fijos(2015,10,10,12,13)
orbitas(posicion_cut, year, month, day)

t, posicion_cut, year, month, day = datos_fijos(2015,10,12,18.75,19.75)
orbitas(posicion_cut, year, month, day)

t, posicion_cut, year, month, day = datos_fijos(2016,'03',16,17.5,18.5)
orbitas(posicion_cut, year, month, day)

t, posicion_cut, year, month, day = datos_fijos(2016,'03',31,12.5,13.5)
orbitas(posicion_cut, year, month, day)

t, posicion_cut, year, month, day = datos_fijos(2016,'04','05',16,17)
orbitas(posicion_cut, year, month, day)

t, posicion_cut, year, month, day = datos_fijos(2017,11,24,12,13)
orbitas(posicion_cut, year, month, day)

plt.legend()
plt.show(block=False)
plt.savefig(f'../outputs/figs_MPB/Orbitas.png', dpi=200)
