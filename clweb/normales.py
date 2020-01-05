import numpy as np
import matplotlib.pyplot as plt
from funciones_plot import set_axes_equal
from plot_orbitas import datos_fijos
from funciones import find_nearest, UTC_to_hdec
from mpl_toolkits.mplot3d import axes3d

"""
Plotea la fig 5 del poster de la AGU2019
"""

def MPB(x0 = 0.78, e = 0.9, L = 0.96):
    theta = np.linspace(0, np.pi *3/4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    r1 = L / (1 + e * np.cos(THETA))

    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)

    #ahora plotea
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.set_xlabel(r'$X_{MSO} (R_m)$', fontsize=14)
    ax.set_ylabel(r'$Y_{MSO} (R_m)$', fontsize=14)
    ax.set_zlabel(r'$Z_{MSO} (R_m)$', fontsize=14)
    ax.set_aspect('equal')
    plot = ax.plot_surface(
        X1, Y1, Z1, rstride=4, cstride=4, alpha=0.5, edgecolor='gray', cmap=plt.get_cmap('Blues_r'))

    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color='#c1440e', linewidth=0.5)
    set_axes_equal(ax)

    # plt.savefig(f'../outputs/figs_MPB/ajuste_{fecha}.png')
    #
    return()

def parametros_elipse(R, x0 = 0.78, e = 0.9, L = 0.96): #R es la posicion del cruce en RM
    ##### conica que toma los parámetros de Vignes
    theta = np.linspace(0, np.pi *3/4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    r = L / (1 + e * np.cos(THETA))

    ####### Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0,0,0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
    r1 = L0 / (1 + e * np.cos(THETA))

    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)

    asc = L0 / (1 - e**2)  #semieje mayor
    bsc = np.sqrt(asc*L0)
    csc = e*asc - x0#donde está centrada. Hay que ver el signo

    return(asc, bsc, csc)

def normal_vignes(R, asc, bsc, csc):
    fig = plt.gcf()
    ax = plt.gca()
    norm_vignes = np.array([(R[0]+csc)*2 / asc**2, R[1]*2 / (bsc)**2, R[2]*2 / (bsc)**2]) #la normal de vignes
    norm_vignes = norm_vignes/np.linalg.norm(norm_vignes) #normalizado
    ax.quiver(R[0], R[1], R[2], norm_vignes[0], norm_vignes[1], norm_vignes[2], color='r', length=0.5) #asi se plotea un vector

    return()

def normal_MVA(R, x3, colour='k', lab = 'MVA normal'):
    fig = plt.gcf()
    ax = plt.gca()
    ax.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color=colour,length=0.5, label=lab)
    return()


MPB()

t, posicion_cut, year, month, day = datos_fijos(2015,10,12,18.75,19.75)
t_medio = UTC_to_hdec('19:19:13')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.956, 0.048, -0.290])
asc,bsc,csc = parametros_elipse(R)

# normal_MVA(R, x3, '2015-oct-12')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2015,10,10,12,13)
t_medio = UTC_to_hdec('12:40:47')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.956, -0.286, 0.070])
asc,bsc,csc = parametros_elipse(R)
# normal_MVA(R, x3, 'b', '2015-oct-10')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2016,'04','05',5,6)
t_medio = UTC_to_hdec('05:16:30')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.815, -0.575, 0.076])
asc,bsc,csc = parametros_elipse(R)
# normal_MVA(R, x3, 'r', '2016-apr-05')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2016,'03',31,12.5,13.5)
t_medio = UTC_to_hdec('13:04:38')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.871, -0.476, -0.117])
asc,bsc,csc = parametros_elipse(R)
# normal_MVA(R, x3, 'g', '2016-mar-31')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2016,'03',16,17.75,18.75)
t_medio = UTC_to_hdec('18:13:50')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.9198,-0.3021,0.2505])
nf = np.array([ 0.8565, -0.0662,  0.5119])
asc,bsc,csc = parametros_elipse(R)

# normal_MVA(R, x3, 'm', '2016-mar-16')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2017,11,24,11.75,12.75)
t_medio = UTC_to_hdec('12:15:27')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.981, -0.032, 0.193])
asc,bsc,csc = parametros_elipse(R)

# normal_MVA(R, x3, 'g', '2017-nov-24')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

plt.show()
