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

    #ahora plotea
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.set_xlabel(r'$X_{MSO} (R_m)$', fontsize=14)
    ax.set_ylabel(r'$Y_{MSO} (R_m)$', fontsize=14)
    ax.set_zlabel(r'$Z_{MSO} (R_m)$', fontsize=14)
    ax.set_aspect('equal')
    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)
    plot = ax.plot_surface(
        X1, Y1, Z1, rstride=4, cstride=4, alpha=0.5, edgecolor='gray', cmap=plt.get_cmap('Blues_r'))

    asc = L / (1 - e**2)  #semieje mayor
    bsc = np.sqrt(asc*L)
    csc = e*asc - x0 #donde está centrada. Hay que ver el signo

    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color='#c1440e', linewidth=0.5)
    set_axes_equal(ax)

    # plt.savefig(f'../outputs/figs_MPB/ajuste_{fecha}.png')
    #
    return(asc, bsc, csc)

def normal_vignes(R, asc, bsc, csc):
    fig = plt.gcf()
    ax = plt.gca()
    norm_vignes = np.array([(R[0]+csc)*2 / asc**2, R[1]*2 / (bsc)**2, R[2]*2 / (bsc)**2]) #la normal de vignes
    norm_vignes = norm_vignes/np.linalg.norm(norm_vignes) #normalizado
    ax.quiver(R[0], R[1], R[2], norm_vignes[0], norm_vignes[1], norm_vignes[2], color='r', length=0.5, label='Fit normal') #asi se plotea un vector

    return()

def normal_MVA(R, x3, colour='k', lab = 'MVA normal'):
    fig = plt.gcf()
    ax = plt.gca()
    ax.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color=colour,length=0.5, label=lab)
    return()

asc, bsc, csc = MPB()

# dates = ['2015-10-12-18.75-19.75-19:19:13', '2015-10-10-12-13', '2016-04-05-5-6', '2016-03-31-12.5-13.5','2016-03-16-17.75-18.75', '2017-11-24-11.75-12.75']
#
# for date in dates:
#     year, month, day, ti, tf, tm = date.split('-')
t, posicion_cut, year, month, day = datos_fijos(2015,10,12,18.75,19.75)
t_medio = UTC_to_hdec('19:19:13')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.956, 0.048, -0.290])

# normal_MVA(R, x3, '2015-oct-12')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2015,10,10,12,13)
t_medio = UTC_to_hdec('12:40:47')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.956, -0.286, 0.070])

# normal_MVA(R, x3, 'b', '2015-oct-10')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2016,'04','05',5,6)
t_medio = UTC_to_hdec('05:16:30')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.815, -0.575, 0.076])

# normal_MVA(R, x3, 'r', '2016-apr-05')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2016,'03',31,12.5,13.5)
t_medio = UTC_to_hdec('13:04:38')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.871, -0.476, -0.117])

# normal_MVA(R, x3, 'g', '2016-mar-31')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2016,'03',16,17.75,18.75)
t_medio = UTC_to_hdec('18:13:50')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.920,-0.302,0.251])

# normal_MVA(R, x3, 'm', '2016-mar-16')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

t, posicion_cut, year, month, day = datos_fijos(2017,11,24,11.75,12.75)
t_medio = UTC_to_hdec('12:15:27')
index = np.where(t == find_nearest(t, t_medio))[0][0]
R = posicion_cut[index,:] #la posicion de la nave en RM
x3 = np.array([0.981, -0.032, 0.193])

# normal_MVA(R, x3, 'g', '2017-nov-24')
normal_MVA(R, x3)
normal_vignes(R, asc, bsc,csc)

# plt.legend()
plt.show()
