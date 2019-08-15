import numpy as np
import matplotlib.pyplot as plt
from funciones_plot import set_axes_equal

def ajuste_conico(posicion, index, orbita, x3, x0 = 0.78, e = 0.9, L = 0.96):
    ##### conica que toma los parámetros de Vignes
    theta = np.linspace(0, np.pi *3/4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    r = L / (1 + e * np.cos(THETA))

    R = posicion[index,:] / 3390 #la posicion de la nave en RM
    ####### Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0,0,0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
    r1 = L0 / (1 + e * np.cos(THETA))

    #ahora plotea
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.set_xlabel(r'$X_{MSO} (R_m)$')
    ax.set_ylabel(r'$Y_{MSO} (R_m)$')
    ax.set_zlabel(r'$Z_{MSO} (R_m)$')
    ax.set_aspect('equal')
    ax.plot(orbita[:,0], orbita[:,1], orbita[:,2], color='green', label='Órbita')
    ax.scatter(R[0], R[1], R[2], label='MAVEN', color='k', s=40)#, marker='x')
    X1 = x0 + r1 * np.cos(THETA)
    Y1 = r1 * np.sin(THETA) * np.cos(PHI)
    Z1 = r1 * np.sin(THETA) * np.sin(PHI)
    plot = ax.plot_surface(
        X1, Y1, Z1, rstride=4, cstride=4, alpha=0.5, edgecolor='none', cmap=plt.get_cmap('Blues_r'))

    asc = L0 / (1 - e**2)  #semieje mayor
    bsc = np.sqrt(asc*L0)
    csc = e*asc - x0 #donde está centrada. Hay que ver el signo

    norm_vignes = np.array([(R[0]+csc)*2 / asc**2, R[1]*2 / (bsc)**2, R[2]*2 / (bsc)**2]) #la normal de vignes
    norm_vignes = norm_vignes/np.linalg.norm(norm_vignes) #normalizado

    ax.quiver(R[0], R[1], R[2], norm_vignes[0], norm_vignes[1], norm_vignes[2], color='b',length=0.5, label='Normal del ajuste') #asi se plotea un vector
    ax.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color='k',length=0.5, label='Normal del MVA')
    # normal_boot = np.array([0.9183, 0.3186, 0.2351])
    # ax.quiver(R[0], R[1], R[2], normal_boot[0], normal_boot[1], normal_boot[2], color='m',length=0.5, label='Normal del bootstrap')

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)
    ax.legend()
    set_axes_equal(ax)

    return(norm_vignes, X1, Y1, Z1)
