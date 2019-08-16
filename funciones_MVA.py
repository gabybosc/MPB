import numpy as np
import matplotlib.pyplot as plt
from funciones_plot import set_axes_equal
from scipy.stats import norm
from matplotlib.mlab import normpdf

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

    return(norm_vignes, X1, Y1, Z1, R)

def plot_velocidades(X1, Y1, Z1, R, norm_vignes, x3, v_media, v_para, v_para_MVA):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1,1,1, projection='3d')
    ax1.set_xlabel('x mso')
    ax1.set_ylabel('y mso')
    ax1.set_zlabel('z mso')
    ax1.set_aspect('equal')
    plot = ax1.plot_surface(
        X1, Y1, Z1, rstride=4, cstride=4, cmap=plt.get_cmap('Blues_r'), alpha=0.5)
    ax1.scatter(R[0], R[1], R[2])

    ax1.quiver(R[0], R[1], R[2], norm_vignes[0], norm_vignes[1], norm_vignes[2], color='g',length=0.5, label='fit normal') #asi se plotea un vector
    ax1.quiver(R[0], R[1], R[2], v_media[0], v_media[1], v_media[2], color='b',length=0.5, label='velocity') #asi se plotea un vector
    ax1.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color='k',length=0.5, label='MVA normal')
    ax1.quiver(R[0], R[1], R[2], v_para[0], v_para[1], v_para[2], color='r',length=0.5, label='v parallel') #asi se plotea un vector
    ax1.quiver(R[0], R[1], R[2], v_para_MVA[0], v_para_MVA[1], v_para_MVA[2], color='m',length=0.5, label='v parallel MVA') #asi se plotea un vector
    ax1.legend()

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    ax1.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)

    set_axes_equal(ax1)

def plot_FLorentz(X1, Y1, Z1, R, J_v, B_upstream, B_downstream, fuerza_mva, x3):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1, projection='3d')
    ax2.set_xlabel(r'$X_{MSO} (R_M)$')
    ax2.set_xlim(left=2, right=0)
    ax2.set_ylabel(r'$Y_{MSO} (R_M)$')
    ax2.set_zlabel(r'$Z_{MSO} (R_M)$')
    ax2.set_aspect('equal')
    ax2.scatter(R[0], R[1], R[2])
    plot = ax2.plot_surface(
        X1, Y1, Z1, rstride=4, cstride=4, alpha=0.5, cmap=plt.get_cmap('Blues_r'))

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    ax2.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)

    ax2.quiver(R[0], R[1], R[2], J_v[0], J_v[1], J_v[2], color='r',length=1E-3, label='Corriente en volumen')
    ax2.quiver(R[0], R[1], R[2], B_upstream[0], B_upstream[1], B_upstream[2], color='b',length=1E-2, label='B upstream')
    ax2.quiver(R[0], R[1], R[2], B_downstream[0], B_downstream[1], B_downstream[2], color='g',length=1E-2, label='B downstream')
    ax2.quiver(R[0], R[1], R[2], fuerza_mva[0], fuerza_mva[1],fuerza_mva[2], color='m',length=2E13, label='Fuerza MVA')
    ax2.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color='k',length=0.5, label='Normal del MVA', linewidths=0.5)

    set_axes_equal(ax2)
    ax2.legend(loc='upper right',bbox_to_anchor=(1.1, 1.05))

def plot_bootstrap(out, out_phi):
    plt.figure()
    plt.subplot(311)
    n, bins, patches = plt.hist(out, 50, density = 1, alpha=0.5)
    (muB, sigmaB) = norm.fit(out)
    y = norm.pdf(bins, muB, sigmaB)
    plt.plot(bins, y)
    plt.xlabel(r'$\langle B_3 \rangle$ (nT)')

    plt.subplot(312)
    n, bins, patches = plt.hist(out_phi[:,0]*57.2958, 50, density=1, alpha=0.5)
    (mu31, sigma31) = norm.fit(out_phi[:,0]*57.2958)
    y = norm.pdf(bins, mu31, sigma31)
    plt.plot(bins, y)
    plt.xlabel(r'$\Delta \phi_{31}$ (º)')

    plt.subplot(313)
    n, bins, patches = plt.hist(out_phi[:,1]*57.2958, 50, density=1, alpha=0.5)
    (mu32, sigma32) = norm.fit(out_phi[:,1]*57.2958)
    y = norm.pdf(bins, mu32, sigma32)
    plt.plot(bins, y)
    plt.xlabel(r'$\Delta \phi_{32}$ (º)')
    plt.tight_layout()

    return(muB, sigmaB, mu31, sigma31, mu32, sigma32)
