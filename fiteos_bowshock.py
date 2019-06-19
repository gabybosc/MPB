import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
#%%
def MAVEN_bowshock(x,y,z):
    #parametros fit MAVEN (Gruesbeck)
    A, B, C, D, E, F, G, H, I = 0.049, 0.157, 0.153, 0.026, 0.012, 0.051, 0.566, -0.031, 0.019
    #A_err, B_err, C_err, D_err, E_err, F_err, G_err, H_err, I_err = 7.8E-6, 1.1E-6, 1.0E-6, 3.8E-6, 2.2E-6, 3.1E-6, 1.0E-6, 2.6E-6, 2.1E-6
    return A*x**2 + B*y**2 + C*z**2 + D*x*y + E*y*z + F*x*z + G*x + H*y + I*z - 1 #defino mi función tal que f(x,y,z)=0

def norm_fit_MAVEN(x,y,z):
    A, B, C, D, E, F, G, H, I = 0.049, 0.157, 0.153, 0.026, 0.012, 0.051, 0.566, -0.031, 0.019
    #A_err, B_err, C_err, D_err, E_err, F_err, G_err, H_err, I_err = 7.8E-6, 1.1E-6, 1.0E-6, 3.8E-6, 2.2E-6, 3.1E-6, 1.0E-6, 2.6E-6, 2.1E-6
    n_fit = np.array([2*A*x + D*y + F*z + G, 2*B*y + D*x + E*z + H, 2*C*z + E*y + F*x + I])/np.linalg.norm(np.array([2*A*x + D*y + F*z + G, 2*B*y + D*x + E*z + H, 2*C*z + E*y + F*x + I]))
    return n_fit


#parametros fit MGS (Vignes)
#X0, epsilon, L = 0.64, 1.03, 2.04
#X0_err, epsilon_err, L_err = 0.02, 0.01, 0.02

def plot_implicit(fn, Rc, n, limites): #para plotear funciones 3D del tipo f(x,y,z) = 0
    xmin, xmax, ymin, ymax, zmin, zmax = limites*3 #limites tiene que ser una tupla (lim_min, lim_max)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, 100) # resolution of the contour
    B = np.linspace(zmin, zmax, 15) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = fn(X,Y,z)
        cset = ax.contour(X, Y, Z+z, [z], zdir='z')
        # [z] defines the only level to plot for this contour for this value of z

    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = fn(X,y,Z)
        cset = ax.contour(X, Y+y, Z, [y], zdir='y')

    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = fn(x,Y,Z)
        cset = ax.contour(X+x, Y, Z, [x], zdir='x')
    
    ax.scatter(Rc[0], Rc[1], Rc[2], color="g", s=100) #plot centro del shock
    
    ax.quiver(Rc[0], Rc[1], Rc[2], n[0], n[1], n[2], length = 2, linewidth = 5, arrow_length_ratio = 0.1, color = 'c', normalize = True, label = 'normal coplanar')
    
    ax.quiver(Rc[0], Rc[1], Rc[2], norm_fit_MAVEN(Rc[0],Rc[1],Rc[2])[0], norm_fit_MAVEN(Rc[0],Rc[1],Rc[2])[1], norm_fit_MAVEN(Rc[0],Rc[1],Rc[2])[2], length = 2, linewidth = 5, arrow_length_ratio = 0.1, color = 'r', normalize = True, label = 'normal fit')    
    
    
    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)
    
    ax.set_xlabel(r'$X_{MSO}$ $[R_M]$', fontsize = 20)
    ax.set_ylabel(r'$Y_{MSO}$ $[R_M]$', fontsize = 20)
    ax.set_zlabel(r'$Z_{MSO}$ $[R_M]$', fontsize = 20)
    plt.tick_params(axis='both', which = 'both', length = 4, width = 2, labelsize = 20)
    plt.legend(loc = 0, fontsize = 20)
    plt.show()
    
#%%
plot_implicit(MAVEN_bowshock, Rc, n, (-2.5, 2.5))
