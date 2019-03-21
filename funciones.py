from __future__ import print_function
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import RectangleSelector
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from mpl_toolkits.mplot3d import axes3d
from matplotlib.image import AxesImage
from numpy.random import rand
import datetime as dt
import calendar

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

def deltat():
    global x1
    global x2
    return(x1, x2)

def deltaB(B):
    B_medio = np.mean(B, axis=0)
    abs_deltaB_para = np.abs(np.dot(B - B_medio, B_medio)) / np.linalg.norm(B_medio)**2   #|deltaB_para / B|

    dot = np.dot(B - B_medio, B_medio / np.linalg.norm(B_medio))
    N = np.zeros((len(dot), len(B_medio)))

    for i in range(len(N)):
        N[i,:] = dot[i] * B_medio/np.linalg.norm(B_medio)
        deltaB_perp = (B - B_medio) - N
        #y ahora necesito el valor abs de perp
        abs_deltaB_perp = np.abs(deltaB_perp) / np.linalg.norm(B_medio)

        return(abs_deltaB_para, abs_deltaB_perp)

def error(lamb, B, M, x, dia):
    phi = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i == j:
                phi[i,j] = 0
            else:
                phi[i,j] = np.sqrt(lamb[2] / (M - 1) * (lamb[i] + lamb[j] - lamb[2])/(lamb[i] - lamb[j])**2) #en radianes

    # print('Matriz de incerteza angular (grados): \n{}'.format(phi  *  57.2958))
    # np.savetxt('phi_%d'%dia[0], phi)

    delta_B3 = np.sqrt(lamb[2] / (M-1) + (phi[2,1] * np.dot(np.mean(B,0), x[1]) )**2 + (phi[2,0] * np.dot(np.mean(B,0), x[0]))**2)
    # print('El error en B3 (nT) = {0:1.3g}'.format(delta_B3))
    return(phi, delta_B3)


#para encontrar el valor mas cercano
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin() #argmin me da el valor del minimo en el array
    return array[idx]

def find_nearest_inicial(array,value): #el valor más cercano pero más grande
    idx = (np.abs(array-value)).argmin() #el indice de la minima dif entre el array y el valor
    if array[idx] < value: #si es menor que el valor que busco:
        idx = np.abs((array - (value + 1/3600))).argmin() #el indice va a ser la min dif entre el array y el valor+1seg
    return array[idx]

def find_nearest_final(array,value): #el valor más cercano pero más chico
    idx = (np.abs(array-value)).argmin() #el indice de la minima dif entre el array y el valor
    if array[idx] > value: #si es mayor que el valor que busco:
        idx = np.abs((array - (value - 1/3600))).argmin() #el indice va a ser la min dif entre el array y el valor-1seg
    return array[idx]

def hodograma(B1, B2, B3, unidad, titulo): #, dia):
    f, (ax1, ax2) = plt.subplots(1, 2, sharex=True) #tienen el mismo eje x
    ax1.plot(B2, B1)
    ax2.plot(B3, B1)
    ax1.scatter(B2[0], B1[0], marker='o', color='r', label = 'inicio')
    ax1.scatter(B2[-1], B1[-1], marker='x',color='r', label = 'fin')
    ax2.scatter(B3[0], B1[0], marker='o', color='r',label = 'inicio')
    ax2.scatter(B3[-1], B1[-1], marker='x',color='r', label = 'fin')
    ax1.set_xlabel('B2 ({})'.format(unidad))
    ax2.set_xlabel('B3 ({})'.format(unidad))
    ax1.set_ylabel('B1 ({})'.format(unidad))
    ax2.set_ylabel('B1 ({})'.format(unidad))
    ax1.grid()
    ax2.grid()
    plt.suptitle(titulo)
    plt.legend()
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.savefig('figuras/hodograma_%d.png'%dia[0])

def line_select_callback(eclick, erelease):
    global x1, x2
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    # t1, t2 = x1, x2

def Mij(B):
    Mij = np.zeros((3,3))
    for i in range(3): #para las tres coordenadas
        for j in range(3):
            Mij[i,j] = np.mean(B[:,i] * B[:,j]) - np.mean(B[:,i]) * np.mean(B[:,j])
    return Mij

def onpick1(event):
    if isinstance(event.artist, Line2D):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        print('X='+str(np.take(xdata, ind)[0])) # Print X point
        # print('Y='+str(np.take(ydata, ind)[0])) # Print Y point

def plot_rect(x, y):
    fig, current_ax = plt.subplots()
    plt.plot(x, y)
    plt.xlabel('Tiempo')
    plt.ylabel('|B| (nT)')
    plt.title('Seleccionar ')
    print("\n      click  -->  release")
    toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # don't use middle button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels',
                                           interactive=True)
    plt.connect('key_press_event', toggle_selector)
    # plt.show()

def plot_select(x,y, E):
    fig, ax = plt.subplots()
    ax.set_title('Seleccionar puntos', picker=True)
    ax.set_xlabel('Tiempo (h)')
    ax.set_ylabel('Flujo (cm⁻² sr⁻¹ s⁻¹)', picker=True)#, bbox=dict(facecolor='red'))
    # line, = ax.semilogy(x, y,picker=5)
    for j in range(len(y[0,:])):
        line = ax.semilogy(x, y[:,j], label='{0:1.4g} eV'.format(E[j]), picker = 5)
    ax.set_ylim(1e4, 4*1e9)
    ax.legend()
    fig.canvas.mpl_connect('pick_event', onpick1)

def rms(x):
    rms = np.sqrt(np.vdot(x, x)/x.size)
    return rms

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
        if event.key in ['A', 'a'] and not toggle_selector.RS.active:
            print(' RectangleSelector activated.')
            toggle_selector.RS.set_active(True)

def unix_to_decimal(t_unix):
    t = np.zeros(np.size(t_unix))
    for i in range(np.size(t)):
        u = datetime.utcfromtimestamp(int(t_unix[i]))
        t[i] = u.second/3600 + u.minute/60 + u.hour #tiempo decimal
    return(t)
    #t me da el mismo array que tengo en el MVA!


###### estas tres (getrem, datenum y araydatenum) las saqué de internet random, son para pasar de hora decimal a datetime
def getrem(input):
    "this function yields the value behind the decimal point"
    import numpy as np
    output=abs(input-np.fix(input))
    return output

def datenum(Yr,Mo=1,Da=1,Hr=0,Mi=0,Se=0,Ms=0):
    "this function works as regular datetime.datetime, but allows for float input"


    #correct faulty zero input
    if Mo<1:
        Mo+=1
    if Da<1:
        Da+=1

    #distribute the year fraction over days
    if  getrem(Yr)>0:
        if calendar.isleap(np.floor(Yr)):
            fac=366
        else:
            fac=365
        Da=Da+getrem(Yr)*fac
        Yr=int(Yr)
    #if months exceeds 12, pump to years
    while int(Mo)>12:
        Yr=Yr+1
        Mo=Mo-12
    #distribute fractional months to days
    if getrem(Mo)>0:
        Da=Da+getrem(Mo)*calendar.monthrange(Yr,int(Mo))[1]
        Mo=int(Mo)
    #datetime input for 28 days always works excess is pumped to timedelta
    if Da>28:
        extraDa=Da-28
        Da=28
    else:
        extraDa=0
    # sometimes input is such that you get 0 day or month values, this fixes this anomaly
    if int(Da)==0:
       Da+=1
    if int(Mo)==0:
       Mo+=1

    #datetime calculation
    mytime=dt.datetime(int(Yr),int(Mo),int(Da))+dt.timedelta(days=extraDa+getrem(Da),hours=Hr,minutes=Mi,seconds=Se,microseconds=Ms)
    return mytime

def araydatenum(*args):
    mydatetimes=[datenum(*[a.squeeze()[x] for a in args]) for x in range(len(args[0].squeeze()))]
    return mydatetimes
