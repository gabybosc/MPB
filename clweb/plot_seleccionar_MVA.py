import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, MultiCursor
from matplotlib.colors import LogNorm
import cdflib as cdf
import datetime as dt
from funciones import find_nearest, deltaB, unix_to_decimal, unix_to_timestamp, Bpara_Bperp
from funciones_plot import onpick1
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import matplotlib.dates as md
import matplotlib.cm as cm
import os as os


"""
Le paso los datos de clweb y elijo los tiempos t1t2t3t4
"""


np.set_printoptions(precision=4)

date_entry = input('Enter a date in YYYY-DDD or YYYY-MM-DD format \n')\

if len(date_entry.split('-')) < 3:
    year, doy = map(int, date_entry.split('-'))
    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date
else:
    year, month, day = map(int, date_entry.split('-'))
    date_orbit = dt.date(year, month, day)

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")


path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
mag = np.loadtxt(path + 'mag.asc')

hh = mag[:,3]
mm = mag[:,4]
ss = mag[:,5]

t = hh + mm/60 + ss/3600 #hdec

M = np.size(t) #el numero de datos

#el campo
B = np.zeros((M, 3))
for i in range(6,9):
    B[:,i-6] = mag[:, i]

Bnorm = mag[:,-1]

#Si quiero elegir manualmente la orbita:
plt.plot(t, Bnorm)
plt.xlabel('t (hdec)')
plt.ylabel('|B|')
plt.title('Orbitas')
plt.show()

ti = float(input("Tiempo inicial = "))
tf = float(input("Tiempo final = "))

while tf < ti:
    print('El tiempo inicial no puede ser mayor al final')
    ti = float(input("Tiempo inicial = "))
    tf = float(input("Tiempo final = "))


B_para, B_perp_norm, j_inicial, j_final = Bpara_Bperp(B, t, ti, tf)
t_plot = t[j_inicial+12:j_final+12]

###############################################################################################SWEA

swea = np.loadtxt(path + 'diff_en_flux.asc')

flux = swea[:,-1]
energia = swea[:,7]

t_swea = swea[:,3] + swea[:,4]/60 + swea[:,5]/3600 #hdec

flux_plot = np.transpose(flux)[::-1]

###############################################################################################SWIA

swia = np.loadtxt(path + 'swia_density.asc')

density = swia[:,-1]

t_swia = swia[:,3] + swia[:,4]/60 + swia[:,5]/3600 #hdec


############################################################################################### LPW
lpw = np.loadtxt(path + 'lpw_density.asc')

e_density = lpw[:,-1]

t_lpw = lpw[:,3] + lpw[:,4]/60 + lpw[:,5]/3600




index = np.array((int(year), int(day)))


happy = False
while not happy:
    val = []
    while len(val) < 4:
        plt.clf()#clear figure
        fig = plt.figure(1, constrained_layout=True)#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
        fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.05,right=0.95, hspace = 0.005, wspace=0.15)
        plt.title('Spacebar when ready to click:')

        ax1 = plt.subplot2grid((3,2),(0,0))
        plt.plot(t_plot, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / B')
        plt.plot(t_plot, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / B')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r'|$\Delta B$|/ B')
        ax1.grid()
        ax1.legend()

        ax4 = plt.subplot2grid((3,2),(1,0), sharex=ax1)
        ax4.plot(t, B[:,0], label='Bx')
        ax4.plot(t, B[:,1], label='By')
        ax4.plot(t, B[:,2], label='Bz')
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.set_ylabel('Bx, By, Bz (nT)')
        ax4.legend()
        ax4.grid()

        ax3 = plt.subplot2grid((3,2),(2,0), sharex=ax1)
        plt.plot(t, Bnorm)
        ax3.grid()
        ax3.set_ylabel('|B| (nT)')
        ax3.set_xlabel('Tiempo (hdec)')

        ax5 = plt.subplot2grid((3,2),(0,1), sharex=ax1)
        # ax5.set_ylabel('Energia', picker=True)#, bbox=dict(facecolor='red'))
        # plt.setp(ax5.get_xticklabels(), visible=False)
        # im = plt.imshow(flux_plot, aspect = 'auto',origin = 'lower', extent=(t_swea[0], t_swea[-1],  energia[-1], energia[0]), cmap='inferno', norm=LogNorm(vmin=1E4, vmax=1E9))
        # divider = make_axes_locatable(ax5)
        # cax = divider.append_axes("top", size="7%", pad="1%")
        # cb = plt.colorbar(im, cax=cax, orientation="horizontal")
        # cax.xaxis.set_ticks_position("top")


        ax7 = plt.subplot2grid((3,2),(1,1), sharex=ax1)
        plt.setp(ax7.get_xticklabels(), visible=False)
        ax7.set_ylabel('Densidad de p+ \n del SW (cm⁻³)')
        plt.plot(t_swia, density)
        ax7.grid()

        ax6 = plt.subplot2grid((3,2),(2,1), sharex=ax1)
        ax6.set_ylabel('Densidad total \n de e- (cm⁻³)')
        ax6.set_xlabel('Tiempo (hdec)')

        plt.semilogy(t_lpw, e_density)
        ax6.grid()

        fig.canvas.mpl_connect('pick_event', onpick1)
        multi = MultiCursor(fig.canvas, (ax1, ax3,ax4,ax5,ax6,ax7), color='black', lw=1)

        zoom_ok = False
        print('\nSpacebar when ready to click:\n')
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print('Click to select MPB: ')
        val = np.asarray(plt.ginput(4))[:,0]
        print('Selected values: ', val)
        outs = np.concatenate((index, val))
        outs = sorted(outs[2:6])

    print('Happy? Keyboard click for yes, mouse click for no.')
    happy = plt.waitforbuttonpress()

plt.show(block=False)

with open('../outputs/t1t2t3t4.txt','a') as file:
    file.write('\n')
    file.write(f'{year}\t{doy}\t')
    for k in outs:
        file.write('{0:1.7}\t'.format(k))
