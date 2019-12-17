import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, MultiCursor
from matplotlib.mlab import normpdf
from scipy.stats import norm
import cdflib as cdf
from datetime import datetime
from funciones import find_nearest, deltaB, unix_to_decimal, datenum, Bpara_Bperp, fechas, tiempos
from importar_datos import importar_mag, importar_lpw, importar_swea, importar_swia
import matplotlib.dates as md
import datetime as dt
import matplotlib.cm as cm
import pdb
from cycler import cycler
import matplotlib as mpl

np.set_printoptions(precision=4)

"""
Este script plotea mag, swea, swia y lpw en la región de interés
Es la fig principal del grl.

"""

year, month, day, doy = fechas()
ti, tf = tiempos()

inbound = input('inbound? y/n\n')


path = f'../../../datos/clweb/{year}-{month}-{day}/'
datos = np.loadtxt('../outputs/t1t2t3t4.txt')
for j in range(len(datos)):
    if datos[j,0] == float(year) and datos[j,1] == float(doy) and int(datos[j,2]) == int(ti):
        i = j

mag, t, B, posicion = importar_mag(year, month, day)
if inbound == 'y':
    t1 = datos[i,2]
    t2 = datos[i,3]
    t3 = datos[i,4]
    t4 = datos[i,5]
    tiempos = np.array([t1,t2,t3,t4])
    t_up = t1 - 0.015
    t_down = t4 + 0.015
    j_inicial = np.where(t == find_nearest(t, t1-0.075))[0][0]
    j_final = np.where(t == find_nearest(t, t4+0.075))[0][0]
    t_cut = t[j_inicial:j_final]

else:
    t1 = datos[i,5]
    t2 = datos[i,4]
    t3 = datos[i,3]
    t4 = datos[i,2]
    tiempos = np.array([t1,t2,t3,t4])
    t_down = t4 - 0.015
    t_up = t1 + 0.015
    j_inicial = np.where(t == find_nearest(t, t4-0.075))[0][0]
    j_final = np.where(t == find_nearest(t, t1+0.075))[0][0]
    t_cut = t[j_inicial:j_final]


Bnorm = np.linalg.norm(B, axis=1)

mag_low = np.loadtxt(path + 'mag_1s.sts', skiprows=160)
tlow = mag_low[:,6]  #el dia decimal
tlow = (tlow -int( doy)) * 24 #para que me de sobre la cantidad de horas

Mlow = np.size(tlow) #el numero de datos
#el campo
Blow = np.zeros((Mlow, 3))
for i in range(7,10):
    Blow[:,i-7] = mag_low[:, i]


B_para, B_perp_norm, t_plot = Bpara_Bperp(Blow, tlow, t[0], t[-1])

########## SWEA

swea, t_swea, energias = importar_swea(year, month, day)

energy = swea[:,7]
JE_total = swea[:,-1]
if inbound == 'y':
    inicio_swea = np.where(t_swea == find_nearest(t_swea, t1-0.075))[0][0]
    fin_swea = np.where(t_swea == find_nearest(t_swea, t4+0.075))[0][0]
else:
    inicio_swea = np.where(t_swea == find_nearest(t_swea, t4-0.075))[0][0]
    fin_swea = np.where(t_swea == find_nearest(t_swea, t1+0.075))[0][0]

######## densidad SWIA

swia, t_swia, density = importar_swia(year, month, day)
if inbound == 'y':
    inicio_swia = np.where(t_swia == find_nearest(t_swia, t1-0.075))[0][0]
    fin_swia = np.where(t_swia == find_nearest(t_swia, t4+0.075))[0][0]
else:
    inicio_swia = np.where(t_swia == find_nearest(t_swia, t4-0.075))[0][0]
    fin_swia = np.where(t_swia == find_nearest(t_swia, t1+0.075))[0][0]

####### densidad electrones
lpw, t_lpw, e_density = importar_lpw(year, month, day)

if inbound == 'y':
    inicio_lpw = np.where(t_lpw == find_nearest(t_lpw, t1-0.075))[0][0]
    fin_lpw = np.where(t_lpw == find_nearest(t_lpw, t4+0.075))[0][0]
else:
    inicio_lpw = np.where(t_lpw == find_nearest(t_lpw, t4-0.075))[0][0]
    fin_lpw = np.where(t_lpw == find_nearest(t_lpw, t1+0.075))[0][0]

############# tiempos UTC
year = int(year)
month = int(month)
day = int(day)

tiempo_mag = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_cut]) #datenum es una función mía
tiempo_swea = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swea[inicio_swea:fin_swea]])
tiempo_swia = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_swia[inicio_swia:fin_swia]])
tiempo_lpw = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_lpw[inicio_lpw:fin_lpw]])
tiempo_low = np.array([np.datetime64(datenum(year, month, day, x)) for x in t_plot]) #datenum es una función mía

tm1 = np.where(t_cut == find_nearest(t_cut, t1))
tm2 = np.where(t_cut == find_nearest(t_cut, t2))
tm3 = np.where(t_cut == find_nearest(t_cut, t3))
tm4 = np.where(t_cut == find_nearest(t_cut, t4))
tm_up = np.where(t_cut == find_nearest(t_cut, t_up))
tm_down = np.where(t_cut == find_nearest(t_cut, t_down))

tiempo_lim = [tiempo_mag[tm1],tiempo_mag[tm2],tiempo_mag[tm3],tiempo_mag[tm4]]

###### funciones para el plot que se repiten

def colores_in(ax, tiempo_mag, tm1, tm4, tm_up, tm_down):
    ax.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[tm1][0], facecolor = '#581845', alpha = 0.2)
    ax.axvspan(xmin = tiempo_mag[tm_up][0], xmax = tiempo_mag[tm1][0], facecolor = '#581845', alpha = 0.5)
    ax.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[tm4][0], facecolor = '#C70039', alpha = 0.35)
    ax.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[-1], facecolor = '#FFC300', alpha = 0.2)
    ax.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[tm_down][0], facecolor = '#FFC300', alpha = 0.5)

def colores_out(ax, tiempo_mag, tm1, tm4, tm_up, tm_down):
    ax.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[tm4][0], facecolor = '#FFC300', alpha = 0.2)
    ax.axvspan(xmin = tiempo_mag[tm_down][0], xmax = tiempo_mag[tm4][0], facecolor = '#FFC300', alpha = 0.5)
    ax.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[tm1][0], facecolor = '#C70039', alpha = 0.35)
    ax.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[-1], facecolor = '#581845', alpha = 0.2)
    ax.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[tm_up][0], facecolor = '#581845', alpha = 0.5)

#####plot grl

mpl.rcParams['axes.prop_cycle'] = cycler('color',['#5F021F','#336699', '#C70039', '#00270F'])
fig = plt.figure(1, figsize=(8,25))#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.12,right=0.95, hspace = 0.0, wspace=0.15)
plt.xticks( rotation=25 )
xfmt = md.DateFormatter('%H:%M')

axz1 = plt.gca()


axz1.xaxis.set_major_formatter(xfmt)
axz1 = plt.subplot2grid((4,1),(0,0))
axz1.plot(tiempo_mag, Bnorm[j_inicial:j_final], label = '|B|')
axz1.plot(tiempo_mag, B[j_inicial:j_final,0], label='Bx MSO')
axz1.plot(tiempo_mag, B[j_inicial:j_final,1], label='By MSO')
axz1.plot(tiempo_mag, B[j_inicial:j_final,2],label='Bz MSO')
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
axz1.set_ylabel('B components and |B| (nT)')
axz1.set_title(f'MAVEN MAG SWEA LPW SWIA {year}-{month}-{day}')

axz2 = plt.subplot2grid((4,1),(1,0), sharex=axz1)
axz2.xaxis.set_major_formatter(xfmt)
axz2.semilogy(tiempo_low, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / |B|')
axz2.semilogy(tiempo_low, B_perp_norm, linewidth=1, label=r'|$\Delta B \perp$| / |B|')
axz2.set_ylabel('Relative variation \n of B')
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)

axz3 = plt.subplot2grid((4,1),(2,0), sharex=axz1)
axz3.xaxis.set_major_formatter(xfmt)
axz3.set_ylabel('Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)')
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[0]
    JE = JE_total[index]
    plt.semilogy(tiempo_swea, JE[inicio_swea:fin_swea], label = f'{energia} eV')
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)

axz4 = plt.subplot2grid((4,1),(3,0), sharex=axz1)
axz4.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_swia, density[inicio_swia:fin_swia])
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
axz4.set_ylabel('SW proton \n density (cm⁻³)')
axz4.set_xlabel('Time (UTC)')


if inbound == 'y':
    for ax in [axz1,axz2,axz3]:
        colores_in(ax, tiempo_mag, tm1, tm4, tm_up, tm_down)
        plt.setp(ax.get_xticklabels(), visible=False)
    axz4.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[tm1][0], facecolor = '#581845', alpha = 0.2, label = 'Magnetosheath')
    axz4.axvspan(xmin = tiempo_mag[tm_up][0], xmax = tiempo_mag[tm1][0], facecolor = '#581845', alpha = 0.5, label = 'Upstream')
    axz4.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[tm4][0], facecolor = '#C70039', alpha = 0.35, label = 'MPB')
    axz4.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[tm_down][0], facecolor = '#FFC300', alpha = 0.5, label = 'Downstream')
    axz4.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[-1], facecolor = '#FFC300', alpha = 0.2, label = 'MPR')
else:
    for ax in [axz1,axz2,axz3]:
        colores_out(ax, tiempo_mag, tm1, tm4, tm_up, tm_down)
        plt.setp(ax.get_xticklabels(), visible=False)
    axz4.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[tm1][0], facecolor = '#FFC300', alpha = 0.2, label = 'MPR')
    axz4.axvspan(xmin = tiempo_mag[tm_down][0], xmax = tiempo_mag[tm4][0], facecolor = '#FFC300', alpha = 0.5, label = 'Downstream')
    axz4.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[tm1][0], facecolor = '#C70039', alpha = 0.35, label = 'MPB')
    axz4.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[tm_up][0], facecolor = '#581845', alpha = 0.5, label = 'Upstream')
    axz4.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[-1], facecolor = '#581845', alpha = 0.2, label = 'Magnetosheath')

for ax in [axz1,axz2,axz3, axz4]:
    ax.set_xlim(tiempo_mag[0],tiempo_mag[-1])
    ax.grid()
    ax.legend()

# plt.tight_layout()
plt.show(block=False)

#########La fig de la tesis sin LPW
mpl.rcParams['axes.prop_cycle'] = cycler('color',['#5F021F','#336699', '#C70039', '#00270F'])
fig = plt.figure(2, figsize=(8,30))#Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.12,right=0.95, hspace = 0.0, wspace=0.15)
plt.xticks( rotation=25 )
xfmt = md.DateFormatter('%H:%M')

axz1 = plt.gca()

axz1.xaxis.set_major_formatter(xfmt)
axz1 = plt.subplot2grid((5,1),(0,0))
axz1.plot(tiempo_mag, B[j_inicial:j_final,0], label='Bx MSO')
axz1.plot(tiempo_mag, B[j_inicial:j_final,1], label='By MSO')
axz1.plot(tiempo_mag, B[j_inicial:j_final,2], label='Bz MSO')
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
axz1.set_ylabel('B components (nT)')
axz1.set_title(f'MAVEN MAG SWEA LPW SWIA {year}-{month}-{day}')

axz2 = plt.subplot2grid((5,1),(1,0), sharex=axz1)
axz2.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_mag, Bnorm[j_inicial:j_final])
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
plt.ylabel('|B| (nT)')

axz3 = plt.subplot2grid((5,1),(3,0), sharex=axz1)
axz3.xaxis.set_major_formatter(xfmt)
axz3.set_ylabel('Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)')
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[0]
    JE = JE_total[index]
    plt.semilogy(tiempo_swea, JE[inicio_swea:fin_swea], label = f'{energia} eV')
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)

# axz4 = plt.subplot2grid((6,1),(4,0), sharex=axz1)
# axz4.xaxis.set_major_formatter(xfmt)
# plt.semilogy(tiempo_lpw, e_density[inicio_lpw:fin_lpw])
# for xc in tiempo_lim:
#     plt.axvline(x = xc, color = 'k', linewidth=1.5)
# axz4.grid()
# axz4.set_ylabel('Total electron \n density (cm⁻³)')

axz5 = plt.subplot2grid((5,1),(4,0), sharex=axz1)
axz5.xaxis.set_major_formatter(xfmt)
plt.plot(tiempo_swia, density[inicio_swia:fin_swia])
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
axz5.set_ylabel('SW proton \n density (cm⁻³)')
axz5.set_xlabel('Time (UTC)')

axz6 = plt.subplot2grid((6,1),(2,0), sharex=axz1)
axz6.xaxis.set_major_formatter(xfmt)
axz6.plot(tiempo_low, B_para, linewidth=1, label=r'|$\Delta B \parallel$| / |B|')
axz6.plot(tiempo_low, B_perp_norm, '-.', linewidth=1, label=r'|$\Delta B \perp$| / |B|')
for xc in tiempo_lim:
    plt.axvline(x = xc, color = 'k', linewidth=1.5)
axz6.set_ylabel('Relative variation \n of B')



if inbound == 'y':
    print('in')
    for ax in [axz1,axz2,axz3, axz4, axz6]:
        colores_in(ax, tiempo_mag, tm1, tm4, tm_up, tm_down)
        plt.setp(ax.get_xticklabels(), visible=False)
    axz5.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[tm1][0], facecolor = '#581845', alpha = 0.2, label = 'Magnetosheath')
    axz5.axvspan(xmin = tiempo_mag[tm_up][0], xmax = tiempo_mag[tm1][0], facecolor = '#581845', alpha = 0.5, label = 'Upstream')
    axz5.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[tm4][0], facecolor = '#C70039', alpha = 0.35, label = 'MPB')
    axz5.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[tm_down][0], facecolor = '#FFC300', alpha = 0.5, label = 'Downstream')
    axz5.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[-1], facecolor = '#FFC300', alpha = 0.2, label = 'MPR')
else:
    print('out')
    for ax in [axz1,axz2,axz3, axz4, axz6]:
        colores_out(ax, tiempo_mag, tm1, tm4, tm_up, tm_down)
        plt.setp(ax.get_xticklabels(), visible=False)
    axz5.axvspan(xmin = tiempo_mag[0], xmax = tiempo_mag[tm1][0], facecolor = '#FFC300', alpha = 0.2, label = 'MPR')
    axz5.axvspan(xmin = tiempo_mag[tm_down][0], xmax = tiempo_mag[tm4][0], facecolor = '#FFC300', alpha = 0.5, label = 'Downstream')
    axz5.axvspan(xmin = tiempo_mag[tm4][0], xmax = tiempo_mag[tm1][0], facecolor = '#C70039', alpha = 0.35, label = 'MPB')
    axz5.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[tm_up][0], facecolor = '#581845', alpha = 0.5, label = 'Upstream')
    axz5.axvspan(xmin = tiempo_mag[tm1][0], xmax = tiempo_mag[-1], facecolor = '#581845', alpha = 0.2, label = 'Magnetosheath')


for ax in [axz1,axz2,axz3, axz4, axz5, axz6]:
    ax.set_xlim(tiempo_mag[0],tiempo_mag[-1])
    ax.grid()
    ax.legend()

# plt.tight_layout()
plt.show(block=False)
