import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from funciones import find_nearest
from funciones_plot import plot_datetime

"""
Plotea sólo los datos de MAG de alta resolución
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

# ti = float(input("Tiempo inicial = "))
# tf = float(input("Tiempo final = "))
# n = int(ti*32*3600)

path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
mag = np.loadtxt(path + 'MAG.asc')
M = len(mag[:,0]) #el numero de datos
B = mag[:, 6:9]
Bnorm = np.linalg.norm(B, axis=1)
Bxyz_paraperp = mag[:,6:9]

hh = mag[:,3]
mm = mag[:,4]
ss = mag[:,5]

t = hh + mm/60 + ss/3600 #hdec

# inicio = np.where(t == find_nearest(t, ti))[0][0]
# fin = np.where(t == find_nearest(t, tf))[0][0]
#
B_norm = np.linalg.norm(B, axis = 1)
# B_cut = B_norm[inicio:fin]
# t_cut = t[inicio:fin]

plt.figure()
plt.plot(t, B_norm, linewidth=0.5)
plt.xlabel('t (hdec)')
plt.ylabel('|B|')
plt.title('MAG hires hdec')

# plt.figure()
# plot_datetime(year, month, day, t_cut, B_cut)
# # plt.plot(t_cut, B_cut)
# plt.xlabel('t (UTC)')
# plt.ylabel('|B|')
# plt.title('MAG hires UTC')

plt.show(block= False)
