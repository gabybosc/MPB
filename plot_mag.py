import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from funciones import find_nearest
from funciones_plot import plot_datetime

"""
Plotea sólo los datos de MAG de alta resolución
"""

np.set_printoptions(precision=4)

date_entry = input('Enter a date in YYYY-DDD format \n')\
# date_entry = '2016-066'
year, doy = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

ti = float(input("Tiempo inicial = "))
tf = float(input("Tiempo final = "))
n = int(ti*32*3600)

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop
mag = np.loadtxt(path + f'MAG_hires/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=n, usecols=(1,6,7,8,9))

dia = mag[:,0]
t = mag[:,1]  #el dia decimal
t = (t - dia) * 24 #hdec

M = np.size(t) #el numero de datos

#el campo
B = np.zeros((M, 3))
for i in range(2,5):
    B[:,i-2] = mag[:, i]

inicio = np.where(t == find_nearest(t, ti))[0][0]
fin = np.where(t == find_nearest(t, tf))[0][0]

B_norm = np.linalg.norm(B, axis = 1)
B_cut = B_norm[inicio:fin]
t_cut = t[inicio:fin]

plt.figure()
plt.plot(t_cut, B_cut)
plt.xlabel('t (hdec)')
plt.ylabel('|B|')
plt.title('MAG hires hdec')

plt.figure()
plot_datetime(year, month, day, t_cut, B_cut)
# plt.plot(t_cut, B_cut)
plt.xlabel('t (UTC)')
plt.ylabel('|B|')
plt.title('MAG hires UTC')

plt.show(block= False)
