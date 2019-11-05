import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from funciones import find_nearest
from funciones_plot import plot_datetime

"""
Plotea sólo los datos de MAG de baja resolución
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

ti = float(input("Tiempo inicial = "))
tf = float(input("Tiempo final = "))

n = 150

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop

mag = np.loadtxt(path + f'MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=n)

dia = mag[:,1]
t = mag[:,6]  #el dia decimal
t = (t - dia) * 24 #hdec

M = np.size(t) #el numero de datos

#el campo
#el campo
B = np.zeros((M, 3))
for i in range(7,10):
    B[:,i-7] = mag[:, i]

#la posición(x,y,z)
posicion = np.zeros((M, 3))
for i in range(11,14):
    posicion[:,i-11] = mag[:, i]


inicio = np.where(t == find_nearest(t, ti))[0][0]
fin = np.where(t == find_nearest(t, tf))[0][0]

B_norm = np.linalg.norm(B, axis = 1)
B_cut = B_norm[inicio:fin]
t_cut = t[inicio:fin]

plt.figure()
plt.plot(t_cut, B_cut)
plt.xlabel('t (hdec)')
plt.ylabel('|B|')
plt.title('MAG lowres hdec')

plt.figure()
plot_datetime(year, month, day, t_cut, B_cut)
# plt.plot(t_cut, B_cut)
plt.xlabel('t (UTC)')
plt.ylabel('|B|')
plt.title('MAG lowres UTC')

plt.figure()
plt.plot(t[inicio:fin], posicion[inicio:fin, 2])
plt.xlabel('t (hdec)')
plt.ylabel('Z_MSO')
plt.title('Z en función del tiempo')

plt.show(block= False)
