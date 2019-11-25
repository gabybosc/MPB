import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import datetime as dt
from shutil import copyfile
from funciones import datenum
import matplotlib.dates as md
import matplotlib.cm as cm


date_entry = '2017-11-24' #input('Enter a date in YYYY-DDD or YYYY-MM-DD format \n')

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
B = mag[:, 6:9]

Bnorm = mag[:,-1]


T = 12.26 - 12.25
Tseg = T * 3600
ws = 2 * np.pi/Tseg
wp = ws - ws/2
# ws = 2*ws
N, Wn = signal.buttord(wp, ws, 40, 50, False)
b,a = signal.butter(N, Wn,'low', False)
Bnorm_filtrado = signal.filtfilt(b, a, Bnorm)

plt.plot(Bnorm, label='sin filtro')
plt.plot(Bnorm_filtrado,linewidth = 1, label = ws)
plt.legend()
plt.show()

# for w in [ws, 2*ws, 3*ws, 5*ws, 10*ws]:
#
#     N, Wn = signal.buttord(wp, w, 40, 50, False)
#     b,a = signal.butter(N, Wn,'low', False)
#     Bnorm_filtrado = signal.filtfilt(b, a, Bnorm)
#
#     plt.plot(Bnorm, label='sin filtro')
#     plt.plot(Bnorm_filtrado,linewidth = 1, label = w)
#     plt.legend()
#     plt.show()

Bx_filtrado = signal.filtfilt(b, a, B[:,0])
By_filtrado = signal.filtfilt(b, a, B[:,1])
Bz_filtrado = signal.filtfilt(b, a, B[:,2])
happy = 'y'# input('If happy press Y\n')

if happy == 'y' or 'Y':
    with open(path + 'mag_filtrado.txt','w') as file:
        file.write(f'Los datos de MAG filtrados para frecuencia wp = {wp}, ws = {ws}.\n')
        file.write(f'Bx  By  Bz  B.\n')
        for i in range(M):
            file.write(f'{Bx_filtrado[i]}\t{By_filtrado[i]}\t{Bz_filtrado[i]}\t{Bnorm_filtrado[i]}\t')
            file.write('\n')
