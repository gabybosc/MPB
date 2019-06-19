import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import glob as glob

np.set_printoptions(precision=4)

#DATOS DE PDS
# date_entry = input('Enter a date in YYYY-DDD format \n')
# year, doy = map(int, date_entry.split('-'))
# date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doy en date
#
# year = date_orbit.strftime("%Y")
# month = date_orbit.strftime("%m")
# day = date_orbit.strftime("%d")
# doy = date_orbit.strftime("%j")

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos
# path = '../../datos/MAG_1s/'
path = glob.glob('../../datos/MAG_1s/subsolares/*.sts')

for i,j in enumerate(path): #loop en todos los archivos .sts para cada año. i me da el índice en la lista, j me da el archivo
    # datos = np.loadtxt(path + f'mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=148) #lee todo y me da todo
    datos = np.loadtxt(j, skiprows = 160)
    n =2
    datos = datos[:-n, :] #borra las ultimas 2 filas, que es ya el dia siguiente (no sé si siempre)

    dia = datos[:,1]
    t = datos[:,6]  #el dia decimal
    t = (t - dia) * 24 #hdec
    year = datos[0,0]

    M = np.size(t) #el numero de datos

    #tengo que asegurarme de que no haya agujeros en mis datos
    # for i in range(M-1):
    #     if t[i+1] - t[i] > 24 * 1.5e-5: #1.1e-5 es 1s, le doy un poco más
    #         print('salto en la linea {} de {} segundos'.format(i+144, (t[i+1] - t[i]) / (24*1.1e-5)))

    #el campo
    B = np.zeros((M, 3))
    for i in range(7,10):
        B[:,i-7] = datos[:, i]

    #la posición(x,y,z)
    posicion = np.zeros((M, 3))
    for i in range(11,14):
        posicion[:,i-11] = datos[:, i]

    #la matriz diaria:
    MD = np.zeros((M, 9))
    MD[:, 0] = t
    for i in range(1,4):
        MD[:, i] = B[:,i-1]
    MD[:,4] = np.sqrt(B[:,0]**2 + B[:,1]**2 + B[:,2]**2)#el modulo de B
    for i in range(5,8):
        MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
    MD[:, 8] = np.sqrt(posicion[:,0]**2 + posicion[:,1]**2 + posicion[:,2]**2) - 3390 #altitud en km

    #Si quiero elegir manualmente la orbita:
    plt.figure()
    plt.plot(t, MD[:,4])
    plt.xlabel('t (hdec)')
    plt.ylabel('|B|')
    plt.ylim([-1, 50])
    plt.title(f'Orbitas del dia {dia[0]}')

plt.show()
