import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import glob as glob

np.set_printoptions(precision=4)

# DATOS DE PDS
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
year = int(input("Year\n"))
path = glob.glob(f"../../../datos/MAG_1s/{year}/*.sts")

for i, j in enumerate(
    path
):  # loop en todos los archivos .sts para cada año. i me da el índice en la lista, j me da el archivo
    # datos = np.loadtxt(path + f'mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=148) #lee todo y me da todo
    datos = np.loadtxt(j, skiprows=160)
    n = 2
    datos = datos[
        :-n, :
    ]  # borra las ultimas 2 filas, que es ya el dia siguiente (no sé si siempre)

    doy = datos[:, 1]
    t = datos[:, 6]  # el dia decimal
    t = (t - doy) * 24  # hdec
    date_orbit = dt.datetime(int(year), 1, 1) + dt.timedelta(
        doy[0] - 1
    )  # para convertir el doty en date

    # year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    # el campo
    B = np.zeros((len(t), 3))
    for i in range(7, 10):
        B[:, i - 7] = datos[:, i]

    # Si quiero elegir manualmente la orbita:
    plt.figure()
    plt.plot(t, np.linalg.norm(B, axis=1))
    plt.xlabel("t (hdec)")
    plt.ylabel("|B|")
    plt.ylim([-1, 80])
    plt.title(f"Orbitas del dia {year}/{month}/{day}")
    plt.show()
