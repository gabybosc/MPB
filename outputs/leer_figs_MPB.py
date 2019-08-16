import pickle
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

fechas = np.loadtxt('../t1t2t3t4.txt')
for i in range(len(fechas)):
# for i in range(1):
    year = int(fechas[i,0])
    doy = int(fechas[i, 1])
    t1 = fechas[i,2]
    t2 = fechas[i,3]
    t3 = fechas[i,4]
    t4 = fechas[i,5]
    tiempos = np.array([t1,t2,t3,t4])

    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date
    year = date_orbit.strftime("%Y")
    doy = date_orbit.strftime("%j")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    fig = pickle.load(open(f'figs_MPB/MPB_y{year}_d{doy}_t{int(t1)}.pkl','rb'))
    ax_master = fig.axes[0]
    for ax in fig.axes:
        if ax is not ax_master:
            ax_master.get_shared_x_axes().join(ax_master, ax)

    plt.show()
