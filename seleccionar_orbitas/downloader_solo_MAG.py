"""
descarga los archivos de mag low res
El problema es que no puede asegurarse de que existan los archivos antes porque la página no tira 404. (Si tira 404, acá está la solución https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head)
"""

import urllib.request
import shutil
import os as os
import numpy as np
import datetime as dt

np.set_printoptions(precision=4)

# year = input('Year\n')
# fechas = np.loadtxt(f'outputs/fechas_buenas_{year}.txt', skiprows = 1)
# lista = []
# for j in range(len(fechas)):
#     lista = np.append(lista, int(fechas[j, 0]))
#
# lista = set(lista) #me saca los repetidos

# directory = f'../../../../../media/gabybosc/datos/MAG_1s/{year}/'
# if not os.path.exists(directory):
#     os.makedirs(directory)

directory = "../../../datos/MAG_1s/"

year = 2015

for doy in [283, 285]:  # range(0,366):
    date_orbit = dt.datetime(int(year), 1, 1) + dt.timedelta(
        int(doy) - 1
    )  # para convertir el doty en date

    # no puedo pasarle así year (más allá de que sería redundante) porque flashea y se me va al año anterior.
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    mag_1s = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"

    with urllib.request.urlopen(mag_1s) as response, open(
        directory + f"mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts", "wb"
    ) as out_file:
        shutil.copyfileobj(response, out_file)
    print(f"{int(doy) * 100/366:.2f}%")
