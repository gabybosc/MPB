"""
descarga los archivos de mag low res
El problema es que no puede asegurarse de que existan los archivos antes porque
la página no tira 404. (Si tira 404, acá está la solución
https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head)
"""

import urllib.request
import shutil
import numpy as np
import datetime as dt
import os as os

np.set_printoptions(precision=4)

# year = input('Year\n')
year = 2015
month = 10
doy = []
day = []

for dia in range(1, 30):
    date_orbit = dt.date(int(year), int(month), int(dia))
    day.append(date_orbit.strftime("%d"))
    doy.append(date_orbit.strftime("%j"))
month = date_orbit.strftime("%m")


def return_url(url):
    """como tiene nombre con una variable, mira todos a ver cuál es el link bueno"""
    prefix = "https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/"
    for num in range(10):
        test = prefix + url + f"{num}.sts"
        result = urllib.request.urlopen(test)
        if result.headers["Content-Length"] is None:
            url = test
            break

    return url


def descargar(fname, urlname):
    """
    Si no existe el file (fname) agarra la url (urlname) y lo descarga
    """

    if not os.path.isfile(fname):
        with urllib.request.urlopen(urlname) as response, open(
            fname,
            "wb",
        ) as out_file:
            shutil.copyfileobj(response, out_file)


directory = f"../../../datos/MAG_1s/{year}/"

for i in range(len(doy)):
    mag_1s = return_url(
        f"maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy[i]}ss1s_{year}{month}{day[i]}_v01_r0"
    )

    path = (
        directory + f"mvn_mag_l2_{year}{doy[i]}ss1s_{year}{month}{day[i]}_v01_r01.sts"
    )
    descargar(path, mag_1s)
