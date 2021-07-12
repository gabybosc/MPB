"""
descarga los archivos de mag low res
El problema es que no puede asegurarse de que existan los archivos antes porque
la página no tira 404. (Si tira 404, acá está la solución
https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head)
Se fija antes de bajarlos que pesen más de un par de kb.
"""

import urllib.request
import shutil
import numpy as np
import datetime as dt
import sys

sys.path.append("..")

from funciones import fechas

np.set_printoptions(precision=4)
#
# meses = {
#     "Jan": "01",
#     "Feb": "02",
#     "Mar": "03",
#     "Apr": "04",
#     "May": "05",
#     "Jun": "06",
#     "Jul": "07",
#     "Aug": "08",
#     "Sep": "09",
#     "Oct": "10",
#     "Nov": "11",
#     "Dec": "12",
# }
#
# shocks = [
#     "25 Dec 2014",
#     "04 Jan 2015",
#     "21 Jan 2015",
#     "03 Feb 2017",
#     "27 Dec 2017",
#     "01 Aug 2018",
# ]

year = int(input("year\n"))

shocks = np.loadtxt(f"calendario_BS_MPB_{year}.txt")

directory = "../../../datos/"


# for i in range(len(shocks)):
#     doy = shocks[i]
#     date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)
#
#     month = date_orbit.strftime("%m")
#     day = date_orbit.strftime("%d")
#     doy = date_orbit.strftime("%j")
#     for j in range(10):
#         mag_1s = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r0{j}.sts"
#
#         req = urllib.request.Request(mag_1s, method="HEAD")
#         f = urllib.request.urlopen(req)
#         f.status
#         size = int(f.headers["Content-Length"])
#
#         if size > 9e5:
#             with urllib.request.urlopen(mag_1s) as response, open(
#                 directory
#                 + f"MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts",
#                 "wb",
#             ) as out_file:
#                 shutil.copyfileobj(response, out_file)
#             print(f"mag dia {doy} es r0{j}")
#             break


for i in range(len(shocks)):
    doy = shocks[i]
    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)

    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    for j in range(10):  # el loop sobre swea
        swea = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swea.calibrated/data/svy_spec/{year}/{month}/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r0{j}.cdf"

        req = urllib.request.Request(swea, method="HEAD")
        f = urllib.request.urlopen(req)
        f.status
        size = int(f.headers["Content-Length"])

        if size > 9e5:
            with urllib.request.urlopen(swea) as response, open(
                directory + f"SWEA/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf",
                "wb",
            ) as out_file:
                shutil.copyfileobj(response, out_file)
            break

    for k in range(10):  # el loop sobre swia
        swia = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/onboard_svy_mom/{year}/{month}/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r0{k}.cdf"
        # onboard porque swica pesa mucho!

        req = urllib.request.Request(swia, method="HEAD")
        f = urllib.request.urlopen(req)
        f.status
        size = int(f.headers["Content-Length"])
        print(size)

        if size > 9e5:
            with urllib.request.urlopen(swia) as response, open(
                directory
                + f"SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf",
                "wb",
            ) as out_file:
                shutil.copyfileobj(response, out_file)
            break


# for i in range(len(shocks)):
#     doy = shocks[i]
#     date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)
#
#     month = date_orbit.strftime("%m")
#     day = date_orbit.strftime("%d")
#     doy = date_orbit.strftime("%j")
#     for j in range(10):
#
#         static = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.static.c/data/c6_32e64m/{year}/{month}/mvn_sta_l2_c6-32e64m_{year}{month}{day}_v02_r0{j}.cdf"
#
#         req = urllib.request.Request(static, method="HEAD")
#         f = urllib.request.urlopen(req)
#         f.status
#         size = int(f.headers["Content-Length"])
#
#         if size > 9e5:
#             with urllib.request.urlopen(static) as response, open(
#                 directory
#                 + f"STATIC/mvn_sta_l2_c6-32e64m_{year}{month}{day}_v01_r01.cdf",
#                 "wb",
#             ) as out_file:
#                 shutil.copyfileobj(response, out_file)
#             print(f"STATIC {doy} es r0{j} con size {size}")
#             break
