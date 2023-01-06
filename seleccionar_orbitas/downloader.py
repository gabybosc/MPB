"""
descarga los archivos de mag hi res (pues si llegue hasta acá es porque ya vi los
low res), lpw, swea y swia de la lista de fechas del año que le pida.
El problema es que no puede asegurarse de que existan los archivos antes porque
la página no tira 404. (Si tira 404, acá está la solución
https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head)
"""

import urllib.request
import shutil
import numpy as np
import datetime as dt
from socket import gethostname
import os as os

# fechas = np.arange(61, 92)  # np.loadtxt("outputs/fechas_MVA_2016.txt")

lst = np.genfromtxt("../outputs/catalogo_grupo1.txt", dtype="str")

# a = [int(x[-2:]) for x in lst]
# fechas = np.unique(a)

for j in lst:
    year, month, day = j.split("-")

    date_orbit = dt.datetime(int(year), int(month), int(day))

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    if gethostname() == "gbosco":
        path = f"../../../../../media/gabybosc/datos/"
    else:
        path = "../../../datos/"

    mag_1s = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"

    # mag_hires = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/highres/{year}/{month}/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts"

    swea = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swea.calibrated/data/svy_spec/{year}/{month}/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf"

    lpw = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.lpw.derived/data/lp-nt/{year}/{month}/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf"

    swia_mom = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/onboard_svy_mom/{year}/{month}/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf"

    # swica = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/coarse_arc_3d/{year}/{month}/mvn_swi_l2_coarsearc3d_{year}{month}{day}_v01_r01.cdf"

    # swifa = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/fine_arc_3d/{year}/{month}/mvn_swi_l2_finearc3d_{year}{month}{day}_v01_r01.cdf"

    p_mag1s = (
        path
        + f"MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"
    )
    # p_maghr = (
    #     path + f"MAG_hires/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts"
    # )
    p_swea = path + f"SWEA/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf"
    p_swiamom = path + f"SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf"
    # p_swifa = path + f"SWIA/mvn_swi_l2_finearc3d_{year}{month}{day}_v01_r01.cdf"
    # p_swica = path + f"SWIA/mvn_swi_l2_coarsearc3d_{year}{month}{day}_v01_r01.cdf"
    p_lpw = path + f"LPW/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf"

    if not os.path.isfile(p_mag1s):
        with urllib.request.urlopen(mag_1s) as response, open(
            p_mag1s,
            "wb",
        ) as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"mag dia {j} listo")

    # if not os.path.isfile(p_maghr):
    #     with urllib.request.urlopen(mag_hires) as response, open(
    #         p_maghr,
    #         "wb",
    #     ) as out_file:
    #         shutil.copyfileobj(response, out_file)
    #     print(f"mag dia {j} listo")

    if not os.path.isfile(p_swea):
        with urllib.request.urlopen(swea) as response, open(p_swea, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"swea dia {j} listo")

    if not os.path.isfile(p_swiamom):
        with urllib.request.urlopen(swia_mom) as response, open(
            p_swiamom,
            "wb",
        ) as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"swia (onboard) dia {j} listo")

    # if not os.path.isfile(p_swica):
    #     with urllib.request.urlopen(swica) as response, open(
    #         p_swica,
    #         "wb",
    #     ) as out_file:
    #         shutil.copyfileobj(response, out_file)
    #     print(f"swica dia {j} listo")

    # if not os.path.isfile(p_swifa):
    #     with urllib.request.urlopen(swifa) as response, open(
    #         p_swifa,
    #         "wb",
    #     ) as out_file:
    #         shutil.copyfileobj(response, out_file)
    #     print(f"swifa dia {j} listo")

    if not os.path.isfile(p_lpw):
        with urllib.request.urlopen(lpw) as response, open(p_lpw, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"lpw dia {j} listo")
