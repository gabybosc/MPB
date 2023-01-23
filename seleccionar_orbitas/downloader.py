"""
descarga los archivos de mag, lpw, swea y swia de la lista de fechas del año que le pida.
chequea si ya lo tengo descargado para no hacerlo dos veces
"""

import urllib.request
import shutil
import numpy as np
import datetime as dt
from socket import gethostname
import os as os

# fechas = np.arange(61, 92)  # np.loadtxt("outputs/fechas_MVA_2016.txt")

lst = np.genfromtxt("../outputs/hoja_grupo2.txt", dtype="str", usecols=0)

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

    """
    Las URLS
    """
    prefix = "https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/"
    mag_1s = (
        prefix
        + f"maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"
    )

    # mag_hires = prefix + f"maven.mag.calibrated/data/ss/highres/{year}/{month}/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts"

    # como tiene nombre con una variable, mira todos a ver cuál es el link bueno
    for num in range(10):
        swea_test = (
            prefix
            + f"maven.swea.calibrated/data/svy_spec/{year}/{month}/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r0{num}.cdf"
        )
        result = urllib.request.urlopen(swea_test)
        if result.headers["Content-Length"] is None:
            swea = swea_test
            break

    for num in range(10):
        lpw_test = (
            prefix
            + f"maven.lpw.derived/data/lp-nt/{year}/{month}/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r0{num}.cdf"
        )
        result_lpw = urllib.request.urlopen(lpw_test)
        if result_lpw.headers["Content-Length"] is None:
            lpw = lpw_test
            break

    swia_mom = (
        prefix
        + f"maven.swia.calibrated/data/onboard_svy_mom/{year}/{month}/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v02_r00.cdf"
    )

    # swica = prefix + f"maven.swia.calibrated/data/coarse_arc_3d/{year}/{month}/mvn_swi_l2_coarsearc3d_{year}{month}{day}_v02_r00.cdf"

    # swifa = prefix + f"maven.swia.calibrated/data/fine_arc_3d/{year}/{month}/mvn_swi_l2_finearc3d_{year}{month}{day}_v02_r00.cdf"
    """
    Los archivos
    """
    p_mag1s = (
        path
        + f"MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"
    )
    # p_maghr = (
    #     path + f"MAG_hires/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts"
    # )
    p_swea = path + f"SWEA/mvn_swe_l2_svyspec_{year}{month}{day}.cdf"
    p_swiamom = path + f"SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}.cdf"
    # p_swifa = path + f"SWIA/mvn_swi_l2_finearc3d_{year}{month}{day}.cdf"
    # p_swica = path + f"SWIA/mvn_swi_l2_coarsearc3d_{year}{month}{day}.cdf"
    p_lpw = path + f"LPW/mvn_lpw_l2_lpnt_{year}{month}{day}.cdf"

    """
    Si el archivo no existe, lo descarga
    """
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
