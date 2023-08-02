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

# lst = np.genfromtxt("../outputs/grupo1.txt", dtype="str", usecols=0)


def return_url(url, extension):
    """como tiene nombre con una variable, mira todos a ver cuál es el link bueno"""
    prefix = "https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/"
    for num in range(10):
        test = prefix + url + f"{num}.{extension}"
        result = urllib.request.urlopen(test)
        if result.headers["Content-Length"] is None:
            url = test
            break

    return url


def descargar(fname, urlname, j, instrumento):
    """
    Si no existe el file (fname) agarra la url (urlname) y lo descarga
    """

    if not os.path.isfile(fname):
        with urllib.request.urlopen(urlname) as response, open(
            fname,
            "wb",
        ) as out_file:
            shutil.copyfileobj(response, out_file)

        # print(f"{instrumento} dia {j} listo")


for grupo in [4, 3, 2, 1]:
    for p in ["para", "perp"]:
        path = f"../bs_mpb/outs_catalogo_previa/grupo{grupo}/"

        fechas = np.load(path + f"fecha_{p}.npy")
        lst = np.unique(fechas)

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

            mag_1s = return_url(
                f"maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r0",
                "sts",
            )
            # mag_hires = return_url(
            #     f"maven.mag.calibrated/data/ss/highres/{year}/{month}/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r0",
            #     "sts",
            # )

            swea = return_url(
                f"maven.swea.calibrated/data/svy_spec/{year}/{month}/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r0",
                "cdf",
            )
            lpw = return_url(
                f"maven.lpw.derived/data/lp-nt/{year}/{month}/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r0",
                "cdf",
            )

            swia_mom = return_url(
                f"maven.swia.calibrated/data/onboard_svy_mom/{year}/{month}/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v02_r0",
                "cdf",
            )

            # swica = return_url(
            #     f"maven.swia.calibrated/data/coarse_arc_3d/{year}/{month}/mvn_swi_l2_coarsearc3d_{year}{month}{day}_v02_r0",
            #     "cdf",
            # )

            # swifa = return_url(
            #     f"maven.swia.calibrated/data/fine_arc_3d/{year}/{month}/mvn_swi_l2_finearc3d_{year}{month}{day}_v02_r0",
            #     "cdf",
            # )

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
            descargar(p_mag1s, mag_1s, j, "mag")
            # descargar(p_maghr, mag_hires, j, "mag")
            descargar(p_swea, swea, j, "swea")
            descargar(p_swiamom, swia_mom, j, "swia")
            # descargar(p_swica, swica, j, "swica")
            # descargar(p_swifa, swifa, j, "swifa")
            # descargar(p_lpw, lpw, j, "lpw")
            print(f"día {j} listo")
