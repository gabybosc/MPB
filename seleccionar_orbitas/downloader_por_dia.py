"""
descarga los archivos de mag hi res (pues si llegue hasta acá es porque ya vi los low res),
lpw, swea y swia de la fecha que quiera.
El problema es que no puede asegurarse de que existan los archivos antes porque
la página no tira 404. Si tira 404, acá está la solución:
https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head
creo que sería mejor que baje el mes entero, si bien ocupa más espacio y tiempo,
ya que no necesita el nombre exacto del archivo.
A veces el problema es que el archivo se llama "_r02" o algo así y por eso no
sirve.
"""

import urllib.request
import shutil
import os
import sys
from socket import gethostname


p = os.path.abspath("../..")
if p not in sys.path:
    sys.path.append(p)  # esto es para poder llamar a funciones desde esta carpeta

from funciones import fechas

year, month, day, doy = fechas()


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


mag_1s = return_url(
    f"maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r0",
    "sts",
)
mag_hires = return_url(
    f"maven.mag.calibrated/data/ss/highres/{year}/{month}/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r0",
    "sts",
)

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

swica = return_url(
    f"maven.swia.calibrated/data/coarse_arc_3d/{year}/{month}/mvn_swi_l2_coarsearc3d_{year}{month}{day}_v02_r0",
    "cdf",
)

swifa = return_url(
    f"maven.swia.calibrated/data/fine_arc_3d/{year}/{month}/mvn_swi_l2_finearc3d_{year}{month}{day}_v02_r0",
    "cdf",
)

if gethostname() == "gbosco":
    path = f"../../../../../media/gabybosc/datos/"
else:
    path = "../../../datos/"


p_mag1s = (
    path + f"MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"
)
p_maghr = path + f"MAG_hires/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts"
p_swea = path + f"SWEA/mvn_swe_l2_svyspec_{year}{month}{day}.cdf"
p_swiamom = path + f"SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}.cdf"
p_swifa = path + f"SWIA/mvn_swi_l2_finearc3d_{year}{month}{day}.cdf"
p_swica = path + f"SWIA/mvn_swi_l2_coarsearc3d_{year}{month}{day}.cdf"
p_lpw = path + f"LPW/mvn_lpw_l2_lpnt_{year}{month}{day}.cdf"


descargar(p_mag1s, mag_1s)
# descargar(p_maghr, mag_hires)
descargar(p_swea, swea)
descargar(p_swiamom, swia_mom)
# descargar(p_swica, swica)
# descargar(p_swifa, swifa)
descargar(p_lpw, lpw)
