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

from MPB.funciones import fechas

year, month, day, doy = fechas()

# mag_hires = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/highres/{year}/{month}/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts"

mag_1s = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/1sec/{year}/{month}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"

swea = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swea.calibrated/data/svy_spec/{year}/{month}/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf"

# lpw = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.lpw.derived/data/lp-nt/{year}/{month}/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf"
#
swia_onboard = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/onboard_svy_mom/{year}/{month}/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf"

# swica = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/coarse_arc_3d/{year}/{month}/mvn_swi_l2_coarsearc3d_{year}{month}{day}_v01_r01.cdf"

# swifa = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/fine_arc_3d/{year}/{month}/mvn_swi_l2_finearc3d_{year}{month}{day}_v01_r01.cdf"


if gethostname() == "gbosco":
    path = f"../../../../../media/gabybosc/datos/"
else:
    path = "../../../datos/"


p_mag1s = (
    path + f"MAG_1s/{year}/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts"
)
p_maghr = path + f"MAG_hires/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts"
p_swea = path + f"SWEA/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf"
p_swiamom = path + f"SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf"
p_swifa = path + f"SWIA/mvn_swi_l2_finearc3d_{year}{month}{day}_v01_r01.cdf"
p_swica = path + f"SWIA/mvn_swi_l2_coarsearc3d_{year}{month}{day}_v01_r01.cdf"
p_lpw = path + f"LPW/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf"


if not os.path.isfile(p_mag1s):
    with urllib.request.urlopen(mag_1s) as response, open(
        p_mag1s,
        "wb",
    ) as out_file:
        shutil.copyfileobj(response, out_file)


# if not os.path.isfile(p_maghr):
#     with urllib.request.urlopen(mag_hires) as response, open(
#         p_maghr,
#         "wb",
#     ) as out_file:
#         shutil.copyfileobj(response, out_file)


if not os.path.isfile(p_swea):
    with urllib.request.urlopen(swea) as response, open(p_swea, "wb") as out_file:
        shutil.copyfileobj(response, out_file)


if not os.path.isfile(p_swiamom):
    with urllib.request.urlopen(swia_onboard) as response, open(
        p_swiamom,
        "wb",
    ) as out_file:
        shutil.copyfileobj(response, out_file)

# if not os.path.isfile(p_swica):
#     with urllib.request.urlopen(swica) as response, open(
#         p_swica,
#         "wb",
#     ) as out_file:
#         shutil.copyfileobj(response, out_file)


# if not os.path.isfile(p_swifa):
#     with urllib.request.urlopen(swifa) as response, open(
#         p_swifa,
#         "wb",
#     ) as out_file:
#         shutil.copyfileobj(response, out_file)


# if not os.path.isfile(p_lpw):
#     with urllib.request.urlopen(lpw) as response, open(p_lpw, "wb") as out_file:
#         shutil.copyfileobj(response, out_file)
