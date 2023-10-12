import os
import requests
from requests.auth import HTTPBasicAuth

import sys

sys.path.append("..")

from funciones import fechas

"""
Esto sólo sirve si uso lo que está en el save llamado "VEX.cl", o sea, si
se guardan las cosas en ese orden.
"""

year, month, day, doy = fechas()

MAG = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj5.asc"
posicion = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj6.asc"
ELS = (
    "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj3_Panel01_ASPERA4_ELS_P48_SC1.asc"
)
IMA = (
    "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj4_Panel01_ASPERA4_IMA_P24_SC1.asc"
)

path = "../../../datos/VEX/"

if not os.path.exists(path):  # crea el directorio si no existe
    os.makedirs(path)

names = [
    f"MAG{year}{month}{day}",
    f"pos{year}{month}{day}",
    f"ELS{year}{month}{day}",
    f"IMA{year}{month}{day}",
]
lst = [MAG, posicion, ELS, IMA]
for i, site_url in enumerate(lst):
    s = requests.Session()
    print("request")
    user = "boscoboinik"
    password = "Maven2014"
    response = s.get(site_url, auth=HTTPBasicAuth(user, password))
    print("response")

    open(path + f"{names[i]}.asc", "wb").write(response.content)
    print("guardar")
