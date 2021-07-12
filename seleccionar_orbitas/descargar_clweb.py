import sys
import os
import requests
from requests.auth import HTTPBasicAuth

sys.path.append("..")

from funciones import fechas

"""
Esto sólo sirve si uso lo que está en el save llamado "todo.cl", o sea, si
se guardan las cosas en ese orden.
"""

year, month, day, doy = fechas()

MAG = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj3.asc"
SWICA = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj4.asc"
SWIFA = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj5.asc"
SWEA = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj6_Panel01_MAVEN_SWEA_P12_SC1.asc"
LPW = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj7.asc"
STATIC = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj8_Panel01_MAVEN_STATIC_P12_SC1.asc"

path = f"../../../datos/clweb/{year}-{month}-{day}/"

if not os.path.exists(path):  # crea el directorio si no existe
    os.makedirs(path)

names = ['MAG', 'SWICA', 'SWIFA', 'SWEA', 'LPW', 'STATIC']
lst = [MAG, SWICA, SWIFA, SWEA, LPW, STATIC]
for i, site_url in enumerate(lst):
    s = requests.Session()
    print('request')
    user = 'boscoboinik'
    password = "Maven2014"
    response = s.get(site_url, auth=HTTPBasicAuth(user, password))
    print('response')

    open(path + f"{names[i]}.asc", "wb").write(response.content)
    print('guardar')
