import os
import requests
from requests.auth import HTTPBasicAuth

"""
Esto es para descargar los momentos de SWIFA y SWICA del clweb según T_n_v_SWIFA/SWICA
"""

instrumento = input("Instrumento (SWICA/SWIFA)\n")

T_full = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj3.asc"
T_corte = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj4.asc"
n_full = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj5.asc"
n_corte = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj6.asc"
v_full = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj7.asc"
v_corte = "http://clweb.irap.omp.eu/cl/boscoboinik/work/obj8.asc"

path = f"../../../datos/clweb/2016-03-16/momentos/{instrumento}/"

if not os.path.exists(path):  # crea el directorio si no existe
    os.makedirs(path)

names = [
    "temp_full",
    "temp_corte",
    "dens_full",
    "dens_corte",
    "vel_full",
    "vel_corte",
]
lst = [T_full, T_corte, n_full, n_corte, v_full, v_corte]
for i, site_url in enumerate(lst):
    s = requests.Session()
    print("request")
    user = "boscoboinik"
    password = "Maven2014"
    response = s.get(site_url, auth=HTTPBasicAuth(user, password))
    print("response")

    open(path + f"{names[i]}.asc", "wb").write(response.content)
    print("guardar")
