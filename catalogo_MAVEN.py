import numpy as np
import csv as csv
import pandas as pd
from funciones import donde

path = "../../Documents/"
cyril = pd.read_csv(path + "CatalogoMAVEN_Cyril.csv")
cyril.rename(
    columns={"yyyy-mm-dd HH:MM:SS": "date", "ThetaBn (deg)": "ThetaBn"}, inplace=True
)
halekas = pd.read_csv(path + "CatalogoMAVEN_Halekas.csv")
halekas.rename(columns={"Unnamed: 0": "date", "Pdyn_proton": "Pdyn"}, inplace=True)

date_c = [fecha[:10] for fecha in cyril["date"]]
time_c = [fecha[11:] for fecha in cyril["date"]]

date_h = [fecha[:10] for fecha in halekas["date"]]
time_h = [fecha[11:] for fecha in halekas["date"]]

"""
Las fechas de Cyril marcan el BS
Las fechas de Halekas marcan el SW
tiene que hacer lo siguiente:
1) elegir una fecha y hora de cyril
2) buscar todas las repeticiones de esa fecha en H
3) ver cuándo coincide la hora de C con la hora de H entre esas repeticiones
"""

lst_h = []  # índice de date y timme donde se cumpla la condición
lst_c = []
# elijo fecha y hora
for i in range(len(date_c)):
    fecha = date_c[i]
    hora = time_c[i]
    # busca todos los f dónde está esa fecha en H
    for idx, f in enumerate(date_h):
        if f == fecha:
            # print(idx, f)
            # print(time_h[idx])
            # # agarra esos índices y mira la hora
            HH = int(hora.split(":")[0])
            hh = int(time_h[idx].split(":")[0])
            # si la hora (HH) es igual a +-1 a alguna de las horas en estos índices, la toma
            if HH == hh or HH == hh + 1 or HH == hh - 1:
                lst_h.append(idx)
                lst_c.append(i)


"""
Ahora que tengo la lista de los valores de c y h donde "coinciden" BS y SW
tengo que hacer una subselección
Para el grupo 1 vamos a empezar por elegir noviembre de 2014
"""

nov2014 = []  # lista de índices de las listas lstc y lsth de noviembre
for i, val in enumerate(lst_c):
    year = date_c[val].split("-")[0]
    month = date_c[val].split("-")[1]
    if year == "2014" and month == "11":
        nov2014.append(i)

# date_c[lst_c[nov2014[0]]] me da la fecha

"""
Ahora que tengo el mes, quiero elegir 10 fechas con Pdyn y Theta variados
primero filtro para quedarme con los que sean quasipara o quasiperp
"""
l = []
for i in range(len(nov2014)):
    c = cyril["ThetaBn"][lst_c[nov2014[i]]]
    if c < 20 or c > 70:  # los que son muy claramente quasipara o quasiperp
        l.append(i)


"""
Podría elegir 10 valores random, que seguramente sean suficientemente espaciados en P.
Eso o mirar los 17 y elegir según la MPB más linda
"""

# if len(l) > 10:
#     for i in l:
#         h = float(halekas["Pdyn"][lst_c[nov2014[i]]])
#         if
