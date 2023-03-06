import numpy as np
import csv as csv
import pandas as pd
from funciones import donde

# path = "../../Documents/"
cyril = pd.read_csv("CatalogoMAVEN_Cyril.csv")
cyril.rename(
    columns={"yyyy-mm-dd HH:MM:SS": "date", "ThetaBn (deg)": "ThetaBn"}, inplace=True
)
halekas = pd.read_csv("CatalogoMAVEN_Halekas.csv")
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

Todo esto está comentado porque con hacerlo una vez alcanza, está guardado el 
array como lst_c.npy y lst_h.npy
"""

# lst_h = []  # índice de date y timme donde se cumpla la condición
# lst_c = []
# # elijo fecha y hora
# for i in range(len(date_c)):
#     fecha = date_c[i]
#     hora = time_c[i]
#     # busca todos los f dónde está esa fecha en H
#     for idx, f in enumerate(date_h):
#         if f == fecha:
#             # print(idx, f)
#             # print(time_h[idx])
#             # # agarra esos índices y mira la hora
#             HH = int(hora.split(":")[0])
#             hh = int(time_h[idx].split(":")[0])
#             # si la hora (HH) es igual a +-1 a alguna de las horas en estos índices, la toma
#             if HH == hh or HH == hh + 1 or HH == hh - 1:
#                 lst_h.append(idx)
#                 lst_c.append(i)

# np.save("outputs/lst_h.npy", np.array(lst_h))

lst_c = np.load("outputs/lst_c.npy")
lst_h = np.load("outputs/lst_h.npy")

"""
Ahora que tengo la lista de los valores de c y h donde "coinciden" BS y SW
tengo que hacer una subselección
Para el grupo 1 vamos a empezar por elegir noviembre de 2014
"""

anio = "2021"  # input("año\n")
# mes = input("mes (MM)\n")

nov2014 = []  # lista de índices de las listas lstc y lsth de noviembre
for i, val in enumerate(lst_c):
    year = date_c[val].split("-")[0]
    month = date_c[val].split("-")[1]
    for mes in range(5, 10):
        if year == anio and int(month) == mes:
            nov2014.append(i)

# date_c[lst_c[nov2014[0]]] me da la fecha

"""
Ahora que tengo el mes, quiero elegir 10 fechas con Pdyn y Theta variados
primero filtro para quedarme con los que sean quasipara o quasiperp
"""
l = []
for i in range(len(nov2014)):
    c = cyril["ThetaBn"][lst_c[nov2014[i]]]
    if c < 30 or c > 70:  # los que son muy claramente quasipara o quasiperp
        l.append(i)

"""
Si son muchas, elige 10 random y después yo miro y buscando la MPB más linda elijo cuál :)
"""

# vv = np.unique([date_c[lst_c[nov2014[i]]] for i in l])  # las fechas sin repetir
# if len(vv) > 10:
#     vv = np.random.choice(vv, 10)


"""
Escribe un archivo con la fecha, hora, theta y pdyn
"""

with open("outputs/catalogo_grupo4.txt", "a") as file:
    for i in l:
        fec = date_c[lst_c[nov2014[i]]]
        tim = time_c[lst_c[nov2014[i]]]
        the = cyril["ThetaBn"][lst_c[nov2014[i]]]
        pre = halekas["Pdyn"][lst_h[nov2014[i]]]
        file.write(f"{fec}\t{tim}\t{the}\t{pre}")
        file.write("\n")


# with open("outputs/catalogo_grupo2.txt", "a") as file:
#     for i in vv:
#         file.write(f"{i}")
#         file.write("\n")
