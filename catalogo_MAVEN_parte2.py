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

catalogo = np.genfromtxt("outputs/catalogo_grupo1.txt", dtype="str")

date_c = [fecha[:10] for fecha in cyril["date"]]
time_c = [fecha[11:] for fecha in cyril["date"]]

date_h = [fecha[:10] for fecha in halekas["date"]]
time_h = [fecha[11:] for fecha in halekas["date"]]

"""
Agarra el sub catálogo que preseleccioné con la primera parte
Me dice qué hora mirar para cada día
"""

l = []  # lista de índices de las listas lstc y lsth de noviembre
for cat in catalogo:
    for i in range(len(date_c)):
        fecha = date_c[i]
        hora = time_c[i]
        # busca todos los f dónde está esa fecha en H
        if cat == fecha:
            for idx, f in enumerate(date_h):
                # agarra esos índices y mira la hora
                if f == cat:
                    HH = int(hora.split(":")[0])
                    hh = int(time_h[idx].split(":")[0])
                    # si la hora (HH) es igual a +-1 a alguna de las horas en estos índices, la toma
                    if HH == hh or HH == hh + 1 or HH == hh - 1:
                        # print("hh y HH", hh, HH)
                        year = fecha.split("-")[0]
                        month = fecha.split("-")[1]
                        c = cyril["ThetaBn"][i]
                        if (
                            c < 20 or c > 70
                        ):  # los que son muy claramente quasipara o quasiperp
                            l.append(fecha + "/" + hora)

with open("outputs/hora_grupo1.txt", "a") as file:
    for i in l:
        file.write(f"{i}")
        file.write("\n")
