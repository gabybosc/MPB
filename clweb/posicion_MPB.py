import numpy as np
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from importar_datos import importar_fila, importar_mag
import sys

sys.path.append("..")
from funciones import fechas, donde

"""
Si no tengo la posición calculada, agarra los t1t2t3t4 y la calcula y la
guarda en el gdocs.
"""


year, month, day, doy = fechas()
hour = int(input("hora\n"))

fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
    year, month, day, doy, hour
)
# hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
t1 = float(hoja_parametros.cell(fila, 6).value)  # ojo que cuenta desde 1 no desde 0
t2 = float(hoja_parametros.cell(fila, 7).value)
t3 = float(hoja_parametros.cell(fila, 8).value)
t4 = float(hoja_parametros.cell(fila, 9).value)

mag, t, B, posicion = importar_mag(year, month, day, t1, t4)

tm = donde(t, np.mean([t1, t2, t3, t4]))
posicion[tm, :]

scope = [
    "https://spreadsheets.google.com/feeds",
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive.file",
    "https://www.googleapis.com/auth/drive",
]

creds = ServiceAccountCredentials.from_json_keyfile_name("../mpb_api.json", scope)

client = gspread.authorize(creds)

hoja_parametros = client.open("MPB").worksheet("Parametros")

hoja_parametros.update_acell(f"AB{fila}", f"{int(posicion[tm, 0])}")
hoja_parametros.update_acell(f"AC{fila}", f"{int(posicion[tm, 1])}")
hoja_parametros.update_acell(f"AD{fila}", f"{int(posicion[tm, 2])}")
