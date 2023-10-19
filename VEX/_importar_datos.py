import numpy as np
import sys
from pathlib import Path
import os as os
from os.path import exists
from glob import glob
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from socket import gethostname


sys.path.append("..")

from funciones import donde, next_available_row

np.set_printoptions(precision=4)


def importar_VEX_mag_AMDA(year, month, day, ti, tf):
    path_B = glob(f"../../../datos/VEX/{day}{month}{year}/*.txt")
    path_pos = glob(f"../../../datos/VEX/{day}{month}{year}_pos/*.txt")
    B = np.genfromtxt(path_B[0], usecols=[1, 2, 3])
    pos = np.genfromtxt(path_pos[0], usecols=[1, 2, 3])
    tt = np.genfromtxt(path_B[0], usecols=0, dtype="str")
    fecha = np.array([x.split("T") for x in tt])
    hora = np.array([x.split(":") for x in fecha[:, 1]])
    hh = np.array([int(x) for x in hora[:, 0]])
    mm = np.array([int(x) for x in hora[:, 1]])
    ss = np.array([float(x) for x in hora[:, 2]])
    t = hh + mm / 60 + ss / 3600  # hdec
    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    return t_cut, B_cut, pos


def importar_MAG_pds(year, doy, ti, tf):
    if gethostname() == "gbosco":
        path = f"../../../../media/gabybosc/datos/VEX/{year}/VEX_MAG_{year}{doy}.tab"
        filt = f"../../../../media/gabybosc/datos/VEX/filtrados/VEX_mag_filtrado_{year}{doy}.gz"
    elif gethostname() == "DESKTOP-2GS0QF2":
        os.chdir(f"G:/")
        path = f"VEX{year}/VEX_MAG_{year}{doy}.tab"
        filt = f"VEX{year}/VEX_mag_filtrado_{year}{doy}.gz"
    else:
        path = f"../../../datos/VEX/{year}/VEX_MAG_{year}{doy}.tab"
        filt = f"../../../datos/VEX/filtrados/VEX_mag_filtrado_{year}{doy}.gz"

    if os.path.exists(path):
        if Path(path).stat().st_size > 1000:
            if os.path.isfile(filt):
                B = np.loadtxt(filt)
                # B = np.genfromtxt(path, skip_header=1, usecols=[5, 6, 7])
            else:
                B = np.genfromtxt(path, skip_header=1, usecols=[5, 6, 7])
            pos = np.genfromtxt(path, skip_header=1, usecols=[8, 9, 10])
            tt = np.genfromtxt(path, skip_header=1, usecols=0, dtype="str")

            hora = np.array([x.split("T")[1] for x in tt])
            t = np.array(
                [
                    int(x.split(":")[0])
                    + int(x.split(":")[1]) / 60
                    + float(x.split(":")[2]) / 3600
                    for x in hora
                ]
            )  # hdec

            inicio = donde(t, ti)
            fin = donde(t, tf)

            t_cut = t[inicio:fin]
            B_cut = B[inicio:fin]
            pos_cut = pos[inicio:fin]
        if (
            gethostname() == "DESKTOP-2GS0QF2"
        ):  # si estoy en la pc tengo que volver al dir original
            os.chdir("C:/Users/RainbowRider/Documents/GitHub/MPB/VEX/")
        else:
            t_cut, B_cut, pos_cut = 0, 0, 0
    return t_cut, B_cut, pos_cut


def importar_MAG_clweb(year, month, day, ti, tf):
    mag = np.loadtxt(f"../../../datos/clweb/VEX_MAG_{year}{month}{day}.asc")

    hh = mag[:, 3]
    mm = mag[:, 4]
    ss = mag[:, 5]

    B = mag[:, 6:9]
    t = hh + mm / 60 + ss / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    return t_cut, B_cut


def importar_ELS_clweb(year, month, day, ti, tf):
    ELS = np.loadtxt(f"../../../datos/clweb/VEX_ASPERA4_ELS_{year}{month}{day}.asc")
    hh = ELS[:, 3]
    mm = ELS[:, 4]
    ss = ELS[:, 5]

    """
    Para cada tiempo hay muchos valores de JE y de energía. Por eso sólo recorto
    el tiempo y dejo el JE entero y en tal caso en el script principal lo recorto.
    """

    t = hh + mm / 60 + ss / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    ELS_cut = ELS[inicio:fin, :]

    return t_cut, ELS_cut


# ################################## godcs
def importar_gdocs():
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive.file",
        "https://www.googleapis.com/auth/drive",
    ]

    creds = ServiceAccountCredentials.from_json_keyfile_name("../mpb_api.json", scope)

    client = gspread.authorize(creds)

    hoja_parametros = client.open("MPB").worksheet("Parametros")
    hoja_MVA = client.open("MPB").worksheet("MVA")
    hoja_Bootstrap = client.open("MPB").worksheet("Bootstrap")
    hoja_Ajuste = client.open("MPB").worksheet("Ajuste")

    return hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste


def poner_fecha(hoja, nr, year, month, day):
    meses = {
        "01": "Jan",
        "02": "Feb",
        "03": "Mar",
        "04": "Apr",
        "05": "May",
        "06": "Jun",
        "07": "Jul",
        "08": "Aug",
        "09": "Sep",
        "10": "Oct",
        "11": "Nov",
        "12": "Dec",
    }
    hoja.update_acell(f"A{nr}", f"{day} {meses[month]} {year}")
    hoja.update_acell(f"C{nr}", "VEX")


def importar_fila(year, month, day):
    hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()

    meses = {
        "Jan": "01",
        "Feb": "02",
        "Mar": "03",
        "Apr": "04",
        "May": "05",
        "Jun": "06",
        "Jul": "07",
        "Aug": "08",
        "Sep": "09",
        "Oct": "10",
        "Nov": "11",
        "Dec": "12",
    }
    fecha = hoja_parametros.col_values(1)[3:]

    fila = None

    for i in range(len(fecha)):
        dd, mm, yy = fecha[i].split(" ")
        mes = meses[mm]
        if yy == str(year) and mes == month and int(dd) == int(day):
            fila = i + 4
            break

    if fila is None:
        fila = next_available_row(hoja_parametros)
        for hojas in [hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste]:
            poner_fecha(hojas, fila, year, month, day)
        print("nueva fila")

    return fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste
