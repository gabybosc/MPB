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

from funciones import donde, next_available_row, doy_to_day

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


def importar_MAG(year, doy, ti, tf):
    """
    Este agarra los datos de pds filtrados y si no existen, busca los de clweb
    """
    year, month, day = doy_to_day(year, doy)

    if gethostname() == "DESKTOP-2GS0QF2":
        os.chdir(f"F:/")
        filt = f"VEX{year}/VEX_mag_filtrado_{year}{doy}.gz"
        magcl = f"VEX{year}/MAG{year}{month}{day}.asc"
        poscl = f"VEX{year}/pos{year}{month}{day}.asc"
    if gethostname() == "gabybosc":
        filt = f"../../../datos/VEX/VEX_mag_filtrado_{year}{doy}.gz"
        magcl = f"../../../datos/VEX/MAG{year}{month}{day}.asc"
        poscl = f"../../../datos/VEX/pos{year}{month}{day}.asc"

    if os.path.exists(filt):
        MAG = np.loadtxt(filt)
        t = MAG[0]
        B = MAG[1:4].T
        pos = MAG[4:].T
        t_pos = 0
        cl = False

    else:
        mag = np.loadtxt(magcl)
        pos = np.loadtxt(poscl)

        hh = mag[:, 3]
        mm = mag[:, 4]
        ss = mag[:, 5]

        B = mag[:, 6:9]
        posicion = pos[:, 6:9]
        t = hh + mm / 60 + ss / 3600  # hdec
        t_pos = pos[:, 3] + pos[:, 4] / 60 + pos[:, 5] / 3600  # hdec
        cl = True

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]

    if type(t_pos) is not int:
        pos_cut = posicion[donde(t_pos, ti): donde(t_pos, tf)]
        tpos_cut = t_pos[donde(t_pos, ti): donde(t_pos, tf)]
    else:
        tpos_cut = t_cut
        pos_cut = pos[inicio:fin]

    if gethostname() == "DESKTOP-2GS0QF2":
        os.chdir("C:/Users/RainbowRider/Documents/GitHub/MPB/VEX/")

    return t_cut, B_cut, pos_cut, cl, tpos_cut  # pos en RV


def importar_MAG_pds(year, doy, ti, tf):
    if gethostname() == "gbosco":
        path = f"../../../../media/gabybosc/datos/VEX/{year}/VEX_MAG_{year}{doy}.tab"
        filt = f"../../../../media/gabybosc/datos/VEX/filtrados/VEX_mag_filtrado_{year}{doy}.gz"
    elif gethostname() == "DESKTOP-2GS0QF2":
        os.chdir(f"G:/")
        path = f"VEX{year}/VEX_MAG_{year}{doy}.tab"
        filt = f"VEX{year}/VEX_mag_filtrado_{year}{doy}.npy"
    else:
        path = f"../../../datos/VEX/{year}/VEX_MAG_{year}{doy}.tab"
        filt = f"../../../datos/VEX/filtrados/VEX_mag_filtrado_{year}{doy}.gz"

    if os.path.exists(path):
        if Path(path).stat().st_size > 1000:
            if os.path.isfile(filt):
                t = np.load(filt)[0]
                B = np.load(filt)[1:]
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


def importar_MAG_clweb(year, doy, ti, tf):
    year, month, day = doy_to_day(year, doy)
    if gethostname() == "DESKTOP-2GS0QF2":
        os.chdir(f"G:/")
        path = f"VEX{year}/"
        mag = np.loadtxt(path + f"MAG{year}{month}{day}.asc")
        pos = np.loadtxt(path + f"pos{year}{month}{day}.asc")[-3:]

    hh = mag[:, 3]
    mm = mag[:, 4]
    ss = mag[:, 5]

    B = mag[:, 6:9]
    t = hh + mm / 60 + ss / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    pos_cut = pos[inicio:fin]

    if gethostname() == "DESKTOP-2GS0QF2":
        os.chdir("C:/Users/RainbowRider/Documents/GitHub/MPB/VEX/")
    return t_cut, B_cut, pos_cut


def importar_ELS_clweb(year, doy, ti, tf):
    year, month, day = doy_to_day(year, doy)
    if gethostname() == "DESKTOP-2GS0QF2":
        os.chdir(f"G:/")
        path = f"VEX{year}/"
        ELS = np.loadtxt(path + f"ELS{year}{month}{day}.asc")
    else:
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

    if gethostname() == "DESKTOP-2GS0QF2":
        os.chdir("C:/Users/RainbowRider/Documents/GitHub/MPB/VEX/")

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
            poner_fecha(hojas, fila, year, str(month), day)
        print("nueva fila")

    return fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste


def importar_t1t2t3t4(year, month, day):
    fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
        year, month, day
    )
    # hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
    t1 = float(hoja_parametros.cell(fila, 6).value)  # ojo que cuenta desde 1 no desde 0
    t2 = float(hoja_parametros.cell(fila, 7).value)
    t3 = float(hoja_parametros.cell(fila, 8).value)
    t4 = float(hoja_parametros.cell(fila, 9).value)

    return t1, t2, t3, t4


def importar_tMVA(year, month, day):
    fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
        year, month, day
    )
    # hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
    ti = float(hoja_MVA.cell(fila, 4).value)  # ojo que cuenta desde 1 no desde 0
    tf = float(hoja_MVA.cell(fila, 5).value)

    return ti, tf
