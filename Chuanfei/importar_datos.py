import numpy as np
import pandas as pd
import os
import sys
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from socket import gethostname

sys.path.append("..")
from funciones import donde


def importar_mag(year, month, day, ti, tf):
    if len(str(int(ti))) == 1:
        t_i = "0" + str(int(ti))
    else:
        t_i = int(ti)

    if len(str(int(tf))) == 1:
        t_f = "0" + str(int(tf))
    else:
        t_f = int(tf)

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/{t_i}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):
            path = f"../../../datos/clweb/{year}-{month}-{day}/{t_f}/"

    if os.path.isfile(path + "mag_filtrado.txt"):
        mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    else:
        mag = np.loadtxt(path + "MAG.asc")
        B = mag[:, 6:9]

    hh = mag[:, 3]
    mm = mag[:, 4]
    ss = mag[:, 5]

    t = hh + mm / 60 + ss / 3600  # hdec
    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    posicion_cut = mag[inicio:fin, 9:12]
    B_cut = B[inicio:fin]

    if len(B_cut) != len(t_cut):
        print("no tenemos la misma cantidad de datos t que de B")
    return mag, t_cut, B_cut, posicion_cut


# #########################################################################SWEA


def importar_swea(year, month, day, ti, tf):
    if len(str(int(ti))) == 1:
        t_i = "0" + str(int(ti))
    else:
        t_i = int(ti)

    if len(str(int(tf))) == 1:
        t_f = "0" + str(int(tf))
    else:
        t_f = int(tf)
    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):
            path = f"../../../datos/clweb/{year}-{month}-{day}/{t_f}/"

    swea = np.loadtxt(path + "SWEA.asc")

    t_swea = np.unique(swea[:, 3] + swea[:, 4] / 60 + swea[:, 5] / 3600)  # hdec

    energias = [50 + i * 50 for i in range(3)]

    energy = swea[:, 7]
    JE_total = swea[:, -1]

    inicio = donde(t_swea, ti)
    fin = donde(t_swea, tf)

    t_cut = t_swea[inicio:fin]

    idx = []
    JE = []
    JE_cut = []
    for energia in energias:
        idx.append(donde(energy, energia))
    for i in range(len(idx)):
        JE.append(JE_total[idx[i]])
        JE_cut.append(JE[i][inicio:fin])

    return (swea, t_cut, energias, JE_cut)


# ######################################################################## SWIA


def importar_swica(year, month, day, ti, tf):
    if len(str(int(ti))) == 1:
        t_i = "0" + str(int(ti))
    else:
        t_i = int(ti)

    if len(str(int(tf))) == 1:
        t_f = "0" + str(int(tf))
    else:
        t_f = int(tf)

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):
            path = f"../../../datos/clweb/{year}-{month}-{day}/{t_f}/"

    swia = np.loadtxt(path + "SWICA.asc")

    t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_cut = swia[inicio:fin, 6]
    vel_cut = swia[inicio:fin, 7:10]

    return swia, t_cut, density_cut, vel_cut


def importar_swifa(year, month, day, ti, tf):

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    swia = np.loadtxt(path + "SWIFA.asc")

    t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_cut = swia[inicio:fin, 6]
    vel_cut = swia[inicio:fin, 7:10]

    return swia, t_cut, density_cut, vel_cut


# ########################################################################## LPW


def importar_lpw(year, month, day, ti, tf):
    if len(str(int(ti))) == 1:
        t_i = "0" + str(int(ti))
    else:
        t_i = int(ti)

    if len(str(int(tf))) == 1:
        t_f = "0" + str(int(tf))
    else:
        t_f = int(tf)

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):
            path = f"../../../datos/clweb/{year}-{month}-{day}/{t_f}/"

    lpw = np.loadtxt(path + "LPW.asc")

    t = lpw[:, 3] + lpw[:, 4] / 60 + lpw[:, 5] / 3600

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    e_density_cut = lpw[inicio:fin, 6]
    flag = lpw[inicio:fin, -1]

    return lpw, t_cut, e_density_cut, flag


# ################################################3 STATIC (solo el de Chris F.)


def importar_STATIC(year, month, day, ti, tf):

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = "../../../datos/STATIC/"

    STATIC = pd.read_csv("STATIC_Ni_2016-03-16.txt")

    t_hdec = []
    for t_UTC in STATIC["Time [UTC]"]:
        (dia, hora) = t_UTC.split("/")
        (h, m, s) = hora.split(":")
        t_hdec.append(int(h) + int(m) / 60 + int(s) / 3600)

    t = np.array(t_hdec)
    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    H_density = STATIC["Density H+ [/cc]"][inicio:fin]
    O_density = STATIC["Density O+ [/cc]"][inicio:fin]
    O2_density = STATIC["Density O2+ [/cc]"][inicio:fin]
    CO2_density = STATIC["Density CO2+ [/cc]"][inicio:fin]

    return STATIC, t_cut, H_density, O_density, O2_density, CO2_density


# #################################### tiempos t1t2t3t4
def importar_t1t2t3t4(year, doy, hour):
    fila = input("Qué fila es en el archivo MPB?\n")
    hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
    # ojo que cuenta desde 1 no desde 0
    t1 = float(hoja_parametros.cell(fila, 6).value)
    t2 = float(hoja_parametros.cell(fila, 7).value)
    t3 = float(hoja_parametros.cell(fila, 8).value)
    t4 = float(hoja_parametros.cell(fila, 9).value)

    """si quiero importarlo del .txt hay que descomentar lo siguiente"""
    # datos = np.loadtxt("../outputs/t1t2t3t4.txt")
    # fecha = datos[:, 0:2].astype(int)
    # hora = datos[:, 2].astype(int)
    # fecha_in = [int(year), int(doy)]
    # for j in range(len(datos)):
    #     if all(fecha[j, 0:2]) == all(fecha_in):
    #         if hora[j] == hour or hora[j] == hour + 1 or hora[j] == hour - 1:
    #             idx = j
    # t1, t2, t3, t4 = datos[idx, 2:]

    return t1, t2, t3, t4


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


##########
def importar_fila(year, month, day, doy, hora):
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
    hh = hoja_parametros.col_values(2)[3:]

    fila = None

    for i in range(len(fecha)):
        dd, mm, yy = fecha[i].split()
        mes = meses[mm]
        if (
            yy == str(year)
            and mes == month
            and int(dd) == int(day)
            and int(hh[i]) == int(hora)
        ):
            fila = i + 4
            break

    if fila is None:
        print("no encuentro la fila")

    return fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste
