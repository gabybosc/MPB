<<<<<<< HEAD
import numpy as np
import os
import sys
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from socket import gethostname

sys.path.append("..")
from funciones import donde, t_clweb, find_nearest


def find_path(year, month, day, t_i, t_f):
    if gethostname() == "gbosco":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):  # si no existe, usa int(tf)
            path = (
                f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/{t_f}/"
            )
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):
            path = f"../../../datos/clweb/{year}-{month}-{day}/{t_f}/"
            if not os.path.exists(path):
                path = f"../../../datos/clweb/{year}-{month}-{day}/"

    return path


def tiempo_limite(ti, tf):
    # no estoy segura de que usar int(ti)/int(tf) sea buena idea, pero veremos
    if len(str(int(ti))) == 1:
        t_i = "0" + str(int(ti))
    else:
        t_i = int(ti)

    if len(str(int(tf))) == 1:
        t_f = "0" + str(int(tf))
    else:
        t_f = int(tf)
    return (t_i, t_f)


def importar_mag(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    if os.path.isfile(path + "MAG_filtrado.npy"):
        mag = np.load(path + "MAG_filtrado.npy")
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    elif os.path.isfile(path + "mag_filtrado.txt"):
        mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    else:
        mag = np.loadtxt(path + "MAG.asc")
        B = mag[:, 6:9]

    t = t_clweb(mag)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    posicion_cut = mag[inicio:fin, 9:12]
    B_cut = B[inicio:fin]

    if len(B_cut) != len(t_cut):
        print("no tenemos la misma cantidad de datos t que de B")
    return mag, t_cut, B_cut, posicion_cut


###########################################################


def importar_VEX_mag(year, month, day, ti, tf):
    if gethostname() == "gbosco":
        path = "../../../../../media/gabybosc/datos/clweb/"
    else:
        path = "../../../datos/clweb/"

    if os.path.isfile(path + "mag_filtrado.txt"):
        mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    else:
        mag = np.loadtxt(path + "VEX_MAG.asc")
        B = mag[:, 6:9]

    t = t_clweb(mag)  # hdec
    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]

    if len(B_cut) != len(t_cut):
        print("no tenemos la misma cantidad de datos t que de B")
    return mag, t_cut, B_cut


# #########################################################################SWEA


def importar_swea(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    swea = np.loadtxt(path + "SWEA.asc")

    """
    Los datos de SWEA se ven de la siguiente manera: para cada tiempo hay muchos
    valores de JE y de energía. Por eso sólo recorto el tiempo y dejo el JE entero
    y en tal caso en el script principal lo recorto.
    """

    t_swea, idx = np.unique(
        swea[:, 3] + swea[:, 4] / 60 + swea[:, 5] / 3600, return_index=True
    )  # hdec
    # al tomar "unique" estoy borrando los t que se repiten para las diferentes energias

    energias = [50 + i * 50 for i in range(3)]

    inicio = donde(t_swea, ti)
    fin = donde(t_swea, tf)  # este corte no sirve para JE ni para energy!!
    energy = swea[:, 7]
    JE_total = swea[:, -1]
    t_cut = t_swea[inicio:fin]

    JE_cut = np.zeros((len(t_cut), len(energias)))
    for i, energia in enumerate(energias):
        index = np.where(energy == find_nearest(energy, energia))[
            0
        ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
        JE = JE_total[index]
        JE_cut[:, i] = JE[inicio:fin]

    return swea, t_cut, JE_cut


# #########################################################################SWIA


def importar_swia(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    if os.path.isfile(path + "SWICA.asc"):  # si no existe uno llamado SWICA, usa SWIA
        swia = np.loadtxt(path + "SWICA.asc")
    else:
        swia = np.loadtxt(path + "SWIA.asc")

    t = t_clweb(swia)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_cut = swia[inicio:fin, 6]
    vel_cut = swia[inicio:fin, 7:10]
    vel_norm = swia[inicio:fin, -1]

    return swia, t_cut, density_cut, vel_cut, vel_norm


def importar_swicfa(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    swica = np.loadtxt(path + "SWICA.asc")
    swifa = np.loadtxt(path + "SWIFA.asc")

    density_c = swica[:, 6]
    density_f = swifa[:, 6]

    t_c = t_clweb(swica)
    t_f = t_clweb(swifa)

    return t_c, t_f, density_c, density_f


# ######################################################################## SWIA


def importar_vel_swica(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)

    path = find_path(year, month, day, t_i, t_f)

    swia = np.loadtxt(path + "SW_vel.asc")

    t = t_clweb(swia)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    vel_cut = swia[inicio:fin, 6:10]
    vel_norm_cut = swia[inicio:fin, -1]

    return swia, t_cut, vel_cut, vel_norm_cut


# ########################################################################## LPW


def importar_lpw(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)

    path = find_path(year, month, day, t_i, t_f)

    lpw = np.loadtxt(path + "LPW.asc")

    t = t_clweb(lpw)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    e_density_cut = lpw[inicio:fin, 6]
    flag = lpw[inicio:fin, -1]

    return lpw, t_cut, e_density_cut, flag


######################


def importar_static(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    static = np.loadtxt(path + "STATIC.asc")

    t = t_clweb(static)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    mass = static[inicio:fin, 7]
    counts = static[inicio:fin, -1]

    return static, t_cut, mass, counts


# #################################### tiempos t1t2t3t4
def importar_t1t2t3t4(year, month, day, hour):
    fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
        year, month, day, hour
    )
    # hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
    t1 = float(hoja_parametros.cell(fila, 6).value)  # ojo que cuenta desde 1 no desde 0
    t2 = float(hoja_parametros.cell(fila, 7).value)
    t3 = float(hoja_parametros.cell(fila, 8).value)
    t4 = float(hoja_parametros.cell(fila, 9).value)

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


###########
def importar_fila(year, month, day, hora):
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
=======
import numpy as np
import os
import sys
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from socket import gethostname

sys.path.append("..")
from funciones import donde, t_clweb, find_nearest


def find_path(year, month, day, t_i, t_f):
    if gethostname() == "gbosco":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):  # si no existe, usa int(tf)
            path = (
                f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/{t_f}/"
            )
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/{t_i}/"
        if not os.path.exists(path):
            path = f"../../../datos/clweb/{year}-{month}-{day}/{t_f}/"
            if not os.path.exists(path):
                path = f"../../../datos/clweb/{year}-{month}-{day}/"

    return path


def tiempo_limite(ti, tf):
    # no estoy segura de que usar int(ti)/int(tf) sea buena idea, pero veremos
    if len(str(int(ti))) == 1:
        t_i = "0" + str(int(ti))
    else:
        t_i = int(ti)

    if len(str(int(tf))) == 1:
        t_f = "0" + str(int(tf))
    else:
        t_f = int(tf)
    return (t_i, t_f)


def importar_mag(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    if os.path.isfile(path + "MAG_filtrado.npy"):
        mag = np.load(path + "MAG_filtrado.npy")
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    elif os.path.isfile(path + "mag_filtrado.txt"):
        mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    else:
        mag = np.loadtxt(path + "MAG.asc")
        B = mag[:, 6:9]

    t = t_clweb(mag)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    posicion_cut = mag[inicio:fin, 9:12]
    B_cut = B[inicio:fin]

    if len(B_cut) != len(t_cut):
        print("no tenemos la misma cantidad de datos t que de B")
    return mag, t_cut, B_cut, posicion_cut


###########################################################


def importar_VEX_mag(year, month, day, ti, tf):
    if gethostname() == "gbosco":
        path = "../../../../../media/gabybosc/datos/clweb/"
    else:
        path = "../../../datos/clweb/"

    if os.path.isfile(path + "mag_filtrado.txt"):
        mag = np.loadtxt(path + "mag_filtrado.txt", skiprows=2)
        B = mag[:, :3]
        mag = np.loadtxt(path + "MAG.asc")

    else:
        mag = np.loadtxt(path + "VEX_MAG.asc")
        B = mag[:, 6:9]

    t = t_clweb(mag)  # hdec
    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]

    if len(B_cut) != len(t_cut):
        print("no tenemos la misma cantidad de datos t que de B")
    return mag, t_cut, B_cut


# #########################################################################SWEA


def importar_swea(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    swea = np.loadtxt(path + "SWEA.asc")

    """
    Los datos de SWEA se ven de la siguiente manera: para cada tiempo hay muchos
    valores de JE y de energía. Por eso sólo recorto el tiempo y dejo el JE entero
    y en tal caso en el script principal lo recorto.
    """

    t_swea, idx = np.unique(
        swea[:, 3] + swea[:, 4] / 60 + swea[:, 5] / 3600, return_index=True
    )  # hdec
    # al tomar "unique" estoy borrando los t que se repiten para las diferentes energias

    energias = [50 + i * 50 for i in range(3)]

    inicio = donde(t_swea, ti)
    fin = donde(t_swea, tf)  # este corte no sirve para JE ni para energy!!
    energy = swea[:, 7]
    JE_total = swea[:, -1]
    t_cut = t_swea[inicio:fin]

    JE_cut = np.zeros((len(t_cut), len(energias)))
    for i, energia in enumerate(energias):
        index = np.where(energy == find_nearest(energy, energia))[
            0
        ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
        JE = JE_total[index]
        JE_cut[:, i] = JE[inicio:fin]

    return swea, t_cut, JE_cut


# #########################################################################SWIA


def importar_swia(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    if os.path.isfile(path + "SWICA.asc"):  # si no existe uno llamado SWICA, usa SWIA
        swia = np.loadtxt(path + "SWICA.asc")
    else:
        swia = np.loadtxt(path + "SWIA.asc")

    t = t_clweb(swia)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_cut = swia[inicio:fin, 6]
    vel_cut = swia[inicio:fin, 7:10]
    vel_norm = swia[inicio:fin, -1]

    return swia, t_cut, density_cut, vel_cut, vel_norm


def importar_swicfa(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    swica = np.loadtxt(path + "SWICA.asc")
    swifa = np.loadtxt(path + "SWIFA.asc")

    density_c = swica[:, 6]
    density_f = swifa[:, 6]

    t_c = t_clweb(swica)
    t_f = t_clweb(swifa)

    return t_c, t_f, density_c, density_f


# ######################################################################## SWIA


def importar_vel_swica(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)

    path = find_path(year, month, day, t_i, t_f)

    swia = np.loadtxt(path + "SW_vel.asc")

    t = t_clweb(swia)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    vel_cut = swia[inicio:fin, 6:10]
    vel_norm_cut = swia[inicio:fin, -1]

    return swia, t_cut, vel_cut, vel_norm_cut


# ########################################################################## LPW


def importar_lpw(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)

    path = find_path(year, month, day, t_i, t_f)

    lpw = np.loadtxt(path + "LPW.asc")

    t = t_clweb(lpw)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    e_density_cut = lpw[inicio:fin, 6]
    flag = lpw[inicio:fin, -1]

    return lpw, t_cut, e_density_cut, flag


######################


def importar_static(year, month, day, ti, tf):
    t_i, t_f = tiempo_limite(ti, tf)
    path = find_path(year, month, day, t_i, t_f)

    static = np.loadtxt(path + "STATIC.asc")

    t = t_clweb(static)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    mass = static[inicio:fin, 7]
    counts = static[inicio:fin, -1]

    return static, t_cut, mass, counts


# #################################### tiempos t1t2t3t4
def importar_t1t2t3t4(year, month, day, hour):
    fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
        year, month, day, hour
    )
    # hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
    t1 = float(hoja_parametros.cell(fila, 6).value)  # ojo que cuenta desde 1 no desde 0
    t2 = float(hoja_parametros.cell(fila, 7).value)
    t3 = float(hoja_parametros.cell(fila, 8).value)
    t4 = float(hoja_parametros.cell(fila, 9).value)

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


###########
def importar_fila(year, month, day, hora):
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
>>>>>>> 89a36a6263086cb4f9fe8802bd1da03877917733
