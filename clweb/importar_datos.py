import numpy as np
import os
import sys
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from socket import gethostname

sys.path.append("..")
from funciones import donde


def importar_mag(year, month, day, ti, tf):

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

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
    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    swea = np.loadtxt(path + "SWEA.asc")

    """
    Los datos de SWEA se ven de la siguiente manera: para cada tiempo hay muchos
    valores de JE y de energía. Por eso sólo recorto el tiempo y dejo el JE entero
    y en tal caso en el script principal lo recorto.
    """

    t_swea = np.unique(swea[:, 3] + swea[:, 4] / 60 + swea[:, 5] / 3600)  # hdec
    # al tomar "unique" estoy borrando los t que se repiten para las diferentes energias

    energias = [50 + i * 50 for i in range(3)]

    inicio = donde(t_swea, ti)
    fin = donde(t_swea, tf)  # este corte no sirve para JE ni para energy!!

    t_cut = t_swea[inicio:fin]

    energy = swea[:, 7]
    JE_total = swea[:, -1]
    #
    # for energia in energias:
    #     index = donde(energy, energia)
    #     JE = JE_total[index]
    #     plt.semilogy(
    #         t_swea[inicio_swea:fin_swea],
    #         JE[inicio_swea:fin_swea],
    #         label=f"{energia} eV",
    #     )
    #

    return swea, t_cut, energias


# #########################################################################SWIA


def importar_swia(year, month, day, ti, tf):
    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    swia = np.loadtxt(path + "SWIA.asc")

    t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_cut = swia[inicio:fin, -1]

    return swia, t_cut, density_cut


def importar_swicfa(year, month, day, ti, tf):
    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    swica = np.loadtxt(path + "SWICA.asc")
    swifa = np.loadtxt(path + "SWIFA.asc")

    density_c = swica[:, -1]
    density_f = swifa[:, -1]

    t_c = swica[:, 3] + swica[:, 4] / 60 + swica[:, 5] / 3600  # hdec
    t_f = swifa[:, 3] + swifa[:, 4] / 60 + swifa[:, 5] / 3600  # hdec

    return t_c, t_f, density_c, density_f


# ######################################################################## SWIA


def importar_swia_vel(year, month, day, ti, tf):

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    swia = np.loadtxt(path + "SWIA_vel.asc")

    t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    density_cut = swia[inicio:fin, 6]
    vel_cut = swia[inicio:fin, 7:10]

    return swia, t_cut, density_cut, vel_cut


def importar_vel_swica(year, month, day, ti, tf):

    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    swia = np.loadtxt(path + "SW_vel.asc")

    t = swia[:, 3] + swia[:, 4] / 60 + swia[:, 5] / 3600  # hdec

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    vel_cut = swia[inicio:fin, 6:10]
    vel_norm_cut = swia[inicio:fin, -1]

    return swia, t_cut, vel_cut, vel_norm_cut


# ########################################################################## LPW


def importar_lpw(year, month, day, ti, tf):
    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    lpw = np.loadtxt(path + "LPW.asc")

    t = lpw[:, 3] + lpw[:, 4] / 60 + lpw[:, 5] / 3600

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    e_density_cut = lpw[inicio:fin, -1]

    return lpw, t_cut, e_density_cut


######################
def importar_static(year, month, day, ti, tf):
    if gethostname() == "magneto2":
        path = f"../../../../../media/gabybosc/datos/clweb/{year}-{month}-{day}/"
    else:
        path = f"../../../datos/clweb/{year}-{month}-{day}/"

    static = np.loadtxt(path + "STATIC.asc")

    t = static[:, 3] + static[:, 4] / 60 + static[:, 5] / 3600

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    mass = static[inicio:fin, 7]
    counts = static[inicio:fin, -1]

    return static, t_cut, mass, counts


# #################################### tiempos t1t2t3t4
def importar_t1t2t3t4(year, doy, hour):
    fila = input("Qué fila es en el archivo MPB?\n")
    hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
    t1 = float(hoja_parametros.cell(fila, 6).value)  # ojo que cuenta desde 1 no desde 0
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
