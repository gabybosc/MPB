import numpy as np
import datetime as dt
from funciones import donde, unix_to_decimal
import cdflib as cdf
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from socket import gethostname


def importar_mag_1s(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    if gethostname() == "magneto2":
        path = f"../../../../media/gabybosc/datos/MAG_1s/{year}/"
    elif gethostname() == "gabybosc":
        path = f"../../datos/MAG_1s/{year}/"
    else:
        path = f"../../../datos/MAG_1s/{year}-{month}-{day}/"

    mag = np.loadtxt(
        path + f"mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts",
        skiprows=160,
    )

    hh = mag[:, 2]
    mm = mag[:, 3]
    ss = mag[:, 4]

    t = hh + mm / 60 + ss / 3600  # hdec

    B = mag[:, 7:10]

    posicion = mag[:, 11:14]

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    posicion_cut = posicion[inicio:fin]
    return mag, t_cut, B_cut, posicion_cut


########################################################################


def importar_mag(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    if gethostname() == "magneto2":
        path = f"../../../../media/gabybosc/datos/MAG_hires/"
    elif gethostname() == "gabybosc":
        path = "../../datos/MAG_hires/"
    else:
        path = f"../../../datos/MAG_hires/"

    mag = np.loadtxt(
        path + f"mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts", skiprows=160
    )

    hh = mag[:, 2]
    mm = mag[:, 3]
    ss = mag[:, 4]

    t = hh + mm / 60 + ss / 3600  # hdec

    B = mag[:, 7:10]

    posicion = mag[:, 11:14]

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    posicion_cut = posicion[inicio:fin]
    mag_cut = mag[inicio:fin]

    return mag_cut, t_cut, B_cut, posicion_cut


# ######################################################################## SWEA
def importar_swea(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    if gethostname() == "magneto2":
        path = f"../../../../media/gabybosc/datos/SWEA/"
    elif gethostname() == "gabybosc":
        path = "../../datos/SWEA/"
    else:
        path = f"../../../datos/SWEA/"

    swea = cdf.CDF(path + f"/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf")

    flux_all = swea.varget("diff_en_fluxes")
    energia = swea.varget("energy")
    t_unix = swea.varget("time_unix")

    t = unix_to_decimal(t_unix)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]

    flux = flux_all[inicio:fin]
    flux_plot = np.transpose(flux)[::-1]

    return swea, t_cut, energia, flux_plot


# ######################################################################## SWIA


def importar_swia(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    if gethostname() == "magneto2":
        path = f"../../../../media/gabybosc/datos/SWIA/"
    elif gethostname() == "gabybosc":
        path = "../../datos/SWIA/"
    else:
        path = f"../../../datos/SWIA/"

    swia = cdf.CDF(path + f"mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf")

    t_unix = swia.varget("time_unix")
    density = swia.varget("density")  # cm⁻³
    temperature = swia.varget("temperature_mso")  # eV
    vel_mso_xyz = swia.varget("velocity_mso")  # km/s

    t_swia = unix_to_decimal(t_unix)
    inicio = donde(t_swia, ti)
    fin = donde(t_swia, tf)

    t_cut = t_swia[inicio:fin]
    density_cut = density[inicio:fin]
    temperature_cut = temperature[inicio:fin]
    vel_mso_cut = vel_mso_xyz[inicio:fin]  # km/s

    return swia, t_cut, density_cut, temperature_cut, vel_mso_cut


# def importar_swicfa(year, month, day, ti, tf):
#     # esto tiene el problema de que no tengo las densidades y eso en swica directo...
#     date_orbit = dt.date(int(year), int(month), int(day))
#     year = date_orbit.strftime("%Y")
#     month = date_orbit.strftime("%m")
#     day = date_orbit.strftime("%d")
#
#     if gethostname() == "magneto2":
#         path = f"../../../../media/gabybosc/datos/SWIA/"
#     elif gethostname() == "gabybosc":
#         path = "../../datos/SWIA/"
#     else:
#         path = f"../../../datos/SWIA/"
#
#     swica = cdf.CDF(path + f"mvn_swi_l2_coarsearc3d_{year}{month}{day}_v01_r01.cdf")
#     swifa = cdf.CDF(path + f"mvn_swi_l2_finearc3d_{year}{month}{day}_v01_r01.cdf")
#
#     t_unix = swia.varget("time_unix")
#     density = swia.varget("density")  # cm⁻³
#     temperature = swia.varget("temperature_mso")  # eV
#     vel_mso_xyz = swia.varget("velocity_mso")  # km/s
#
#     t_swia = unix_to_decimal(t_unix)
#     inicio = donde(t_swia, ti)
#     fin = donde(t_swia, tf)
#
#     t_cut = t_swia[inicio:fin]
#     density_cut = density[inicio:fin]
#     temperature_cut = temperature[inicio:fin]
#     vel_mso_cut = vel_mso_xyz[inicio:fin]  # km/s
#
#     return swia, t_cut, density_cut, temperature_cut, vel_mso_cut

# ###################################################################### LPW
def importar_lpw(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    if gethostname() == "magneto2":
        path = f"../../../../media/gabybosc/datos/LPW/"
    elif gethostname() == "gabybosc":
        path = "../../datos/LPW/"
    else:
        path = f"../../../datos/LPW/"

    lpw = cdf.CDF(path + f"mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf")

    t_unix = lpw.varget("time_unix")
    e_density = lpw.varget("data")[:, 3]

    t = unix_to_decimal(t_unix)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    e_density_cut = e_density[inicio:fin]

    return lpw, t_cut, e_density_cut


###########################
def importar_static(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    if gethostname() == "magneto2":
        path = f"../../../../media/gabybosc/datos/STATIC/"
    elif gethostname() == "gabybosc":
        path = "../../datos/STATIC/"
    else:
        path = f"../../../datos/STATIC/"

    static = cdf.CDF(path + f"mvn_sta_l2_c6-32e64m_{year}{month}{day}_v02_r01.cdf")

    t_unix = static.varget("time_unix")
    mass = static.varget("mass_arr")
    energy = static.varget("energy")

    t = unix_to_decimal(t_unix)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t[inicio:fin]
    mass_cut = mass[inicio:fin]
    energy_cut = energy[inicio:fin]

    return static, t_cut, mass_cut, energy_cut


###########################
def importar_t1t2t3t4(year, month, day, doy, hour):
    fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
        year, month, day, doy, hour
    )
    # hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_gdocs()
    t1 = float(hoja_parametros.cell(fila, 6).value)  # ojo que cuenta desde 1 no desde 0
    t2 = float(hoja_parametros.cell(fila, 7).value)
    t3 = float(hoja_parametros.cell(fila, 8).value)
    t4 = float(hoja_parametros.cell(fila, 9).value)

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


###################################
def importar_gdocs():
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive.file",
        "https://www.googleapis.com/auth/drive",
    ]

    creds = ServiceAccountCredentials.from_json_keyfile_name("mpb_api.json", scope)

    client = gspread.authorize(creds)

    hoja_parametros = client.open("MPB").worksheet("Parametros")
    hoja_MVA = client.open("MPB").worksheet("MVA")
    hoja_Bootstrap = client.open("MPB").worksheet("Bootstrap")
    hoja_Ajuste = client.open("MPB").worksheet("Ajuste")

    return hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste


###########
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
