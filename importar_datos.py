import numpy as np
import datetime as dt
from funciones import find_nearest, unix_to_decimal
import cdflib as cdf


def importar_mag_1s(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    path = f"../../../../media/gabybosc/datos/MAG_1s/{year}/"

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

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    posicion_cut = posicion[inicio:fin]
    return mag, t_cut, B_cut, posicion_cut


###############################################################################################


def importar_mag(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    path = f"../../../../media/gabybosc/datos/MAG_hires/"

    mag = np.loadtxt(
        path + f"mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts", skiprows=160
    )

    hh = mag[:, 2]
    mm = mag[:, 3]
    ss = mag[:, 4]

    t = hh + mm / 60 + ss / 3600  # hdec

    B = mag[:, 7:10]

    posicion = mag[:, 11:14]

    inicio = np.where(t == find_nearest(t, ti))[0][0]
    fin = np.where(t == find_nearest(t, tf))[0][0]

    t_cut = t[inicio:fin]
    B_cut = B[inicio:fin]
    posicion_cut = posicion[inicio:fin]
    return mag, t_cut, B_cut, posicion_cut


# ##############################################################################################SWEA
# def importar_swea(year, month, day, ti, tf):
#     date_orbit = dt.date(int(year), int(month), int(day))
#     year = date_orbit.strftime("%Y")
#     month = date_orbit.strftime("%m")
#     day = date_orbit.strftime("%d")
#     doy = date_orbit.strftime("%j")
#
#     path = f'../../../../media/gabybosc/datos/SWEA/'
#
#     swea = cdf.CDF(path + f'SWEA/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf.cdf')
#
#     energy = swea[:, 7]
#     JE_total = swea[:, -1]
#
#     t = np.unique(swea[:,3] + swea[:,4]/60 + swea[:,5]/3600) #hdec
#
#     energias = [50 + i*50 for i in range(3)]
#
#     inicio = np.where(t == find_nearest(t, ti))[0][0]
#     fin = np.where(t == find_nearest(t, tf))[0][0]
#
#     return(swea, t, energias)
# ##############################################################################################SWIA


def importar_swia(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    path = f"../../../../media/gabybosc/datos/SWIA/"

    swia = cdf.CDF(path + f"mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf")

    t_unix = swia.varget("time_unix")
    density = swia.varget("density")  # cm⁻³
    temperature = swia.varget("temperature_mso")  # eV
    vel_mso_xyz = swia.varget("velocity_mso")  # km/s

    t_swia = unix_to_decimal(t_unix)
    inicio = np.where(t_swia == find_nearest(t_swia, ti))[0][0]
    fin = np.where(t_swia == find_nearest(t_swia, tf))[0][0]

    t_cut = t_swia[inicio:fin]
    density_cut = density[inicio:fin]
    temperature_cut = temperature[inicio:fin]
    vel_mso_cut = vel_mso_xyz[inicio:fin]  # km/s

    return swia, t_cut, density_cut, temperature_cut, vel_mso_cut


# ############################################################################################## LPW
# def importar_lpw(year, month, day, ti, tf):
#     date_orbit = dt.date(int(year), int(month), int(day))
#     year = date_orbit.strftime("%Y")
#     month = date_orbit.strftime("%m")
#     day = date_orbit.strftime("%d")
#     doy = date_orbit.strftime("%j")
#
#     path = f'../../../../media/gabybosc/datos/LPW/'
#
#     swia = cdf.CDF(path + f'mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf')
#
#     e_density = lpw[:,-1]
#
#     t = lpw[:,3] + lpw[:,4]/60 + lpw[:,5]/3600
#
#     inicio = np.where(t == find_nearest(t, ti))[0][0]
#     fin = np.where(t == find_nearest(t, tf))[0][0]
#
#     t_cut = t[inicio:fin]
#     e_density_cut = e_density[inicio:fin]
#
#     return(lpw, t_cut, e_density_cut)
