import numpy as np
import matplotlib.dates as md
import datetime as dt
import calendar
from scipy.spatial import KDTree

"""
Acá vienen las funciones que son útiles para muchas cosas, no sólo para este análisis.
También están las del análisis de la MPB.
"""


def angulo(v1, v2):
    """Calcula el ángulo (en radianes) entre dos vectores, si no están normalizados los normaliza"""
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)
    angle = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0))
    return angle


def ancho_mpb(t1, t2, t3, t4, normal, vel):
    """Devuelve el ancho minimo y máximo, en kilómetros"""
    deltat_14 = (t4 - t1) * 3600
    deltat_23 = (t3 - t2) * 3600

    v_para = np.dot(vel, normal) * normal  # km/s
    x_14 = np.linalg.norm(v_para * deltat_14)  # en km
    x_23 = np.linalg.norm(v_para * deltat_23)

    return x_14, x_23


def autovectores(M_ij):
    """Devuelve los autovectores del MVA"""
    # ahora quiero los autovectores y autovalores
    [lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

    # Los ordeno de mayor a menor
    idx = lamb.argsort()[::-1]
    lamb = lamb[idx]
    x = x[:, idx]
    # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
    x1 = x[:, 0]
    x2 = x[:, 1]
    x3 = x[:, 2]

    if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
        x3 = -x3
    if any(np.cross(x1, x2) - x3) > 0.01:
        # print("Cambio el signo de x1 para que los av formen terna derecha")
        x1 = -x1

    avec = [x1, x2, x3]
    aval = lamb
    return avec, aval


def deltaB(B):
    """B es un array de Nx3, una matriz."""
    B_medio = np.mean(B, axis=0)
    Bnorm = np.linalg.norm(B_medio)

    prod_interno = np.dot(B - B_medio, B_medio) / Bnorm
    abs_deltaB_para = np.abs(prod_interno) / Bnorm  # |deltaB_para / B|, de Nx1

    N = np.zeros((len(prod_interno), len(B_medio)))

    for i in range(len(N)):
        N[i, :] = prod_interno[i] * B_medio / Bnorm
        deltaB_perp = (B - B_medio) - N
        # y ahora necesito el valor abs de perp
        abs_deltaB_perp = np.abs(deltaB_perp) / Bnorm  # de Nx3
    # ante dudas leer la tesis de lic
    return abs_deltaB_para, abs_deltaB_perp


def Bpara_Bperp(B, t, ti, tf):
    """
    Devuelve Bpara Bperp y el tiempo para plotear
    """
    j_inicial = donde(t, ti)
    j_final = donde(t, tf)

    # Lo hago en ventanas de Mf-Mi, moviendose de a j (1s en baja resolución).
    B_para = np.zeros(j_final - j_inicial)
    B_perp = np.zeros((j_final - j_inicial, 3))
    B_perp_norm = np.zeros(j_final - j_inicial)
    for j in range(j_inicial, j_final - 25):
        Mi = j
        Mf = j + 25
        M_delta = 12  # overlap de 12
        B_delta = B[Mi:Mf]
        deltaB_para, deltaB_perp = deltaB(B_delta)
        B_para[j - j_inicial] = deltaB_para[M_delta]
        B_perp[j - j_inicial, :] = deltaB_perp[M_delta, :]
        B_perp_norm[j - j_inicial] = np.linalg.norm(deltaB_perp[M_delta, :])

        t_plot = t[j_inicial:j_final]

    return B_para, B_perp_norm, t_plot


def corrientes(normal, Bup, Bdown, ancho_mpb):
    """Toma la normal, el campo up/downstream (en nT) y el ancho de la mpb (en km)
    y devuelve j_s y j_v en mA/m y nA/m²"""
    mu = 4 * np.pi * 1e-7  # Henry/m
    js = np.cross(normal, (Bup - Bdown)) / mu  # nA/m
    jv = js / (1000 * ancho_mpb)  # nA/m²

    return js * 1e-6, jv


def day_to_doy(year, month, day):
    date_orbit = dt.date(int(year), int(month), int(day))

    year = date_orbit.strftime("%Y")
    doy = date_orbit.strftime("%j")

    return year, doy


def doy_to_day(year, doy):
    date_orbit = dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    return year, month, day


def diezmar(largo, corto):
    """le doy un array largo que quiero que mida como el corto y me devuelve los
    índices tal que len(largo[idx]) == len(corto). Ambos tienen que tener
    valores similares (e.g. ambos ser arrays temporales)."""
    idx = np.zeros(len(corto))
    for i in range(len(corto)):
        idx[i] = donde(largo, corto[i])
    idx = idx.astype(int)
    return idx


def donde(array, valor):
    """Me dice dónde en un array está el valor más parecido a un valor dado.
    Solo funciona para 1D"""
    resultado = np.where(array == find_nearest(array, valor))[0][0]
    return resultado


def donde_kD(array, valor):
    """Me dice dónde en un array está el valor más parecido a un valor dado.
    Funciona para k dimensiones"""
    kdtree = KDTree(array)
    distancia, indice = kdtree.query(valor)
    resultado = array[indice[0]]
    return resultado


def error(lamb, B, x):
    phi = np.zeros((3, 3))
    M = len(B)
    for i in range(3):
        for j in range(3):
            if i == j:
                phi[i, j] = 0
            else:
                phi[i, j] = np.sqrt(
                    lamb[2]
                    / (M - 1)
                    * (lamb[i] + lamb[j] - lamb[2])
                    / (lamb[i] - lamb[j]) ** 2
                )  # en radianes

    delta_B3 = np.sqrt(
        lamb[2] / (M - 1)
        + (phi[2, 1] * np.dot(np.mean(B, 0), x[1])) ** 2
        + (phi[2, 0] * np.dot(np.mean(B, 0), x[0])) ** 2
    )
    return phi, delta_B3


def find_nearest(array, value):
    """Busca el valor más cercano a uno dado en un array"""
    idx = (np.abs(array - value)).argmin()
    # argmin me da el valor del minimo en el flattened array
    return array[idx]


def find_nearest_inicial(array, value):
    """Busca el valor (V) más cercano a uno dado (value) en un array pidiendo que
    V > value"""
    idx = (np.abs(array - value)).argmin()
    # el indice de la minima dif entre el array y el valor
    if array[idx] < value:
        idx = np.abs((array - (value + 1 / 3600))).argmin()

    return array[idx]


def find_nearest_final(array, value):
    """Busca el valor (V) más cercano a uno dado (value) en un array pidiendo que
    V < value"""
    idx = (np.abs(array - value)).argmin()
    # el indice de la minima dif entre el array y el valor
    if array[idx] > value:
        idx = np.abs((array - (value - 1 / 3600))).argmin()
    return array[idx]


def fechas():
    date_entry = input("Enter a date in YYYY-DDD or YYYY-MM-DD format \n")

    if len(date_entry.split("-")) < 3:
        year, doy = map(int, date_entry.split("-"))
        date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(
            doy - 1
        )  # para convertir el doty en date
    else:
        year, month, day = map(int, date_entry.split("-"))
        date_orbit = dt.date(year, month, day)

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    return year, month, day, doy


def tiempos(string=" "):
    print(string)
    tii = input("Tiempo inicial hh:mm:ss o hdec\n")
    tff = input("Tiempo final hh:mm:ss o hdec\n")
    while tff < tii:
        print("t final no puede ser menor a t inicial. \n")
        tii = input("Tiempo inicial hh:mm:ss o hdec\n")
        tff = input("Tiempo final hh:mm:ss o hdec\n")

    if ":" in tii:
        ti_MVA = UTC_to_hdec(tii)
    else:
        ti_MVA = float(tii)
    if ":" in tff:
        tf_MVA = UTC_to_hdec(tff)
    else:
        tf_MVA = float(tff)

    return ti_MVA, tf_MVA


def Mij(B):
    """Calcula la matriz Mij para un array de Nx3."""
    M_ij = np.zeros((3, 3))
    for i in range(3):  # para las tres coordenadas
        for j in range(3):
            M_ij[i, j] = np.mean(B[:, i] * B[:, j]) - np.mean(B[:, i]) * np.mean(
                B[:, j]
            )
    return M_ij


def next_available_row(sheet):
    """Devuelve la próxima fila vacía en una spreadsheet"""
    str_list = list(filter(None, sheet.col_values(1)))
    return str(len(str_list) + 1)


def proyecciones(B):
    M = len(B)
    Bnorm = np.linalg.norm(B, axis=1)

    n_p = int(M / 2)

    M_ij = Mij(B)

    # ahora quiero los autovectores y autovalores
    [lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

    # Los ordeno de mayor a menor
    idx = lamb.argsort()[::-1]
    lamb = lamb[idx]
    x = x[:, idx]
    # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
    x1 = x[:, 0]
    x2 = x[:, 1]
    x3 = x[:, 2]

    if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
        x3 = -x3
    if any(np.cross(x1, x2) - x3) > 0.01:
        x1 = -x1

    # las proyecciones
    B1 = np.dot(B, x1)
    B2 = np.dot(B, x2)
    B3 = np.dot(B, x3)

    return B1, B2, B3


def rms(x):
    """Calcula el root mean square"""
    sol = np.sqrt(np.vdot(x, x) / x.size)
    return sol


def SZA(posicion, index):
    """El SZA en grados"""
    SZA = angulo(posicion[index, :], [1, 0, 0]) * 180 / np.pi
    return SZA


def unix_to_decimal(t_unix):
    """Le doy un tiempo en unix y me lo pasa a hora decimal."""
    t = np.zeros(np.size(t_unix))
    for i in range(np.size(t)):
        u = dt.datetime.utcfromtimestamp(int(t_unix[i]))
        t[i] = u.second / 3600 + u.minute / 60 + u.hour  # tiempo decimal
    return t


def unix_to_timestamp(t_unix):
    """Le doy un tiempo en unix y me lo pasa a hora UTC."""
    n = len(t_unix)
    timestamps = np.linspace(t_unix[0], t_unix[-1], n)
    dates = [dt.datetime.utcfromtimestamp(ts) for ts in timestamps]  # me lo da en UTC
    datenums = md.date2num(dates)
    return datenums


def UTC_to_hdec(t_UTC):
    """Convierte de UTC a hdec"""
    t = t_UTC.split(":")
    if len(t) == 3:
        t_hdec = int(t[0]) + int(t[1]) / 60 + float(t[2]) / 3600
    elif len(t) == 2:
        t_hdec = int(t[0]) + int(t[1]) / 60

    return t_hdec


def hdec_to_UTC(hdec):
    """Convierte de hdec a utc"""
    h = str(int(hdec)).zfill(2)
    m = str(int((hdec % 1) * 60)).zfill(2)
    s = str(int((((hdec % 1) * 60) % 1) * 60)).zfill(2)
    UTC = f"{h}:{m}:{s}"
    return UTC


def t_clweb(dataset):
    """
    Me convierte el t que está en hh mm ss en diferentes columnas en un t hdec
    """
    hh = dataset[:, 3]
    mm = dataset[:, 4]
    ss = dataset[:, 5]

    t = hh + mm / 60 + ss / 3600  # hdec
    return t


def getrem(ins):
    """this function yields the value behind the decimal point"""
    output = abs(ins - np.fix(ins))
    return output


def datenum(Yr, Mo=1, Da=1, Hr=0, Mi=0, Se=0, Ms=0):
    """this function works as regular datetime.datetime, but allows for float input"""

    # correct faulty zero input
    if Mo < 1:
        Mo += 1
    if Da < 1:
        Da += 1

    # distribute the year fraction over days
    if getrem(Yr) > 0:
        if calendar.isleap(np.floor(Yr)):
            fac = 366
        else:
            fac = 365
        Da = Da + getrem(Yr) * fac
        Yr = int(Yr)
    # if months exceeds 12, pump to years
    while int(Mo) > 12:
        Yr = Yr + 1
        Mo = Mo - 12
    # distribute fractional months to days
    if getrem(Mo) > 0:
        Da = Da + getrem(Mo) * calendar.monthrange(Yr, int(Mo))[1]
        Mo = int(Mo)
    # datetime input for 28 days always works excess is pumped to timedelta
    if Da > 28:
        extraDa = Da - 28
        Da = 28
    else:
        extraDa = 0
    # sometimes input is such that you get 0 day or month values, this fixes this anomaly
    if int(Da) == 0:
        Da += 1
    if int(Mo) == 0:
        Mo += 1

    # datetime calculation
    mytime = dt.datetime(int(Yr), int(Mo), int(Da)) + dt.timedelta(
        days=extraDa + getrem(Da), hours=Hr, minutes=Mi, seconds=Se, microseconds=Ms
    )
    return mytime


def array_datenums(year, month, day, t):
    year = int(year)
    month = int(month)
    day = int(day)
    timestamps = np.array(
        [np.datetime64(datenum(year, month, day, x)) for x in t]
    )  # datenum es una función mía
    return timestamps
