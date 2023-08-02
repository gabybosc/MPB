import numpy as np
import csv as csv
import pandas as pd
import sys

sys.path.append("..")
from funciones import donde

# path = "../../Documents/"
jacob = pd.read_csv("CatalogoMAVEN_Jacob.csv")
# jacob.rename(
#     columns={"yyyy-mm-dd HH:MM:SS": "date", "ThetaBn (deg)": "ThetaBn"}, inplace=True
# )

date = [fecha[:10] for fecha in jacob["TIME"]]
time = [fecha[11:] for fecha in jacob["TIME"]]
x = jacob["XMSO [km]"]
y = jacob["YMSO [km]"]
z = jacob["ZMSO [km]"]
theta = jacob["thetaBN [deg]"]
beta = jacob["beta"]
# pdyn = jacob["Pdyn"]  # está vacía


"""
mira el catalogo de Jacob
quiero elegir 10 fechas con Pdyn y Theta variados
La lista de pdyn está vacía así que uso beta
primero filtro para quedarme con los que sean quasipara o quasiperp
grupo1 = 2014/11 2014/12 2015/01 2015/2 2015/3 2022/4 2022/5 2022/06 2022/7 2022/8
grupo2 = 2018/7 2018/8 2018/09 2018/10 2018/11 2018/12 2020/6 2020/07 2020/8 2020/9 2020/10 2020/11
grupo3 = 2015/9 2015/10 2015/11 2015/12 2016/1
grupo4 = 2017/7 2017/8  2017/09 2017/10 2017/11 2017/12 2019/06 2019/7 2019/8 2019/09 2019/10 2019/11 2021/05 2021/6 2021/7 2021/08 2021/9
"""

# # for mes in ["01", "02", "03","04", "05", "06", "07", "08", "09", "10", "11", "12"]:
# grupo4 = []
# anio = "2017"
# for mes in ["07", "08", "09", "10", "11", "12"]:
#     for i, val in enumerate(date):
#         year = val.split("-")[0]
#         month = val.split("-")[1]
#         if year == anio and month == mes:
#             grupo4.append(i)
# anio = "2019"
# for mes in ["06", "07", "08", "09", "10", "11"]:
#     for i, val in enumerate(date):
#         year = val.split("-")[0]
#         month = val.split("-")[1]
#         if year == anio and month == mes:
#             grupo4.append(i)
# anio = "2021"
# for mes in ["05", "06", "07", "08", "09"]:
#     for i, val in enumerate(date):
#         year = val.split("-")[0]
#         month = val.split("-")[1]
#         if year == anio and month == mes:
#             grupo4.append(i)

# np.save("grupo1.npy", grupo1)
# np.save("grupo2.npy", grupo2)
# np.save("grupo3.npy", grupo3)
# np.save("grupo4.npy", grupo4)

"""
lo mejor va a ser hacer cuatro listas, una de para, una de perp, 
una de thetas unicos y una de beta unicos y después comparar las listas
"""


for n in [1, 2, 3, 4]:
    grupo = np.load(
        f"outs_catalogo_previa/grupo{n}.npy"
    )  # idx del grupo en la data de jacob
    print("len grupo ", len(grupo))

    """
    Hace dos listas, con los índices de los cruces paralelos y perpendiculares
    """
    l_para = []
    l_perp = []
    for i in grupo:
        c = theta[i]
        if c < 30:
            l_para.append(i)
        elif c > 70:
            l_perp.append(i)

    print("len para ", len(l_para))
    print("len perp ", len(l_perp))

    """
    lista de theta unicos
    """

    th_unique, idx_th = np.unique(theta[grupo], return_index=True)
    print("len theta unicos ", len(idx_th))

    beta_unique, idx_beta = np.unique(beta[grupo], return_index=True)
    print("len pdyn unicos ", len(idx_beta))

    if len(idx_th) and len(idx_beta) == len(grupo):
        i = 0
        for l in [l_para, l_perp]:
            posicion = np.transpose([x[l], y[l], z[l]])
            if i == 0:
                g = "para"
            elif i == 1:
                g = "perp"
            fecha = np.array([date[j] for j in l])
            hora = np.array([time[j] for j in l])
            np.save(f"outs_catalogo_previa/grupo{n}/fecha_{g}.npy", fecha)
            np.save(f"outs_catalogo_previa/grupo{n}/hora_{g}.npy", hora)
            np.save(f"outs_catalogo_previa/grupo{n}/beta_{g}.npy", beta[l])
            np.save(f"outs_catalogo_previa/grupo{n}/theta_{g}.npy", theta[l])
            np.save(f"outs_catalogo_previa/grupo{n}/posicion_{g}.npy", posicion)
            i += 1
