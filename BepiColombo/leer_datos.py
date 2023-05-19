import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

sys.path.append("../")

from funciones import UTC_to_hdec, donde


# path = "../../../../media/gabybosc/datos/vso_ob/vso/"
# path = "../../../datos/vso_ob/vso/"
def importar_bepi(ti, tf):
    path = "../../../datos/vso_ob/"

    # datos = np.loadtxt(path + "mag_der_sc_ob_a001_vso_00000_20210807.tab", dtype=str)
    timedata = np.loadtxt(
        path + "mag_der_sc_ob_a001_vso_00000_20210810.tab", usecols=0, dtype=str
    )
    data = np.loadtxt(
        path + "mag_der_sc_ob_a001_vso_00000_20210810.tab",
        usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10],
        delimiter=",",
    )

    # fecha = timedata[0].split("T")[0]
    hora = [timedata[i].split("T")[1].split("Z")[0] for i in range(len(timedata))]
    hdec = np.array([UTC_to_hdec(hh) for hh in hora])

    x = [data[i][0] for i in range(len(data))]
    y = [data[i][1] for i in range(len(data))]
    z = [data[i][2] for i in range(len(data))]
    pos = np.transpose([x, y, z])

    Bx = [data[i][-6] for i in range(len(data))]
    By = [data[i][-5] for i in range(len(data))]
    Bz = [data[i][-4] for i in range(len(data))]
    B = np.transpose([Bx, By, Bz])

    t_cut = hdec[donde(hdec, ti) : donde(hdec, tf)]
    B_cut = B[donde(hdec, ti) : donde(hdec, tf)]
    pos_cut = pos[donde(hdec, ti) : donde(hdec, tf)]

    return t_cut, B_cut, pos_cut


# plt.plot(t_cut, B_cut)
# plt.show()
# plt.plot(t_cut, np.linalg.norm(B_cut, axis=1))
# plt.show()
# plt.plot(pos_cut[:, 0], np.sqrt(pos_cut[:, 1] ** 2 + pos_cut[:, 2] ** 2))
# plt.show()
