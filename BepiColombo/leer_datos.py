import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

sys.path.append("../")
from funciones import UTC_to_hdec, donde

# path = "../../../../media/gabybosc/datos/vso_ob/vso/"
path = "../../../datos/vso_ob/vso/"

# datos = np.loadtxt(path + "mag_der_sc_ob_a001_vso_00000_20210807.tab", dtype=str)
timedata = np.loadtxt(
    path + "mag_der_sc_ob_a001_vso_00000_20210810.tab", usecols=0, dtype=str
)
data = np.loadtxt(
    path + "mag_der_sc_ob_a001_vso_00000_20210810.tab",
    usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10],
    delimiter=",",
)

"2021-08-07T00:00:06.000000Z, 1/0693014405:00000,   -2335501.76,    -907639.57,       1312.82,     -0.842,     -5.043,     -4.298,   -6.16,   -6.13,   19.47"

fecha = timedata[0].split("T")[0]
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

t_cut = hdec[donde(hdec, 13) : donde(hdec, 15)]
B_cut = B[donde(hdec, 13) : donde(hdec, 15)]
pos_cut = pos[donde(hdec, 13) : donde(hdec, 15)]

plt.plot(t_cut, B_cut)
plt.plot(t_cut, np.linalg.norm(B_cut, axis=1))
plt.plot(pos_cut[:, 0], np.sqrt(pos_cut[:, 1] ** 2 + pos_cut[:, 2] ** 2))
plt.show()
