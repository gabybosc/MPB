import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path = "../../../../media/gabybosc/datos/vso_ob/vso/"

# datos = np.loadtxt(path + "mag_der_sc_ob_a001_vso_00000_20210807.tab", dtype=str)
timedata = np.loadtxt(
    path + "mag_der_sc_ob_a001_vso_00000_20210807.tab", usecols=0, dtype=str
)
data = np.loadtxt(
    path + "mag_der_sc_ob_a001_vso_00000_20210807.tab",
    usecols=[2, 3, 4, 5, 6, 7, 8, 9, 10],
    delimiter=",",
)

"2021-08-07T00:00:06.000000Z, 1/0693014405:00000,   -2335501.76,    -907639.57,       1312.82,     -0.842,     -5.043,     -4.298,   -6.16,   -6.13,   19.47"

fecha = timedata[0].split("T")[0]
hora = timedata[0].split("T")[1].split("Z")[0]

Bx = [data[i][-6] for i in range(len(data))]
By = [data[i][-5] for i in range(len(data))]
Bz = [data[i][-4] for i in range(len(data))]
# B_norm = [data[i][-1] for i in range(len(data))]

plt.plot(Bx)
plt.plot(By)
plt.plot(Bz)
# plt.plot(B_norm)
plt.show()
