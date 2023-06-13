import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")
from funciones import donde, t_clweb

path = "../../../datos/clweb/2015-01/"

mag = np.loadtxt(path + "MAG.asc")
swica = np.loadtxt(path + "SWICA.asc")
swifa = np.loadtxt(path + "SWIFA.asc")
# swea = np.loadtxt(path + "SWEA.asc")

dic = {"mag": mag, "swica": swica, "swifa": swifa}

for d in dic:
    data = dic[str(d)]
    year = data[:, 0]
    month = data[:, 1]
    day = data[:, 2]
    idx = [
        0,
    ]  # para que siempre cuente al día 1
    for i in range(1, len(day)):
        if day[i] != day[i - 1]:
            idx.append(i)

    for j in range(len(idx) - 1):
        inicio = idx[j]
        fin = idx[j + 1] - 1
        np.save(path + f"{d}_dia{int(day[inicio])}", data[inicio:fin])


"""
mag = np.loadtxt(path + "MAG.asc")

year = mag[:, 0]
month = mag[:, 1]
day = mag[:, 2]
B = mag[:, 6:9]
t = t_clweb(mag)

plt.plot(t, B)
plt.show()

idx = [
    0,
]  # para que siempre cuente al día 1
for i in range(1, len(day)):
    if day[i] != day[i - 1]:
        idx.append(i)

for j in range(len(idx) - 1):
    inicio = idx[j]
    fin = idx[j + 1] - 1
    np.save(path + f"dia{int(day[inicio])}", mag[inicio:fin])
    # # si quiero chequear que el corte sea correcto
    # B = mag[inicio:fin, 6:9]
    # t = t_clweb(mag[inicio:fin])

    # plt.plot(t, B)
    # plt.show()

"""
