import numpy as np
import glob as glob
import sys
from _old_fit_venus import fit_Xu
from socket import gethostname
import os
import matplotlib.pyplot as plt

sys.path.append("..")
from funciones import donde, angulo, doy_to_day

"""
Agarra los cruces de la MPB de VEX y va a devolverme una lista 
con el SZA de cada uno. Compara con el fit de Xu 2021.
"""
year = 2014
# path = glob.glob("../../../VEX.txt")
if gethostname() == "DESKTOP-2GS0QF2":
    os.chdir(f"G:/VEX{year}/")
    path = glob.glob("*.gz") + glob.glob("pos*.asc")  # los filtrados
    # path = "VEX_mag_filtrado_2014026.gz"
    # path = "pos20140806.asc"
    # os.chdir("C:/Users/RainbowRider/Documents/GitHub/MPB/VEX/")
x_Xu, yz_Xu = fit_Xu()

# calendario = np.zeros((len(path), 3))
calendario = []
timetable = []
angulos = []
# for j in [path]:
for k, j in enumerate(path):
    if "pos" in j:
        pos = np.loadtxt(j)
        hh = pos[:, 3]
        mm = pos[:, 4]
        ss = pos[:, 5]
        orbita = pos[:, 6:9]
        t = hh + mm / 60 + ss / 3600  # hdec
        date = f"{int(pos[1, 0])}-{str(int(pos[1, 1])).zfill(2)}-{str(int(pos[1, 2])).zfill(2)}"
    else:
        MAG = np.loadtxt(j)
        posicion = MAG[4:].T
        t = MAG[0]
        dd = j.split("_")[-1].split(".")[0]
        date = doy_to_day(dd[:4], dd[4:])
        date = f"{date[0]}-{date[1]}-{date[2]}"
        orbita = posicion / 6050  # en RV

    XX = orbita[:, 0]
    YZ = np.sqrt(orbita[:, 1] ** 2 + orbita[:, 2] ** 2)

    """
    Vamos a tirar todos los puntos donde la órbita esté lejos del planeta
    y me quedo solo con el dayside:
    Me quedo con yz < 2, 0 < x < 2
    """
    idx = [i for i in range(len(XX)) if 0 < XX[i] < 2]
    idx2 = [i for i in range(len(YZ)) if YZ[i] < 2]

    indice = list(set(idx) & set(idx2))  # los índices que estén en ambas listas

    x_cut = XX[indice]
    yz_cut = YZ[indice]

    """Para comparar, necesito que las dos listas tengan igual length"""
    if len(x_cut) > 0 and len(yz_cut) > 0:  # si no hay cruce dayside no sigue
        calendario.append(date)
        a = [donde(x_cut, x_Xu[i]) for i in range(len(x_Xu))]  # len(a) = len(Xu)
        pos = np.transpose([x_cut[a], yz_cut[a]])
        pos_xu = np.transpose([x_Xu, yz_Xu])
        # plt.plot(pos, pos_xu)
        # plt.show()

        """
        pos y pos_xu tienen la misma longitud. Como las órbitas son bien portadas
        (valores siempre positivos), puedo simplemente hacer la resta entre las
        normas y ver cuándo es mínima.
        Me quedo con ese punto entonces.
        """

        idx = np.linalg.norm(pos - pos_xu, axis=1).argmin()

        def SZA(posicion, index):
            SZA = angulo(posicion[index, :], [1, 0]) * 180 / np.pi
            return SZA

        sza_mpb = SZA(pos, idx)

        timetable.append(t[idx])
        angulos.append(sza_mpb)

if gethostname() == "DESKTOP-2GS0QF2":
    os.chdir("C:/Users/RainbowRider/Documents/GitHub/MPB/VEX/")

with open(f"../outputs/orbitas_VEX{year}.txt", "a") as file:
    for i in range(len(calendario)):
        file.write(f"{calendario[i]}\t{timetable[i]:1.3g}\t{angulos[i]:1.3g}")
        file.write("\n")

ii = [i for i in range(len(angulos)) if angulos[i] < 65]
with open(f"../outputs/VEX{year}_menor65.txt", "a") as file:
    for i in ii:
        file.write(f"{calendario[i]}\t{timetable[i]:1.3g}\t{angulos[i]:1.3g}")
        file.write("\n")
