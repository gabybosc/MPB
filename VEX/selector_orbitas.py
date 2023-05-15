import numpy as np
import glob as glob
import sys
from fit_venus import fit_Xu

sys.path.append("..")
from funciones import donde, angulo


"""
Agarra los cruces de la MPB de VEX y va a devolverme una lista 
con el SZA de cada uno. Compara con el fit de Xu 2021.
"""

# path = glob.glob("../../../VEX.txt")

path = glob.glob("../../../datos/VEX/2009/*.tab")
# Ajuste de Xu de la MPB:
x_Xu, yz_Xu = fit_Xu()


# calendario = np.zeros((len(path), 3))
calendario = []
timetable = []
angulos = []
for k, j in enumerate(path):
    posicion = np.loadtxt(j, skiprows=1, usecols=[8, 9, 10])
    orbita = posicion / 6050  # en RV

    timedata = np.loadtxt(j, skiprows=1, usecols=0, dtype=str)
    time = np.array([timedata[i].split("T")[1] for i in range(len(timedata))])
    date = timedata[0].split("T")[0]

    XX = orbita[:, 0]
    YZ = np.sqrt(orbita[:, 1] ** 2 + orbita[:, 2] ** 2)

    # calendario[k, 0] = date  # guarda la fecha

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

        # calendario[k, 1] = time[idx]
        # calendario[k, 2] = sza_mpb
        timetable.append(time[idx])
        angulos.append(sza_mpb)


with open("../outputs/orbitas_VEX.txt", "a") as file:
    for i in range(len(calendario)):
        file.write(f"{calendario[i]}\t{timetable[i]}\t{angulos[i]}")
        file.write("\n")
