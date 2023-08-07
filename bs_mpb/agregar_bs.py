import numpy as np
import sys
from os.path import exists

sys.path.append("..")
from funciones import donde, UTC_to_hdec

grupo = input("grupo?\n")
# p = input("para/perp?\n")
# grupo = 1
p = "perp"
path = f"../outputs/grupo{grupo}/jacob_dayside.txt"
final = f"../outputs/grupo{grupo}/jacob_dayside_bs.txt"

datos = np.genfromtxt(path, dtype="str", skip_header=1)
fechas = np.load(f"outs_catalogo_previa/grupo{grupo}/fecha_{p}.npy")
horas = np.load(f"outs_catalogo_previa/grupo{grupo}/hora_{p}.npy")

if not exists(final):
    with open(final, "w") as file:
        file.write("date\tt_bs\tMPB_min\tMPB\tMPB_max\tflag\ttheta\tbeta\n")

for i in range(len(datos[:, 0])):
    for j in range(len(fechas)):
        if fechas[j] == datos[i, 0]:
            bs_utc = UTC_to_hdec(horas[j])
            mpb_utc = UTC_to_hdec(datos[i, 1])
            if bs_utc - 1 < mpb_utc < bs_utc + 1:
                with open(final, "a") as file:
                    file.write(
                        f"{datos[i, 0]}\t{horas[j]}\t{datos[i, 1]}\t{datos[i, 2]}\t{datos[i, 3]}\t{datos[i, 4]}\t{datos[i, 5]}\t{datos[i, 6]}\n"
                    )
