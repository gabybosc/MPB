import numpy as np
import sys
from os.path import exists

sys.path.append("..")
from funciones import donde, UTC_to_hdec

# grupo = input("grupo?\n")
# grupo = 1

for grupo in [1, 2, 3, 4]:
    path = f"../outputs/grupo{grupo}/jacob_dayside.txt"
    final = f"../outputs/grupo{grupo}/jacob_dayside_bs.txt"

    datos = np.genfromtxt(path, dtype="str", skip_header=1)
    fechas_para = np.load(f"outs_catalogo_previa/grupo{grupo}/fecha_para.npy")
    horas_para = np.load(f"outs_catalogo_previa/grupo{grupo}/hora_para.npy")
    fechas_perp = np.load(f"outs_catalogo_previa/grupo{grupo}/fecha_perp.npy")
    horas_perp = np.load(f"outs_catalogo_previa/grupo{grupo}/hora_perp.npy")

    fechas = np.concatenate((fechas_para, fechas_perp))
    horas = np.concatenate((horas_para, horas_perp))

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
