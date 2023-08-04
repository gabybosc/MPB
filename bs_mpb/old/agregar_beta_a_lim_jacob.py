import numpy as np
import sys
from os.path import exists

sys.path.append("..")
from funciones import donde


for g in [3]:
    cat = np.genfromtxt(
        f"../outputs/grupo{g}/limites_mpb_jacob_para.txt", skip_header=1, dtype="str"
    )
    beta = np.load(f"outs_catalogo_previa/grupo{g}/pdyn_perp.npy")
    theta = np.load(f"outs_catalogo_previa/grupo{g}/theta_perp.npy")
    fecha = np.load(f"outs_catalogo_previa/grupo{g}/fecha_perp.npy")
    hora = np.load(f"outs_catalogo_previa/grupo{g}/hora_perp.npy")

    filepath = f"../outputs/grupo{g}/limites_mpb_jacob_para2.txt"
    if not exists(filepath):
        with open(filepath, "w") as file:
            file.write("date\tMPB_min\tMPB\tMPB_max\tflag\ttheta\tbeta\n")

    for i in range(len(cat)):
        for j in range(len(fecha)):
            if cat[i, 0] == fecha[j]:
                tcat = int(cat[i, 1].split(":")[:1][0])  # hora del cat
                tnpy = int(hora[j].split(":")[:1][0])  # hora del npy
                if tcat == tnpy or tcat == tnpy - 1 or tcat == tnpy + 1:
                    # print(cat[i, 0], tcat, fecha[j], tnpy)
                    with open(filepath, "a") as file:
                        file.write(
                            f"{cat[i, 0]}\t{cat[i, 1]}\t{cat[i, 2]}\t{cat[i, 3]}\t{cat[i, 4]}\t{theta[j]}\t{beta[j]}\n"
                        )
