import numpy as np

# import matplotlib.pyplot as plt
from os.path import exists
import sys

sys.path.append("..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde, SZA

# grupo = input("grupo\n")
for grupo in [2, 3, 4]:
    lista = np.genfromtxt(
        f"../outputs/grupo{grupo}/jacob_todo.txt", skip_header=1, dtype=str
    )

    """
    date	MPB_min	MPB	MPB_max	flag	theta	beta
    2015-09-08	07:56:40	07:58:16	08:02:38	1	9.86699	6.0438123
    """

    # pos_mpb = []
    # newdates = []
    final = []

    for l in lista:
        year, month, day = l[0].split("-")

        t_mpb_min = UTC_to_hdec(l[1])
        t_mpb_max = UTC_to_hdec(l[3])
        t_mpb = UTC_to_hdec(l[2])

        ti = t_mpb_min - 0.2
        tf = t_mpb_max + 0.2

        if ti < 0:
            ti = 0
        if tf > 24:
            tf = 24

        mag, t, B, pos = importar_mag_1s(year, month, day, ti, tf)
        idx_mpb = donde(t, t_mpb)

        sza_mpb = SZA(pos, idx_mpb)
        if sza_mpb < 85:
            # if pos[idx_mpb][0] > 0:
            # pos_mpb.append(pos[idx_mpb] / 3390)
            # newdates.append((year, month, day))
            final.append(l)

    # pos_mpb = np.transpose(pos_mpb)
    filepath = f"../outputs/grupo{grupo}/jacob_dayside.txt"

    if not exists(filepath):
        with open(filepath, "w") as file:
            file.write("date\tMPB_min\tMPB\tMPB_max\tflag\ttheta\tbeta\n")

    with open(filepath, "a") as file:
        for f in final:
            file.write(f"{f[0]}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}\t{f[5]}\t{f[6]}\n")
