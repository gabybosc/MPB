import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("..")
from importar_datos import importar_mag_1s
from funciones import UTC_to_hdec, donde, SZA

grupo = input("grupo\n")
lista = np.genfromtxt(
    f"../outputs/grupo{grupo}/limites_mpb_jacob.txt", skip_header=1, dtype=str
)
print(len(lista))


pos_bs = []
pos_mpb = []
newdates = []
final = []

for l in lista:
    year, month, day = l[0].split("-")

    t_bs = UTC_to_hdec(l[1])
    t_mpb = UTC_to_hdec(l[2])

    if t_bs < t_mpb:
        ti = t_bs - 0.2
        tf = t_mpb + 0.2
    else:
        ti = t_mpb - 0.2
        tf = t_bs + 0.2
    if ti < 0:
        ti = 0
    if tf > 24:
        tf = 24

    mag, t, B, pos = importar_mag_1s(year, month, day, ti, tf)
    idx_bs = donde(t, t_bs)
    idx_mpb = donde(t, t_mpb)
    sza_bs = SZA(pos, idx_bs)
    sza_mpb = SZA(pos, idx_mpb)

    # if pos[idx_mpb][0] > 0:
    if sza_mpb < 90:
        pos_bs.append(pos[idx_bs] / 3390)
        pos_mpb.append(pos[idx_mpb] / 3390)
        newdates.append((year, month, day))
        final.append(l)


pos_bs = np.transpose(pos_bs)
pos_mpb = np.transpose(pos_mpb)

print(len(final))
with open(f"../outputs/grupo{grupo}/jacob_final.txt", "a") as file:
    for f in final:
        file.write(f"{f[0]}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}")
        file.write("\n")
