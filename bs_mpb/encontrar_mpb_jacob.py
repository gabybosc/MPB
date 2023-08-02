import numpy as np
import matplotlib.pyplot as plt
import sys
from cycler import cycler
from os.path import exists
from figura_encontrar import plot_encontrar

sys.path.append("..")
from funciones import Bpara_Bperp, UTC_to_hdec, donde
from importar_datos import importar_mag_1s, importar_swea, importar_swia


"""
Permite elegir dos tiempos (límites) + uno central para la MPB y los guarda en un txt
Tiene flag
empiezo haciendo los quasipara de todos los grupos
"""

np.set_printoptions(precision=4)
grupo = 2
p = "para"  # para o perp

"""
salteados xq swea anda mal (después se pueden ver aparte):
grupo 1 quasipara num=3, 4, 14, 15 

"""

path = f"outs_catalogo_previa/grupo{grupo}/"

fecha = np.load(path + f"fecha_{p}.npy")
hora_bs = np.load(path + f"hora_{p}.npy")

catalogo = np.genfromtxt(f"../outputs/grupo{grupo}/bs_mpb_final.txt", dtype="str")

plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

# num = int(input("numero de lista\n"))

for num in range(0, len(fecha)):
    print(num)
    year, month, day = fecha[num].split("-")

    t_bs = UTC_to_hdec(hora_bs[num])

    if t_bs > 23 or t_bs < 1:  # no cuenta los que están cerca del cambio de día
        continue
        # num = num + 1
        # print(num)
        # year, month, day = fecha[num].split("-")
        # t_bs = UTC_to_hdec(hora_bs[num])

    ti = t_bs - 1.5  # mira +- 1.5h respecto del BS
    if ti < 0:
        ti = 0.2
    tf = t_bs + 1.5
    if tf > 24:
        tf = 24

    mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
    B_norm = np.linalg.norm(B, axis=1)
    B_para, B_perp_norm, tpara = Bpara_Bperp(B, t, ti + 0.2, tf - 0.2)

    swia, t_swia, i_density, i_temp, vel_mso = importar_swia(year, month, day, ti, tf)

    swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
    energias = [50 + i * 50 for i in range(3)]
    if type(t_swea) != int:
        JE_pds = np.zeros((len(t_swea), len(energias)))

        for i, e in enumerate(energias):
            j = donde(energia, e)
            JE_pds[:, i] = flux_plot[j]
    else:
        JE_pds = 0

    val_MPB = plot_encontrar(
        "MPB",
        fecha[num],
        tpara,
        B_para,
        B_perp_norm,
        t,
        B,
        t_bs,
        B_norm,
        t_swia,
        vel_mso,
        i_density,
        t_swea,
        JE_pds,
        energias,
    )

    flag_MPB = None
    while flag_MPB == None:
        flag = input("MPB confiable? y/n\n")
        if flag == "y":
            flag_MPB = 1
        elif flag == "n":
            flag_MPB = 0

    filepath = f"../outputs/grupo{grupo}/limites_mpb_jacob.txt"

    if not exists(filepath):
        with open(filepath, "w") as file:
            file.write("date\tMPB_min\tMPB\tMPB_max\tflag\n")

    with open(filepath, "a") as file:
        file.write(
            f"{fecha[num]}\t{val_MPB[0]}\t{val_MPB[1]}\t{val_MPB[2]}\t{flag_MPB}\n"
        )
