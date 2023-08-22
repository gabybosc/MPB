import numpy as np
import matplotlib.pyplot as plt
import sys
from cycler import cycler
from os.path import exists
from figura_encontrar import plot_encontrar_con_orbita

sys.path.append("..")
from funciones import Bpara_Bperp, UTC_to_hdec, donde, fechas
from importar_datos import importar_mag_1s, importar_swea, importar_swia


"""
Para revisar los cruces que tienen flag 0
"""

np.set_printoptions(precision=4)

grupo = 2  # input("grupo\n")

path = f"outs_catalogo_previa/grupo{grupo}/"

dates = fechas()
time = int(input("hora (HH)\n"))
d_str = str(dates[0] + "-" + dates[1] + "-" + dates[2])

catalogo = np.genfromtxt(f"../outputs/grupo{grupo}/jacob_dayside.txt", dtype="str")

fecha = catalogo[:, 0]
hora_min = catalogo[:, 1]
hora_mpb = catalogo[:, 2]
hora_max = catalogo[:, 3]

for ff in range(len(fecha)):
    if d_str == fecha[ff]:
        if time == int(hora_mpb[ff].split(":")[0]):
            num = ff

plt.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

year, month, day = fecha[num].split("-")
t_mpb = (
    UTC_to_hdec(hora_min[num]),
    UTC_to_hdec(hora_mpb[num]),
    UTC_to_hdec(hora_max[num]),
)

ti = t_mpb[1] - 1
if ti < 0:
    ti = 0.2
tf = t_mpb[1] + 1
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

idx_mpb = donde(t, t_mpb[1])

val_MPB = plot_encontrar_con_orbita(
    fecha[num],
    tpara,
    B_para,
    B_perp_norm,
    t,
    B,
    idx_mpb,
    t_mpb,
    B_norm,
    t_swia,
    posicion,
    i_density,
    t_swea,
    JE_pds,
    energias,
)

flag_MPB = None
while flag_MPB == None:
    flag = input("Sobreescribir MPB? y/n\n")
    if flag == "y":
        flag_MPB = 1
    elif flag == "n":
        flag_MPB = 0

filepath = f"../outputs/grupo{grupo}/jacob_dayside_revised.txt"

if not exists(filepath):
    with open(filepath, "w") as file:
        file.write("date\tMPB_min\tMPB\tMPB_max\tflag\ttheta\tbeta\n")

with open(filepath, "a") as file:
    file.write(
        f"{fecha[num]}\t{val_MPB[0]}\t{val_MPB[1]}\t{val_MPB[2]}\t{flag_MPB}\t{catalogo[num, -2]}\t{catalogo[num, -1]}\n"
    )
