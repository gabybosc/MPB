import numpy as np
from os.path import exists


for grupo in [1, 2, 3, 4]:
    mpb_path = f"../../Documents/sofibis/jacob_dayside_grupo{grupo}.txt"
    bs_path = f"../../Documents/sofibis/jacob_dayside_bs_limites_grupo{grupo}.txt"
    final = f"../outputs/grupo{grupo}/FINAL_ESTA_SI.txt"

    mpb = np.genfromtxt(mpb_path, dtype="str", skip_header=1)
    fecha = mpb[:, 0]
    mpb_min = mpb[:, 1]
    mpb_max = mpb[:, 3]
    hora_mpb = mpb[:, 2]
    flag_mpb = mpb[:, 4]

    bs = np.genfromtxt(bs_path, dtype="str", skip_header=1)
    fecha_bs = bs[:, 0]
    bs_min = bs[:, 1]
    bs_max = bs[:, 3]
    hora_bs = bs[:, 2]
    flag_bs = bs[:, 4]

    if not exists(final):
        with open(final, "w") as file:
            file.write(
                "date_bs\tdate_mpb\tBS_min\tBS\tBS_max\tflag_BS\tMPB_min\tMPB\tMPB_max\tflag_MPB\ttheta\tbeta\n"
            )

    for i in range(len(mpb[:, 0])):
        with open(final, "a") as file:
            file.write(
                f"{fecha_bs[i]}\t{fecha[i]}\t{bs_min[i]}\t{hora_bs[i]}\t{bs_max[i]}\t{flag_bs[i]}\t{mpb_min[i]}\t{hora_mpb[i]}\t{mpb_max[i]}\t{flag_mpb[i]}\t{mpb[i, 5]}\t{mpb[i, 6]}\n"
            )
