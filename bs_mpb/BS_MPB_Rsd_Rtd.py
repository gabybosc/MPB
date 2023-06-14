import numpy as np
import Vignesfit_functions as fvig


"""
Calculo las standoff y terminator distance (uso funciones que están en Vignes_functions <-- ahí está la definición del terminator para chequear).
"""


g = input("número de grupo\n")

path = f"../outputs/grupo{g}/"

R_BS = np.transpose(np.load(path + "pos_bs.npy"))
R_MPB = np.transpose(np.load(path + "pos_mpb.npy"))

Rpolar_BS = np.load(path + "pos_polar_bs.npy")

Rpolar_MPB = np.load(path + "pos_polar_mpb.npy")


N_events = len(R_BS[:, 0])


# ADJUST L VIGNES TO FIT BS & MPB CROSSINGS


L_BS = np.empty(N_events)

L_MPB = np.empty(N_events)

for i in range(N_events):
    L_BS[i] = fvig.fit_L(Rpolar_BS[i, 0], Rpolar_BS[i, 1])

    L_MPB[i] = fvig.fit_L(Rpolar_MPB[i, 0], Rpolar_MPB[i, 1], 0.90)


# CALCULATE VIGNES STANDOFF AND TERMINATOR DISTANCE FOR BS & MPB CROSSINGS


Rsd_BS = np.empty(N_events)

Rtd_BS = np.empty(N_events)

Rsd_MPB = np.empty(N_events)

Rtd_MPB = np.empty(N_events)

for i in range(N_events):
    Rsd_BS[i] = fvig.Rsd(L_BS[i])

    Rtd_BS[i] = fvig.Rtd(L_BS[i], R_BS[i, 1], R_BS[i, 2])

    Rsd_MPB[i] = fvig.Rsd(L_MPB[i], 0.78, 0.90)

    Rtd_MPB[i] = fvig.Rtd(L_MPB[i], R_MPB[i, 1], R_MPB[i, 2], 0.78, 0.90)


np.save(path + "Rsd_bs.npy", Rsd_BS)
np.save(path + "Rtd_bs.npy", Rtd_BS)
np.save(path + "Rsd_mpb.npy", Rsd_MPB)
np.save(path + "Rtd_mpb.npy", Rtd_MPB)
