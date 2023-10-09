import numpy as np
import Vignesfit_functions as fvig


"""
Necesita antes haber corrido BS_MPB_positions 
Calculo las standoff y terminator distance (uso funciones que están en Vignes_functions <-- ahí está la definición del terminator para chequear).
"""


def Rsd_Rtd(pos, pos_polar, x0, eps):
    """
    Encuentra la standoff distance y terminator distance dado un array de posiciones
    x0, eps son los parámetros de la cónica
    """

    N_events = len(pos)
    L = np.empty(N_events)

    for i in range(N_events):
        L[i] = fvig.fit_L(pos_polar[i, 0], pos_polar[i, 1], eps)

    # CALCULATE VIGNES STANDOFF AND TERMINATOR DISTANCE

    Rsd = np.empty(N_events)
    Rtd = np.empty(N_events)

    for i in range(N_events):
        Rsd[i] = fvig.Rsd(L[i], x0, eps)

        Rtd[i] = fvig.Rtd(L[i], x0, eps)

    return Rsd, Rtd


g = input("número de grupo\n")

path = f"../outputs/grupo{g}/"

R_BS = np.load(path + "pos_bs.npy")
R_MPB = np.load(path + "pos_mpb.npy")

Rpolar_BS = np.load(path + "pos_polar_bs.npy")

Rpolar_MPB = np.load(path + "pos_polar_mpb.npy")


N_events = len(R_BS)


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

    Rtd_BS[i] = fvig.Rtd(L_BS[i])

    Rsd_MPB[i] = fvig.Rsd(L_MPB[i], 0.78, 0.90)

    Rtd_MPB[i] = fvig.Rtd(L_MPB[i], 0.78, 0.90)


np.save(path + "Rsd_bs.npy", Rsd_BS)
np.save(path + "Rtd_bs.npy", Rtd_BS)
np.save(path + "Rsd_mpb.npy", Rsd_MPB)
np.save(path + "Rtd_mpb.npy", Rtd_MPB)


Rsd_BS_f, Rtd_BS_f = Rsd_Rtd(R_BS, Rpolar_BS, 0.64, 1.03)
Rsd_MPB_f, Rtd_MPB_f = Rsd_Rtd(R_MPB, Rpolar_MPB, 0.78, 0.9)

Rsd_BS_f == Rsd_BS
Rtd_BS_f == Rtd_BS
Rsd_MPB_f == Rsd_MPB
Rtd_MPB_f == Rtd_MPB

# Lf_bs = L_func(R_BS, Rpolar_BS)
# Lf_mpb = L_func(R_MPB, Rpolar_MPB)

# Lf_mpb == L_MPB


# def L_func(pos, pos_polar, eps):
#     """
#     Encuentra la standoff distance y terminator distance dado un array de posiciones
#     x0, eps son los parámetros de la cónica
#     """

#     N_events = len(pos[:, 0])
#     L = np.empty(N_events)

#     for i in range(N_events):
#         L[i] = fvig.fit_L(pos_polar[i, 0], pos_polar[i, 1], eps)
#     return L


# def Rsd(pos, L, x0, eps):
#     # CALCULATE VIGNES STANDOFF AND TERMINATOR DISTANCE

#     Rsd = np.empty(N_events)
#     Rtd = np.empty(N_events)

#     for i in range(N_events):
#         Rsd[i] = fvig.Rsd(L[i], x0, eps)

#         Rtd[i] = fvig.Rtd(L[i], pos[i, 1], pos[i, 2], x0, eps)

#     return Rsd, Rtd
