import numpy as np
import func_position as fpos
import sys

sys.path.append("../../")

"""
busco las posiciones en cartesianas MSO para los tiempos de cruces y después las paso a polares 
"""


# import BS and MPB crossing times

g = input("número de grupo\n")

path = f"../outputs/grupo{g}/"


# find spacecraft position at crossing times
pos_bs = np.transpose(np.load(path + "pos_bs.npy"))
pos_mpb = np.transpose(np.load(path + "pos_mpb.npy"))
r0_MPB = np.array([0.78, 0, 0])
r0_BS = np.array([0.64, 0, 0])


def polarizar(pos, r0):
    Rpolar = np.empty_like(pos)

    for i in range(len(pos)):
        # BS

        x, y, z = pos[i, 0], pos[i, 1], pos[i, 2]

        rho, theta, phi = fpos.cartesian2polar(x, y, z, r0)

        Rpolar[i, 0] = rho
        Rpolar[i, 1] = theta
        Rpolar[i, 2] = phi

    return Rpolar


Rpolar_BS = np.empty_like(pos_bs)

Rpolar_MPB = np.empty_like(pos_mpb)


# focus of conic section from Vignes model

for i in range(len(pos_bs)):
    # BS

    x_BS, y_BS, z_BS = pos_bs[i, 0], pos_bs[i, 1], pos_bs[i, 2]

    rho_BS, theta_BS, phi_BS = fpos.cartesian2polar(x_BS, y_BS, z_BS, r0_BS)

    Rpolar_BS[i, 0] = rho_BS

    Rpolar_BS[i, 1] = theta_BS

    Rpolar_BS[i, 2] = phi_BS

    # MPB

    x_MPB, y_MPB, z_MPB = pos_mpb[i, 0], pos_mpb[i, 1], pos_mpb[i, 2]

    rho_MPB, theta_MPB, phi_MPB = fpos.cartesian2polar(x_MPB, y_MPB, z_MPB, r0_MPB)

    Rpolar_MPB[i, 0] = rho_MPB

    Rpolar_MPB[i, 1] = theta_MPB

    Rpolar_MPB[i, 2] = phi_MPB

# np.save(path + "pos_polar_bs.npy", Rpolar_BS)
# np.save(path + "pos_polar_mpb.npy", Rpolar_MPB)

Rpolar_BS_f = polarizar(pos_bs, r0_BS)
Rpolar_MPB_f = polarizar(pos_mpb, r0_MPB)
