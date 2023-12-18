import numpy as np


def importar_posiciones(path):
    date = np.load(path + "date.npy")
    t_bs = np.load(path + "BS_t.npy")
    t_mpb = np.load(path + "MPB_t.npy")
    times = np.vstack((t_bs, t_mpb))
    pos_bs = np.vstack(
        (
            np.load(path + "BS_x.npy").astype(float),
            np.load(path + "BS_y.npy").astype(float),
            np.load(path + "BS_z.npy").astype(float),
        )
    )
    pos_mpb = np.vstack(
        (
            np.load(path + "MPB_x.npy").astype(float),
            np.load(path + "MPB_y.npy").astype(float),
            np.load(path + "MPB_z.npy").astype(float),
        )
    )
    Rsd = np.vstack(
        (
            np.load(path + "Rsd_BS.npy").astype(float),
            np.load(path + "Rsd_MPB.npy").astype(float),
        )
    )

    return date, times, pos_bs, pos_mpb, Rsd


def importar_params(path):
    beta = np.load(path + "beta_p.npy").astype(float)
    cone_angle = np.load(path + "cone_angle.npy").astype(float)
    Mfms = np.load(path + "Mfms.npy")
    Ls = np.load(path + "Ls.npy")
    pdyn = np.load(path + "pdyn.npy")

    return beta, cone_angle, Mfms, Ls, pdyn
