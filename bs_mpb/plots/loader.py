import numpy as np


def importar_posiciones(path):
    date = np.load(path + "date.npy")
    t_bs = np.load(path + "BS_t.npy")
    t_mpb = np.load(path + "MPB_t.npy")
    times = np.vstack((t_bs, t_mpb))
    pos_bs = np.transpose(
        np.vstack(
            (
                np.load(path + "BS_x.npy").astype(float),
                np.load(path + "BS_y.npy").astype(float),
                np.load(path + "BS_z.npy").astype(float),
            )
        )
    )
    pos_mpb = np.transpose(
        np.vstack(
            (
                np.load(path + "MPB_x.npy").astype(float),
                np.load(path + "MPB_y.npy").astype(float),
                np.load(path + "MPB_z.npy").astype(float),
            )
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
    beta = np.load(path + "beta_p.npy", allow_pickle=True).astype(float)
    cone_angle = np.load(path + "cone_angle.npy").astype(float)
    Mfms = np.load(path + "Mfms.npy")
    Ls = np.load(path + "Ls.npy")
    pdyn = np.load(path + "pdyn.npy")

    return beta, cone_angle, Mfms, Ls, pdyn


def importar_tiempos(path):
    date = np.load(path + "date.npy")
    BS_t_min = np.load(path + "BS_t_min.npy")
    BS_t = np.load(path + "BS_t.npy")
    BS_t_max = np.load(path + "BS_t_max.npy")
    MPB_t_min = np.load(path + "MPB_t_min.npy")
    MPB_t_max = np.load(path + "MPB_t_max.npy")
    MPB_t = np.load(path + "MPB_t.npy")

    t_bs = np.vstack((BS_t_min, BS_t, BS_t_max))
    t_mpb = np.vstack((MPB_t_min, MPB_t, MPB_t_max))

    return date, t_bs, t_mpb
