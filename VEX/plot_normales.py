import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import axes3d  # lo usa aunque esté marcado como que no!
from _fit import fit_3d  # , plot_3d
from _importar_datos import importar_MAG

sys.path.append("..")
from funciones import find_nearest, UTC_to_hdec

"""
plotea todas las normales en 3D
"""


def plot_3d(fig, x, y, z):
    # nmva = [0.517, 0.103, 0.850]
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_xlabel(r"$X_{VSO} (R_V)$")
    ax.set_ylabel(r"$Y_{VSO} (R_V)$")
    ax.set_zlabel(r"$Z_{VSO} (R_V)$")
    ax.plot_wireframe(x, y, z, color="gray", alpha=0.5, linewidth=0.5)

    # Ajustar la escala de los ejes para que sean iguales
    max_range = (
            np.array([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]).max() / 2.0
    )
    mid_x = (x.max() + x.min()) * 0.5
    mid_y = (y.max() + y.min()) * 0.5
    mid_z = (z.max() + z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    plt.legend()
    return ax


def plot_normales(ax, R, nmva, nfit):
    ax.quiver(
        R[0],
        R[1],
        R[2],
        nfit[0],
        nfit[1],
        nfit[2],
        color="#ffa600",
        length=0.5,
        label="Ajuste",
    )
    ax.quiver(
        R[0],
        R[1],
        R[2],
        nmva[0],
        nmva[1],
        nmva[2],
        color="#de425b",
        length=0.5,
        label="MVA",
    )


x3 = [
    np.array([0.655, 0.639, 0.403]),
    np.array([0.935, 0.349, 0.055]),
    np.array([0.581, -0.189, 0.792]),
    np.array([0.730, 0.251, 0.636]),
    np.array([0.679, 0.325, 0.658]),
    np.array([0.770, 0.292, 0.567]),
    np.array([0.865, 0.297, 0.404]),
    np.array([0.517, 0.103, 0.85]),
]
dates = [
    [2008, 281],  # , 10, "07"],
    [2008, 50],  # "02", 19],
    [2008, 62],  # "03", "02"],
    [2008, 53],  # "02", 22],
    [2008, 271],  # , "09", 27],
    [2014, 109],  # , "04", 19],
    [2014, 116],  # , "04", 26],
    [2008, 302],  # , 10, 28],
]

tmean = [7.7865, 4.5269, 4.3899, 4.5427, 7.4257, 2.3766, 2.3257, 8.5479]
nfit = [
    np.array([0.590, 0.164, 0.791]),
    np.array([0.695, 0.318, 0.645]),
    np.array([0.619, 0.093, 0.78]),
    np.array([0.559, 0.205, 0.803]),
    np.array([0.524, 0.283, 0.803]),
    np.array([0.791, 0.247, 0.560]),
    np.array([0.869, 0.115, 0.481]),
    np.array([0.587, -0.116, 0.801]),
]
nmva = [
    np.array([0.655, 0.639, 0.403]),
    np.array([0.935, 0.349, 0.055]),
    np.array([0.581, -0.189, 0.792]),
    np.array([0.730, 0.251, 0.636]),
    np.array([0.679, 0.325, 0.658]),
    np.array([0.770, 0.292, 0.567]),
    np.array([0.865, 0.297, 0.404]),
    np.array([0.517, 0.103, 0.850]),
]

fig = plt.figure()
x, y, z = fit_3d()
ax = plot_3d(fig, x, y, z)
for i in range(len(x3)):
    t, B, pos, cl, tpos = importar_MAG(
        dates[i][0], dates[i][1], tmean[i] - 0.2, tmean[i] + 0.2
    )
    index = np.where(tpos == find_nearest(tpos, tmean[i]))[0][0]
    if np.linalg.norm(pos[index, :]) < 200:
        R = pos[index, :]
    else:
        R = pos[index, :] / 6050  # la posicion de la nave en RV
    plot_normales(ax, R, nmva[i], nfit[i])
    if i == 0:
        plt.legend()

plt.show()

"""
fecha           & tiempo    & h23   & h14   & lon i & rg    & jv     & FH                        & pdyn      & pmag \\
19 feb 2008     & 04:41:36  & 112   & 581   & $71$  & $80$  & 449    & $2.51 \times 10^{-14}$    & $1.12$    & $0.67$\\
22 feb 2008     & 04:32:33  & 264   & 539   & $117$ & $100$ & 52     & $2.26 \times 10^{-15}$    & $0.07$    & $0.61$\\
2 mar 2008      & 04:23:23  & 178   & 474   & $42$  & $25$  & 121    & $5.36 \times 10^{-15}$    & $0.60$    & $0.83$\\
27 sept 2008    & 07:25:32  & 337   & 547   & $57$  & $125$ & 88     & $2.50 \times 10^{-15}$    & $0.39$    & $0.57$\\
7 oct 2008      & 07:47:11  & 77    & 343   & $60$  & $54$  & 128    & $4.15 \times 10^{-15}$    & $0.67$    & $0.42$ \\
28 oct 2008     & 08:32:52  & 185   & 380   & $52$  & $274$ & 143    & $6.01 \times 10^{-15}$    & $0.54$    & $0.74$\\
19 abr 2014     & 02:22:35  & 252   & 459   & $117$ & -\footnotemark[1] & 56    & $1.81 \times 10 ^{-15}$   & $0.16$    & $0.65$\\
26 abr 2014     & 02:19:32  & 296   & 586   & $121$ & -\footnotemark[1] & 72    & $3.45 \times 10^{-15}$    & $1.55$    & $1.10$

fecha           & tiempo    & sza   & tv    & tB    & h23   & h14   & lon i & rg    \\
19 feb 2008     & 04:41:36  & 51.1  & 154   & 86    & 112   & 581   & $71$  & $80$  \\
22 feb 2008     & 04:32:33  & 62.1  & 125   & 102   & 264   & 539   & $117$ & $100$ \\
2 mar 2008      & 04:23:23  & 57.5  & 108   & 101   & 178   & 474   & $42$  & $25$  \\
27 sept 2008    & 07:25:32  & 64.8  & 122   & 82    & 337   & 547   & $57$  & $125$ \\
7 oct 2008      & 07:47:11  & 59.8  & 131   & 75    & 77    & 343   & $60$  & $54$  \\
28 oct 2008     & 08:32:52  & 60.0  & 102   & 96    & 185   & 380   & $52$  & $274$ \\
19 abr 2014     & 02:22:35  & 41.9  & 110   & 91    & 252   & 459   & $117$ & -\footnotemark[1] \\
26 abr 2014     & 02:19:32  & 33.0  & 116   & 92    & 296   & 586   & $121$ & -\footnotemark[1] 

fecha           & tiempo    & h23    & js   & jv     & FH                        & pdyn      & pmag \\
19 feb 2008     & 04:41:36  & 112    & 50.4 & 449    & $2.51 \times 10^{-14}$    & $1.12$    & $0.67$\\
22 feb 2008     & 04:32:33  & 264    & 13.7 & 52     & $2.26 \times 10^{-15}$    & $0.07$    & $0.61$\\
2 mar 2008      & 04:23:23  & 178    & 21.5 & 121    & $5.36 \times 10^{-15}$    & $0.60$    & $0.83$\\
27 sept 2008    & 07:25:32  & 337    & 29.7 & 88     & $2.50 \times 10^{-15}$    & $0.39$    & $0.57$\\
7 oct 2008      & 07:47:11  & 77     & 9.9  & 128    & $4.15 \times 10^{-15}$    & $0.67$    & $0.42$ \\
28 oct 2008     & 08:32:52  & 185    & 26.4 & 143    & $6.01 \times 10^{-15}$    & $0.54$    & $0.74$\\
19 abr 2014     & 02:22:35  & 252    & 14.2 & 56     & $1.81 \times 10 ^{-15}$   & $0.16$    & $0.65$\\
26 abr 2014     & 02:19:32  & 296    & 21.5 & 72     & $3.45 \times 10^{-15}$    & $1.55$    & $1.10$

"""
