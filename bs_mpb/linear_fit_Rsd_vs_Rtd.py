import numpy as np
import matplotlib.pyplot as plt


"""
Necesita antes haber corrido BS_MPB_positions y BS_MPB_Rsd_Rtd
"""

g = input("número de grupo\n")

path = f"../outputs/grupo{g}/"

R_BS = np.load(path + "pos_bs.npy")
R_MPB = np.load(path + "pos_mpb.npy")

Rpolar_BS = np.load(path + "pos_polar_bs.npy")
Rpolar_MPB = np.load(path + "pos_polar_mpb.npy")
Rsd_BS = np.load(path + "Rsd_bs.npy")
Rtd_BS = np.load(path + "Rtd_bs.npy")
Rsd_MPB = np.load(path + "Rsd_mpb.npy")
Rtd_MPB = np.load(path + "Rtd_mpb.npy")
newdates = np.load(path + "newdates.npy")


# CREATE LINEAR REGRESSION


# obtain linear fit slope and origin

m_Rsd, m_Rtd = np.polyfit(Rsd_MPB, Rtd_MPB, 1, cov=False)

b_Rsd, b_Rtd = np.polyfit(Rsd_BS, Rtd_BS, 1, cov=False)


names = [i[0] + "-" + i[1] + "-" + i[2] + "h" + str(i[3])[:3] for i in newdates]
# calculate reduced chi square


def chired(ydata, ymod, param_fit):
    """
    Calcula el chi cuadrado reducido de un ajuste lineal respecto a los parámetros
    (ej: un ajuste lineal, que tiene param_fit = 2)
    """
    sigma_sq = np.mean(ydata**2) - (np.mean(ydata)) ** 2

    chi = np.sum(((ydata - ymod) ** 2) / sigma_sq)

    chi_red = chi / (len(ydata) - param_fit)

    return chi_red


chi_Rsd = np.round(chired(Rsd_BS, m_Rsd * Rsd_MPB + b_Rsd, 2), 3)

chi_Rtd = np.round(chired(Rtd_BS, m_Rtd * Rtd_MPB + b_Rtd, 2), 3)


# PLOT


fig = plt.figure(1)
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))
scatter_sd = ax1.scatter(Rsd_MPB, Rtd_MPB)  # measured data
ax2.scatter(Rsd_BS, Rtd_BS)  # measured data
# ax1.plot(
#     Rsd_MPB, m_Rsd * Rsd_MPB + b_Rsd, color="red", label=r"$\chi^2=${}".format(chi_Rsd)
# )  # model

# scatter_td = ax2.scatter(Rtd_MPB, Rtd_BS)  # measured data
# ax2.plot(
#     Rtd_MPB, m_Rtd * Rtd_MPB + b_Rtd, color="red", label=r"$\chi^2=${}".format(chi_Rtd)
# )  # model
ax1.set_xlabel(r"standoff distance MPB [$R_M$]")
ax1.legend(loc=0)
ax2.set_ylabel(r"terminator distance BS [$R_M$]")
ax2.legend(loc=0)


plt.show()
