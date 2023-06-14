import numpy as np
import matplotlib.pyplot as plt


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

m_Rsd, b_Rsd = np.polyfit(Rsd_MPB, Rsd_BS, 1, cov=False)

m_Rtd, b_Rtd = np.polyfit(Rtd_MPB, Rtd_BS, 1, cov=False)


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


# standoff distance

fig, ax = plt.subplots()
scatter_sd = ax.scatter(Rsd_MPB, Rsd_BS)  # measured data
ax.plot(
    Rsd_MPB, m_Rsd * Rsd_MPB + b_Rsd, color="red", label=r"$\chi^2=${}".format(chi_Rsd)
)  # model
ax.set_xlabel(r"standoff distance MPB [$R_M$]")
ax.set_ylabel(r"standoff distance BS [$R_M$]")
ax.legend(loc=0)
# plt.savefig(path_main+path_Rsd_Rtd+r'/plots/Group{}_standoff_polyfit.png'.format(num_group), dpi=300)


# terminator distance

fig, ax = plt.subplots()
scatter_td = ax.scatter(Rtd_MPB, Rtd_BS)  # measured data
ax.plot(
    Rtd_MPB, m_Rtd * Rtd_MPB + b_Rtd, color="red", label=r"$\chi^2=${}".format(chi_Rtd)
)  # model
ax.set_xlabel(r"terminator distance MPB [$R_M$]")
ax.set_ylabel(r"terminator distance BS [$R_M$]")
ax.legend(loc=0)
# plt.savefig(path_main+path_Rsd_Rtd+r'/plots/Group{}_terminator_polyfit.png'.format(num_group), dpi=300)


annot = ax.annotate(
    "",
    xy=(0, 0),
    xytext=(20, 20),
    textcoords="offset points",
    bbox=dict(boxstyle="round", fc="w"),
    arrowprops=dict(arrowstyle="->"),
)
annot.set_visible(False)


def update_annot(ind, scatter):
    pos = scatter.get_offsets()[ind["ind"][0]]

    annot.xy = pos
    text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont_bs, ind_bs = scatter_td.contains(event)
        cont_mpb, ind_mpb = scatter_sd.contains(event)
        if cont_bs:
            update_annot(ind_bs, scatter_td)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        elif cont_mpb:
            update_annot(ind_mpb, scatter_sd)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
