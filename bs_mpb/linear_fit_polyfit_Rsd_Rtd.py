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


fig = plt.figure(1)
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.1, wspace=0.15
)
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))
scatter_sd = ax1.scatter(Rsd_MPB, Rsd_BS)  # measured data
ax1.plot(
    Rsd_MPB, m_Rsd * Rsd_MPB + b_Rsd, color="red", label=r"$\chi^2=${}".format(chi_Rsd)
)  # model

scatter_td = ax2.scatter(Rtd_MPB, Rtd_BS)  # measured data
ax2.plot(
    Rtd_MPB, m_Rtd * Rtd_MPB + b_Rtd, color="red", label=r"$\chi^2=${}".format(chi_Rtd)
)  # model
ax1.set_xlabel(r"standoff distance MPB [$R_M$]")
ax1.set_ylabel(r"standoff distance BS [$R_M$]")
ax1.legend(loc=0)
ax2.set_xlabel(r"terminator distance MPB [$R_M$]")
ax2.set_ylabel(r"terminator distance BS [$R_M$]")
ax2.legend(loc=0)


annot = ax1.annotate(
    "",
    xy=(0, 0),
    xytext=(20, 20),
    textcoords="offset points",
    bbox=dict(boxstyle="round", fc="w"),
    arrowprops=dict(arrowstyle="->"),
)
annot.set_visible(False)

annot = ax2.annotate(
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
    if event.inaxes == ax1:
        cont_sd, ind_sd = scatter_sd.contains(event)
        if cont_sd:
            update_annot(ind_sd, scatter_sd)
            annot.set_visible(True)
            scatter_sd.set_facecolor("C0")  # Reset facecolor of all points in ax1
            fig.canvas.draw_idle()
            if ind_sd["ind"][0] < len(Rsd_BS):
                ind_td = {"ind": [ind_sd["ind"][0]]}
                update_annot(ind_td, scatter_td)
                scatter_td.set_facecolor(
                    ["red" if i in ind_td["ind"] else "C0" for i in range(len(Rtd_BS))]
                )  # Highlight the corresponding point in ax2
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()
    elif event.inaxes == ax2:
        cont_td, ind_td = scatter_td.contains(event)
        if cont_td:
            update_annot(ind_td, scatter_td)
            annot.set_visible(True)
            scatter_td.set_facecolor("C0")  # Reset facecolor of all points in ax2
            fig.canvas.draw_idle()
            if ind_td["ind"][0] < len(Rtd_MPB):
                ind_sd = {"ind": [ind_td["ind"][0]]}
                update_annot(ind_sd, scatter_sd)
                scatter_sd.set_facecolor(
                    ["red" if i in ind_sd["ind"] else "C0" for i in range(len(Rsd_MPB))]
                )  # Highlight the corresponding point in ax1
        else:
            if vis:
                annot.set_visible(False)
                scatter_sd.set_facecolor("C0")  # Reset facecolor of all points in ax1
                fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
