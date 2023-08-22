import numpy as np
import matplotlib.pyplot as plt


newdates = [
    ("2014", "12", "24"),
    ("2014", "12", "28"),
    ("2014", "12", "29"),
    ("2015", "01", "01"),
    ("2015", "01", "01"),
    ("2015", "01", "02"),
    ("2015", "01", "07"),
    ("2015", "01", "09"),
    ("2015", "02", "12"),
    ("2015", "02", "26"),
]

x_bs = np.array(
    [
        1.23198112,
        1.27474897,
        2.05484307,
        0.97637552,
        1.07587404,
        1.32274189,
        1.09654926,
        1.11391711,
        1.44729145,
        1.23296755,
    ]
)
yz_bs = np.array(
    [
        0.36175449,
        0.82431184,
        1.01575546,
        1.09738858,
        0.95457265,
        1.26839685,
        1.26582335,
        1.24230988,
        1.43468644,
        1.49841115,
    ]
)

x_mpb = np.array(
    [
        0.8074472,
        0.23314956,
        0.58728083,
        0.49102419,
        0.36582596,
        0.41338201,
        0.42607463,
        0.43074661,
        0.32540531,
        0.33611239,
    ]
)

yz_mpb = np.array(
    [
        0.90665237,
        1.12559293,
        1.16655237,
        1.21803787,
        1.22664709,
        1.27315433,
        1.28039496,
        1.28971479,
        1.35886177,
        1.39878099,
    ]
)
names = [i[0] + i[1] + i[2] for i in newdates]

fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.set_xlim(0, 3)
ax.set_ylim(0, 2.5)
circle = plt.Circle((0, 0), 1, color="#c1440e", clip_on=True)
ax.add_artist(circle)
ax.set_title("MAVEN MSO coordinates", fontsize=16)
ax.set_xlabel(r"$X_{MSO}$ ($R_M$)", fontsize=14)
ax.set_ylabel(r"$(Y²_{MSO} + Z²_{MSO} )^{1/2}$ ($R_M$)", fontsize=14)

scatter_bs = ax.scatter(x_bs, yz_bs)
scatter_mpb = ax.scatter(x_mpb, yz_mpb)

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
    if scatter == scatter_bs:
        pos = scatter_bs.get_offsets()[ind["ind"][0]]
    elif scatter == scatter_mpb:
        pos = scatter_mpb.get_offsets()[ind["ind"][0]]
    else:
        return

    annot.xy = pos
    text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont_bs, ind_bs = scatter_bs.contains(event)
        cont_mpb, ind_mpb = scatter_mpb.contains(event)
        if cont_bs:
            update_annot(ind_bs, scatter_bs)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        elif cont_mpb:
            update_annot(ind_mpb, scatter_mpb)
            annot.set_visible(True)

            # Highlight the corresponding point on scatter_bs
            ind_bs_highlight = np.where(
                (x_bs == x_mpb[ind_mpb]) & (yz_bs == yz_mpb[ind_mpb])
            )
            scatter_bs.set_facecolor("none")  # Reset all points to default color
            scatter_bs.set_edgecolor("blue")  # Set default edge color for all points
            scatter_bs.set_edgecolor(
                "red", ind_bs_highlight
            )  # Highlight the corresponding point

            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)

                # Reset the colors of scatter_bs
                scatter_bs.set_facecolor("none")
                scatter_bs.set_edgecolor("blue")

                fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
