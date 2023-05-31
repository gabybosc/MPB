import matplotlib.pyplot as plt
import numpy as np

"""
from https://stackoverflow.com/questions/7908636/how-to-add-hovering-annotations-to-a-plot
cuando tengo un scatter plot y le paso el cursor me devuelve la etiqueta que le haya puesto
"""

# np.random.seed(1)

# x = np.random.rand(15)
# y = np.random.rand(15)
names = np.array(list("ABCDEFGHIJKLMNO" * 4 + "AAAAAA"))

fig, ax = plt.subplots()

sc = plt.scatter(x_mpb, yz_mpb)
plt.scatter(x_bs, yz_bs)


annot = ax.annotate(
    "",
    xy=(0, 0),
    xytext=(20, 20),
    textcoords="offset points",
    bbox=dict(boxstyle="round", fc="w"),
    arrowprops=dict(arrowstyle="->"),
)
annot.set_visible(False)


def update_annot(ind):
    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            update_annot(ind)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
