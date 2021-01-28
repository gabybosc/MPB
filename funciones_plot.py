import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.widgets import RectangleSelector
from matplotlib.lines import Line2D
from funciones import array_datenums
from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()

"""
Acá las funciones para hacer plots en general.
"""


def axisEqual3D(ax):
    extents = np.array([getattr(ax, "get_{}lim".format(dim))() for dim in "xyz"])
    sz = extents[:, 1] - extents[:, 0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize / 2
    for ctr, dim in zip(centers, "xyz"):
        getattr(ax, "set_{}lim".format(dim))(ctr - r, ctr + r)


def deltat():
    global x1
    global x2
    return x1, x2


def hodograma(B1, B2, B3):
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), sharex=True)

    ax1.plot(B2, B1, zorder=1)
    ax2.plot(B3, B1, zorder=1)
    ax1.scatter(B2[0], B1[0], s=50, zorder=2, marker="o", color="r", label="start")
    ax1.scatter(B2[-1], B1[-1], s=50, zorder=2, marker="x", color="r", label="end")
    ax2.scatter(B3[0], B1[0], s=50, zorder=2, marker="o", color="r", label="start")
    ax2.scatter(B3[-1], B1[-1], s=50, zorder=2, marker="x", color="r", label="end")
    ax1.set_xlabel(f"B2 (nT)", fontsize=16)
    ax2.set_xlabel(f"B3 (nT)", fontsize=16)
    ax1.set_ylabel(f"B1 (nT)", fontsize=16)
    ax2.set_ylabel(f"B1 (nT)", fontsize=16)
    ax1.grid()
    ax2.grid()
    ax1.tick_params(axis="both", which="major", labelsize=14)
    ax2.tick_params(axis="both", which="major", labelsize=14)
    plt.suptitle("MAVEN MAG MVA 10-10-2015", fontsize=18)
    plt.legend(fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])


def imshow_UTC(
    year, month, day, t, heatmap, eje_y, colormap="inferno", minimo=0, maximo=20
):
    """ Le das una fecha en np.datetime64 (UTC) y te grafica el imshow.    """
    timestamps = array_datenums(year, month, day, t)
    t_graph = md.date2num(timestamps)

    plt.subplots_adjust(bottom=0.2)
    plt.xticks(rotation=25)
    ax = plt.gca()
    xfmt = md.DateFormatter("%H:%M:%S")
    ax.xaxis.set_major_formatter(xfmt)
    plt.imshow(
        heatmap,
        aspect="auto",
        origin="lower",
        extent=(t_graph[0], t_graph[-1], eje_y[0], eje_y[-1]),
        cmap=colormap,
        vmin=minimo,
        vmax=maximo,
    )
    plt.colorbar()


def line_select_callback(eclick, erelease):
    global x1, x2
    "eclick and erelease are the press and release events"
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    # t1, t2 = x1, x2


def make_colormap(vmin, vmax, *args):
    """    args = list of colors, in order from vmin to vmax    """
    colors_map_default = [
        "black",
        "darkblue",
        "mediumseagreen",
        "green",
        "red",
        "orange",
        "yellow",
    ]

    if not args:
        args = colors_map_default

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    bins = np.linspace(vmin, vmax, len(args))

    colors = [[norm(bins[i]), args[i]] for i in range(len(args))]
    colormap = mpl.colors.LinearSegmentedColormap.from_list("", colors)

    return colormap


def onpick1(event):
    if isinstance(event.artist, Line2D):
        thisline = event.artist
        xdata = thisline.get_xdata()
        # ydata = thisline.get_ydata()
        ind = event.ind
        print("X=" + str(np.take(xdata, ind)[0]))  # Print X point
        # print('Y='+str(np.take(ydata, ind)[0])) # Print Y point


def plot_datetime(
    year,
    month,
    day,
    t,
    y,
    colour="C0",
    estilo_linea="-",
    ancho_linea=1,
    transparencia=1,
):
    timestamps = array_datenums(year, month, day, t)
    t_graph = md.date2num(timestamps)

    plt.subplots_adjust(bottom=0.2)
    plt.xticks(rotation=25)
    ax = plt.gca()
    xfmt = md.DateFormatter("%H:%M:%S")
    ax.xaxis.set_major_formatter(xfmt)
    plt.plot(
        t_graph,
        y,
        linestyle=estilo_linea,
        color=colour,
        linewidth=ancho_linea,
        alpha=transparencia,
    )


def plot_rect(x, y):
    fig, current_ax = plt.subplots()
    plt.plot(x, y)
    plt.xlabel("Tiempo")
    plt.ylabel("|B| (nT)")
    plt.title("Seleccionar ")
    print("\n      click  -->  release")
    toggle_selector.RS = RectangleSelector(
        current_ax,
        line_select_callback,
        drawtype="box",
        useblit=True,
        button=[1, 3],  # don't use middle button
        minspanx=5,
        minspany=5,
        spancoords="pixels",
        interactive=True,
    )
    plt.connect("key_press_event", toggle_selector)


def plot_select(x, y, E):
    fig, ax = plt.subplots()
    ax.set_title("Seleccionar puntos", picker=True)
    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel("Flujo (cm⁻² sr⁻¹ s⁻¹)", picker=True)  # , bbox=dict(facecolor='red'))
    # line, = ax.semilogy(x, y,picker=5)
    for j in range(len(y[0, :])):
        ax.semilogy(x, y[:, j], label="{0:1.4g} eV".format(E[j]), picker=5)
    ax.set_ylim(1e4, 4 * 1e9)
    ax.legend()
    fig.canvas.mpl_connect("pick_event", onpick1)


def scatter_datetime(
    year, month, day, t, y, colour="C0", marcador="o", tamanio=1, transparencia=1,
):
    timestamps = array_datenums(year, month, day, t)
    t_graph = md.date2num(timestamps)

    plt.subplots_adjust(bottom=0.2)
    plt.xticks(rotation=25)
    ax = plt.gca()
    xfmt = md.DateFormatter("%H:%M:%S")
    ax.xaxis.set_major_formatter(xfmt)
    plt.scatter(
        t_graph,
        y,
        marker=marcador,
        color=colour,
        linewidths=tamanio,
        alpha=transparencia,
    )


def set_axes_equal(ax):  # creo que ya no anda: usar equal_axes
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    Input
      ax: a matplotlib axis, e.g., as output from plt.gca()."""

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def equal_axes(ax, x, y, z):
    """
    Le damos un ax y los datos y va a crear un cubo transparente que incluya todo.
    Entonces al tener una bounding box desaparece el problema.
    Hay que pasarle x,y,z de los datos más grandes (por ejemplo, la órbita entera
    en lugar del planeta, etc.)
    """
    max_range = np.array(
        [x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]
    ).max()
    Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (
        x.max() + x.min()
    )
    Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (
        y.max() + y.min()
    )
    Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (
        z.max() + z.min()
    )
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], "w")


def toggle_selector(event):
    print(" Key pressed.")
    if event.key in ["Q", "q"] and toggle_selector.RS.active:
        print(" RectangleSelector deactivated.")
        toggle_selector.RS.set_active(False)
        if event.key in ["A", "a"] and not toggle_selector.RS.active:
            print(" RectangleSelector activated.")
            toggle_selector.RS.set_active(True)


def zoom_scroll(ax, base_scale=0.2):
    def zoom_fun(event):
        current_x = ax.get_xlim()
        current_y = ax.get_ylim()
        current_xrange = (current_x[1] - current_x[0]) * 0.5
        current_yrange = (current_y[1] - current_y[0]) * 0.5
        xdata = event.xdata
        ydata = event.ydata
        if event.button == "up":
            scale_factor = 1 / base_scale
        elif event.button == "down":
            scale_factor = base_scale
        else:  # por si pasa algo raro
            scale_factor = 1
            print(event.button)
        # nuevos limites
        ax.set_xlim = [
            xdata - current_xrange * scale_factor,
            xdata + current_xrange * scale_factor,
        ]
        ax.set_ylim = [
            ydata - current_yrange * scale_factor,
            ydata + current_yrange * scale_factor,
        ]
        plt.draw()

    fig = ax.get_figure()
    fig.canvas.mpl_connect("scroll_event", zoom_fun)

    return zoom_fun
