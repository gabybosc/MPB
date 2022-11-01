import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.widgets import MultiCursor
import matplotlib as mpl
from cycler import cycler
import sys
from glob import glob

sys.path.append("..")
from funciones import datenum, donde, Bpara_Bperp, fechas, tiempos
from funciones_plot import onpick1

np.set_printoptions(precision=4)


def importar_VEX_mag_AMDA(year, month, day, ti, tf):
    path = glob(f"../../../datos/VEX MAG/{day}{month}{year}/*.txt")
    B = np.genfromtxt(path[0], usecols=[1, 2, 3])
    tt = np.genfromtxt(path[0], usecols=0, dtype="str")

    fecha = np.array([x.split("T") for x in tt])
    hora = np.array([x.split(":") for x in fecha[:, 1]])
    hh = np.array([int(x) for x in hora[:, 0]])
    mm = np.array([int(x) for x in hora[:, 1]])
    ss = np.array([float(x) for x in hora[:, 2]])
    t = hh + mm / 60 + ss / 3600  # hdec

    # hay que chequear por si hay NaN
    # me da los índices de los NaN sin repetir
    idx = np.unique(np.argwhere(np.isnan(B))[:, 0])
    B_del = np.delete(B, idx, axis=0)
    t_del = np.delete(t, idx, axis=0)

    inicio = donde(t, ti)
    fin = donde(t, tf)

    t_cut = t_del[inicio:fin]
    B_cut = B_del[inicio:fin]
    return t_cut, B_cut


# year, month, day, doy = fechas()
# ti, tf = tiempos()

year, month, day = 2007, 11, 21
ti, tf = 1, 3

index = np.array((int(year), int(day)))
t, B = importar_VEX_mag_AMDA(year, month, day, ti, tf)

Bnorm = np.linalg.norm(B, axis=1)

Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)

# ############ tiempos UTC


def tiempos_UTC(yy, mm, dd, t):
    tt = np.array([np.datetime64(datenum(yy, mm, dd, x)) for x in t])
    return tt


mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)


happy = False
while not happy:
    val = []
    while len(val) < 4:
        plt.clf()  # clear figure
        fig = plt.figure(
            1, constrained_layout=True
        )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
        fig.subplots_adjust(
            top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
        )
        # xfmt = md.DateFormatter("%H:%M")
        ax1 = plt.gca()

        plt.title("Spacebar when ready to click:")
        ax1 = plt.subplot2grid((3, 1), (0, 0))
        ax1.plot(t, Bnorm, linewidth=0.5)
        ax1.set_ylabel("|B| (nT)")
        ax1.set_title(f"VEX MAG {year}-{month}-{day}")
        ax1.grid()

        ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
        ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=0.5)
        ax2.plot(t, B[:, 1], label="By VSO", linewidth=0.5)
        ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=0.5)
        ax2.set_ylabel("B components (nT)")
        ax2.legend(loc="upper left")
        ax2.grid()

        ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
        ax3.plot(tpara, Bpara, linewidth=0.5, label="B ||")
        ax3.plot(tpara, Bperp, linewidth=0.5, label="B perp")
        ax3.set_ylabel("variación de Bpara perp")
        ax3.set_xlabel("Tiempo (hdec)")
        ax3.set_ylim([-0.1, 5])
        ax3.legend(loc="upper left")
        ax3.grid()

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(fig.canvas, (ax1, ax2, ax3), color="black", lw=1)

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select MPB: ")
        val = np.asarray(plt.ginput(4))[:, 0]
        print("Selected values: ", val)
        outs = sorted(val)

    print("Happy? Keyboard click for yes, mouse click for no.")
    happy = plt.waitforbuttonpress()

plt.show()

yy = int(year)
mm = int(month)
dd = int(day)
tiempo_mag = tiempos_UTC(yy, mm, dd, t)
tiempo_paraperp = tiempos_UTC(yy, mm, dd, tpara)
MPB = tiempos_UTC(yy, mm, dd, outs)

fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M")
for ax in [ax1, ax2, ax3]:
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xlim(tiempo_mag[0], tiempo_mag[-1])
    ax.grid()
    ax.legend(loc="upper left")

ax1 = plt.gca()

ax1 = plt.subplot2grid((3, 1), (0, 0))
ax1.plot(tiempo_mag, Bnorm, linewidth=0.5)
ax1.set_ylabel("|B| (nT)")
ax1.set_title(f"VEX MAG {year}-{month}-{day}")

ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax2.plot(tiempo_mag, B[:, 0], label="Bx VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 1], label="By VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 2], label="Bz VSO", linewidth=0.5)
ax2.set_ylabel("B components (nT)")

ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
ax3.plot(tiempo_paraperp, Bpara, linewidth=0.5, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_paraperp, Bperp, "-.", linewidth=0.5, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Relative variation \n of B")
ax3.set_xlabel("Tiempo (UTC)")
ax3.set_ylim([-0.1, 5])


for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
for ax in [ax1, ax2, ax3]:
    ax.axvspan(xmin=MPB[0], xmax=MPB[3], facecolor="#581845", alpha=0.5)
    ax.axvspan(xmin=MPB[1], xmax=MPB[2], facecolor="#FFC300", alpha=0.5)
    ax.grid()
    ax.legend(loc="upper right")
    # en un radio de 10 min de la MPB
    ax.set_xlim([MPB[0] - np.timedelta64(10, "m"), MPB[-1] + np.timedelta64(10, "m")])

# plt.tight_layout()

plt.show()

# buenas órbitas: SZA no tan alto, el campo en SW no es Bx
# 21 nov 2007
# 14 abr 2007
