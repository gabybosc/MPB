import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.widgets import MultiCursor
import matplotlib as mpl
from cycler import cycler
import sys
from _importar_datos import importar_MAG, importar_ELS_clweb
from _old_fit_venus import plot_orbita

sys.path.append("..")
from funciones import datenum, donde, Bpara_Bperp, fechas, tiempos, find_nearest
from funciones_plot import onpick1

np.set_printoptions(precision=4)

"""
comparar con los datos de 1Hz a ver si están bien calibrados
puedo o hacer un avg o un downsampling
"""

year, month, day, doy = fechas()

# plt.figure()
# plt.plot(t, B)
# plt.show()

ti, tf = 0, 24  # tiempos()
t, B, pos, cl, tpos = importar_MAG(year, doy, ti, tf)
if cl == True:
    Bpara, Bperp, tpara = Bpara_Bperp(B, t, ti, tf)  # si son datos de clweb 1s
else:
    # para datos de PDS filtrados y diezmados
    Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)
Bnorm = np.linalg.norm(B, axis=1)
pos_RV = pos / 6050

t_els, ELS = importar_ELS_clweb(year, doy, ti, tf)
energy = ELS[:, 7]
JE_total = ELS[:, -1]
energias = [25 + i * 25 for i in range(4)]


# ############ tiempos UTC


def tiempos_UTC(yy, mm, dd, t):
    tt = np.array([np.datetime64(datenum(yy, mm, dd, x)) for x in t])
    return tt


def altitude(SZA):
    alt = 0.11 * SZA**2 - 0.22 * SZA + 389
    return alt / 6050


def fit():
    sza = np.linspace(0, np.pi, 100)
    alt = 1 + altitude(sza * 180 / np.pi)

    y_alt = np.array([alt[i] * np.sin(sza[i]) for i in range(len(alt))])
    x_alt = np.array([alt[i] * np.cos(sza[i]) for i in range(len(alt))])

    yz = y_alt[x_alt >= 0]
    xx = x_alt[x_alt >= 0]

    return xx, yz


mpl.rcParams["axes.prop_cycle"] = cycler(
    "color",
    ["#003f5c", "#ffa600", "#de425b", "#68abb8", "#f3babc", "#6cc08b", "#cacaca"],
)

xx, yz = fit()
orbita = np.sqrt(pos_RV[:, 1] ** 2 + pos_RV[:, 2] ** 2)

plot_orbita(pos_RV, orbita, xx, yz)
plt.show()
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
        ax1 = plt.subplot2grid((2, 2), (0, 0))
        ax1.plot(t, Bnorm, linewidth=0.5)
        ax1.set_ylabel("|B| (nT)")
        ax1.set_title(f"VEX MAG {year}-{month}-{day}")
        ax1.grid()

        ax2 = plt.subplot2grid((2, 2), (1, 0), sharex=ax1)
        ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=0.5)
        ax2.plot(t, B[:, 1], label="By VSO", linewidth=0.5)
        ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=0.5)
        ax2.set_ylabel("B components (nT)")
        ax2.legend(loc="upper left")
        ax2.grid()

        ax3 = plt.subplot2grid((2, 2), (0, 1), sharex=ax1)
        ax3.plot(tpara, Bpara, linewidth=0.5, label="B ||")
        ax3.plot(tpara, Bperp, linewidth=0.5, label="B perp")
        ax3.set_ylabel("variación de Bpara perp")
        ax3.set_xlabel("Tiempo (hdec)")
        ax3.set_ylim([-0.1, 5])
        ax3.legend(loc="upper left")
        ax3.grid()

        ax4 = plt.subplot2grid((2, 2), (1, 1), sharex=ax1)
        for energia in energias:
            index = np.where(energy == find_nearest(energy, energia))[
                0
            ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
            plt.semilogy(
                t_els[index],
                JE_total[index],
                label=f"{energia} eV",
                linewidth=0.5,
            )
        ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")
        ax4.legend(loc="center right")
        ax4.grid()

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(fig.canvas, [ax1, ax2, ax3, ax4], color="black", lw=1)

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
tiempo_els = tiempos_UTC(yy, mm, dd, t_els)
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

ax1 = plt.subplot2grid((4, 1), (0, 0))
ax1.plot(tiempo_mag, Bnorm, linewidth=0.5)
ax1.set_ylabel("|B| (nT)")
ax1.set_title(f"VEX MAG {year}-{month}-{day}")

ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
ax2.plot(tiempo_mag, B[:, 0], label="Bx VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 1], label="By VSO", linewidth=0.5)
ax2.plot(tiempo_mag, B[:, 2], label="Bz VSO", linewidth=0.5)
ax2.set_ylabel("B components (nT)")

ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
ax3.plot(tiempo_paraperp, Bpara, linewidth=0.5, label=r"|$\Delta B \parallel$| / |B|")
ax3.plot(tiempo_paraperp, Bperp, "-.", linewidth=0.5, label=r"|$\Delta B \perp$| / |B|")
ax3.set_ylabel("Relative variation \n of B")
ax3.set_xlabel("Tiempo (UTC)")
ax3.set_ylim([-0.1, 5])

ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[
        0
    ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
    plt.semilogy(
        tiempo_els[index],
        JE_total[index],
        label=f"{energia} eV",
        linewidth=0.5,
    )
ax4.legend(loc="center right")
ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")

for ax in [ax1, ax2, ax3]:
    plt.setp(ax.get_xticklabels(), visible=False)
for ax in [ax1, ax2, ax3, ax4]:
    ax.axvspan(xmin=MPB[1], xmax=MPB[2], facecolor="#79B953", alpha=0.5)
    ax.axvspan(xmin=MPB[0], xmax=MPB[1], facecolor="#cdcdcd", alpha=0.7)
    ax.axvspan(xmin=MPB[2], xmax=MPB[3], facecolor="#cdcdcd", alpha=0.7)
    ax.grid()
    ax.legend(loc="upper right")
    # en un radio de 10 min de la MPB
    ax.set_xlim([MPB[0] - np.timedelta64(10, "m"), MPB[-1] + np.timedelta64(10, "m")])

# plt.tight_layout()

plt.show()

# buenas órbitas: SZA no tan alto, el campo en SW no es Bx
# 21 nov 2007
# 14 abr 2007
with open("../outputs/VEX_times.txt", "a") as file:
    file.write(f"{year}\t{doy}\t")
    for k in range(len(outs)):
        file.write(f"{outs[k]}\t")
    file.write("\n")
