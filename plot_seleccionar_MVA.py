import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from matplotlib.colors import LogNorm
from funciones import Bpara_Bperp, fechas, tiempos, next_available_row
from funciones_plot import onpick1
from importar_datos import importar_mag_1s, importar_lpw, importar_swea, importar_swia
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gspread
from oauth2client.service_account import ServiceAccountCredentials


"""
Le paso un día entero y elijo la órbita y luego elijo los tiempos t1t2t3t4,
que va a guardar en un archivo aparte / en gdocs.
Para eso, grafica |B|, BxByBz, swea, swia y lpw.
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = tiempos()

mag, t, B, posicion = importar_mag_1s(year, month, day, ti, tf)
swea, t_swea, energia, flux_plot = importar_swea(year, month, day, ti, tf)
swia, t_swia, i_density, i_temp, vel_mso = importar_swia(year, month, day, ti, tf)
lpw, t_lpw, e_density = importar_lpw(year, month, day, ti, tf)

B_norm = np.linalg.norm(B, axis=1)

B_para, B_perp_norm, t_plot = Bpara_Bperp(
    B, t, ti + 0.2, tf - 0.2
)  # estos ti, tf tienen que ser menores que el total de datos

index = np.array((int(year), doy))


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
        plt.title("Spacebar when ready to click:")

        ax1 = plt.subplot2grid((3, 2), (0, 0))
        plt.plot(t_plot, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
        plt.plot(
            t_plot, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / B"
        )
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r"|$\Delta B$|/ B")
        ax1.grid()
        ax1.legend()

        ax4 = plt.subplot2grid((3, 2), (1, 0), sharex=ax1)
        ax4.plot(t, B)
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.set_ylabel("Bx, By, Bz (nT)")
        ax4.legend(["Bx", "By", "Bz"])
        ax4.grid()

        ax3 = plt.subplot2grid((3, 2), (2, 0), sharex=ax1)
        plt.plot(t, B_norm)
        ax3.grid()
        ax3.set_ylabel("|B| (nT)")
        ax3.set_xlabel("Tiempo (hdec)")

        ax5 = plt.subplot2grid((3, 2), (0, 1), sharex=ax1)
        ax5.set_ylabel("Energia", picker=True)  # , bbox=dict(facecolor='red'))
        plt.setp(ax5.get_xticklabels(), visible=False)
        im = plt.imshow(
            flux_plot,
            aspect="auto",
            origin="lower",
            extent=(t_swea[0], t_swea[-1], energia[-1], energia[0]),
            cmap="inferno",
            norm=LogNorm(vmin=1e4, vmax=1e9),
        )
        divider = make_axes_locatable(ax5)
        cax = divider.append_axes("top", size="7%", pad="1%")
        cb = plt.colorbar(im, cax=cax, orientation="horizontal")
        cax.xaxis.set_ticks_position("top")

        ax7 = plt.subplot2grid((3, 2), (1, 1), sharex=ax1)
        plt.setp(ax7.get_xticklabels(), visible=False)
        ax7.set_ylabel("Densidad de p+ \n del SW (cm⁻³)")
        plt.plot(t_swia, i_density)
        ax7.grid()

        ax6 = plt.subplot2grid((3, 2), (2, 1), sharex=ax1)
        ax6.set_ylabel("Densidad total \n de e- (cm⁻³)")
        ax6.set_xlabel("Tiempo (hdec)")
        plt.semilogy(t_lpw, e_density)
        ax6.grid()

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(
            fig.canvas, (ax1, ax3, ax4, ax5, ax6, ax7), color="black", lw=1
        )

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select MPB: ")
        val = np.asarray(plt.ginput(4))[:, 0]
        print("Selected values: ", val)
        outs = np.concatenate((index, val))
        outs = sorted(outs[2:6])

    print("Happy? Keyboard click for yes, mouse click for no.")
    happy = plt.waitforbuttonpress()

with open("outputs/t1t2t3t4.txt", "a") as file:
    file.write(f"{year}\t{doy}\t")
    for k in range(len(outs)):
        file.write(f"{outs[k]}\t")
    file.write("\n")

meses = {
    1: "Jan",
    2: "Feb",
    3: "Mar",
    4: "Apr",
    5: "May",
    6: "Jun",
    7: "Jul",
    8: "Aug",
    9: "Sep",
    10: "Oct",
    11: "Nov",
    12: "Dec",
}

scope = [
    "https://spreadsheets.google.com/feeds",
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive.file",
    "https://www.googleapis.com/auth/drive",
]

creds = ServiceAccountCredentials.from_json_keyfile_name("mpb_api.json", scope)

client = gspread.authorize(creds)

hoja_parametros = client.open("MPB").worksheet("Parametros")
hoja_MVA = client.open("MPB").worksheet("MVA")
hoja_Bootstrap = client.open("MPB").worksheet("Bootstrap")
hoja_Ajuste = client.open("MPB").worksheet("Ajuste")

nr = next_available_row(hoja_parametros)

mes = meses[int(month)]
hoja_parametros.update_acell(f"A{nr}", f"{int(day)} {mes} {int(year)}")
hoja_parametros.update_acell(f"B{nr}", f"{int(float(outs[0]))}")

cell_times = hoja_parametros.range(f"F{nr}:I{nr}")
for i, cell in enumerate(cell_times):
    cell.value = round(float(outs[i]), 4)
hoja_parametros.update_cells(cell_times)

hoja_MVA.update_acell(f"A{nr}", f"{int(day)} {mes} {int(year)}")
hoja_MVA.update_acell(f"B{nr}", f"{int(float(outs[0]))}")
hoja_Bootstrap.update_acell(f"A{nr}", f"{int(day)} {mes} {int(year)}")
hoja_Bootstrap.update_acell(f"B{nr}", f"{int(float(outs[0]))}")
hoja_Ajuste.update_acell(f"A{nr}", f"{int(day)} {mes} {int(year)}")
hoja_Ajuste.update_acell(f"B{nr}", f"{int(float(outs[0]))}")


plt.show(block=False)
