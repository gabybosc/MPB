import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
from Titan.Analisis import importar_bepi
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import sys

sys.path.append("..")
from funciones_plot import onpick1
from funciones import Bpara_Bperp, next_available_row

"""
Le paso un día entero y elijo la órbita y luego elijo los tiempos t1t2t3t4,
que va a guardar en un archivo aparte / en gdocs.
Para eso, grafica |B|, BxByBz, swea, swia y lpw.
"""

np.set_printoptions(precision=4)

ti = 13.5
tf = ti + 1
t, B, posicion = importar_bepi(ti, tf)

B_norm = np.linalg.norm(B, axis=1)
B_para, B_perp_norm, t_plot = Bpara_Bperp(B, t, ti + 0.2, tf - 0.2)

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

        ax1 = plt.subplot2grid((3, 1), (0, 0))
        plt.plot(t_plot, B_para, linewidth=1, label=r"|$\Delta B \parallel$| / B")
        plt.plot(
            t_plot, B_perp_norm, "-.", linewidth=1, label=r"|$\Delta B \perp$| / B"
        )
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylabel(r"|$\Delta B$|/ B")
        ax1.grid()
        ax1.legend()

        ax4 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
        ax4.plot(t, B)
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.set_ylabel("Bx, By, Bz (nT)")
        ax4.legend(["Bx", "By", "Bz"])
        ax4.grid()

        ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
        plt.plot(t, B_norm)
        ax3.grid()
        ax3.set_ylabel("|B| (nT)")
        ax3.set_xlabel("Tiempo (hdec)")

        fig.canvas.mpl_connect("pick_event", onpick1)
        multi = MultiCursor(fig.canvas, (ax1, ax3, ax4), color="black", lw=1)

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

plt.show(block=False)

with open("t1t2t3t4.txt", "a") as file:
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

creds = ServiceAccountCredentials.from_json_keyfile_name("../mpb_api.json", scope)

client = gspread.authorize(creds)

hoja_parametros = client.open("MPB").worksheet("Parametros")
hoja_MVA = client.open("MPB").worksheet("MVA")
hoja_Bootstrap = client.open("MPB").worksheet("Bootstrap")
hoja_Ajuste = client.open("MPB").worksheet("Ajuste")

nr = next_available_row(hoja_parametros)

day, month, year = "10", "08", "2021"
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

# 18.9111
# 18.92222
