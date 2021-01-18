import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from funciones import fechas, donde
from funciones_plot import onpick1
from importar_datos import importar_mag_1s, importar_swia, importar_fila

np.set_printoptions(precision=4)

"""
Ex densidades.py
Calcula el giroradio térmico y la long inercial. Hay que seleccionar el rango
en el cual tomo la densidad (el intervalo upstream es elegido a mano)
"""

# #########CONSTANTES
mp = 1.67e-27  # masa del proton en kg
# mp = 1.5e-10 #masa del proton en joules/c^2
kB = 1.38e-23  # cte de Boltzmann en J/K
q_e = 1.602e-19  # carga del electrón en C

# ########## DATOS
year, month, day, doy = fechas()
ti = int(input("Hora del cruce (HH)\n"))
tf = ti + 1
in_out = input("Inbound? (y/n)\n")

fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
    year, month, day, doy, ti
)

if in_out == "y":
    t1 = float(hoja_parametros.cell(fila, 6).value)
    t2 = float(hoja_parametros.cell(fila, 7).value)
    t3 = float(hoja_parametros.cell(fila, 8).value)
    t4 = float(hoja_parametros.cell(fila, 9).value)

else:
    t1 = float(hoja_parametros.cell(fila, 9).value)
    t2 = float(hoja_parametros.cell(fila, 8).value)
    t3 = float(hoja_parametros.cell(fila, 7).value)
    t4 = float(hoja_parametros.cell(fila, 6).value)


mag, t_mag, B, posicion = importar_mag_1s(year, month, day, ti, tf)
swia, t_swia, density, temperature, vel_mso_xyz = importar_swia(
    year, month, day, ti, tf
)

# quiero diezmar el tiempo y el campo para que sean los mismos que tiene swia
idx = np.zeros(len(t_swia))
for i in range(len(idx)):
    idx[i] = donde(t_mag, t_swia[i])
idx = idx.astype(int)

t_diezmado = t_mag[idx]  # lo diezmó

B_cut = B[idx]
posicion_cut = posicion[idx]


happy = False
while not happy:
    val = []
    while len(val) < 2:
        plt.clf()

        fig = plt.figure(1, constrained_layout=True)
        fig.subplots_adjust(
            top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
        )
        plt.title("Click to select upstream region:")

        ax1 = plt.subplot2grid((1, 1), (0, 0))
        for xc in [t1, t2, t3, t4]:
            plt.axvline(x=xc, color="g", linestyle="-.", linewidth=1)
        ax1.grid()
        ax1.legend()
        plt.plot(t_swia, density, label="density")
        ax1.set_ylabel(r"ion density (cm⁻³)")
        ax1.set_xlabel("time (hdec)")

        fig.canvas.mpl_connect("pick_event", onpick1)
        cursor = Cursor(ax1, useblit=True, color="k", linewidth=1)

        zoom_ok = False
        print("\nSpacebar when ready to click:\n")
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress(-1)
        print("Click to select upstream region: ")
        val = np.asarray(plt.ginput(2))[:, 0]
        print("Selected values: ", val)
        outs = sorted(val)

    print("Happy? Keyboard click for yes, mouse click for no.")
    happy = plt.waitforbuttonpress()

plt.show()


ti_swia = donde(t_swia, outs[0])
tf_swia = donde(t_swia, outs[1])

ti_mag = donde(t_diezmado, outs[0])
tf_mag = donde(t_diezmado, outs[1])

inicio_mpb = donde(t_diezmado, t1)
fin_mpb = donde(t_diezmado, t4)

ancho_mpb = np.linalg.norm(posicion_cut[inicio_mpb] - posicion_cut[fin_mpb])  # km

####################

density_mean = np.empty(tf_swia - ti_swia)  # upstream
paso = 20  # cada paso son 4 segundos.
for i in range(tf_swia - ti_swia):
    density_mean[i] = np.mean(
        density[ti_swia + i - paso : ti_swia + i]
    )  # toma desde atrás del ti así no se mete en la MPB nunca

# Elijo la densidad promedio media
idx_d = int(np.abs(np.argmin(density_mean) + np.argmax(density_mean)) / 2)
density_avg = density_mean[idx_d]  # cm

ion_length = 2.28e07 / np.sqrt(density_avg) * 1e-5  # km
print(f"Ion inertial length = {ion_length:1.3g} km")


B_avg = np.empty((len(idx), 3))
B_avg_normalized = np.empty((len(idx), 3))
temp_para_xyz = np.empty((len(idx), 3))

for i in range(len(idx) - 1):
    B_avg[i, :] = np.mean(B_cut[i : i + 30, :], axis=0) * 1e-5  # lo paso a gauss
    B_avg_normalized[i, :] = B_avg[i, :] / np.linalg.norm(B_avg[i, :])  # adimensional
    temp_para_xyz[i, :] = (
        np.dot(B_avg_normalized[i, :], temperature[i, :]) * B_avg_normalized[i, :]
    )  # eV

temp_para = np.linalg.norm(temp_para_xyz, axis=1)  # eV
temp_perp = np.linalg.norm(temperature - temp_para_xyz, axis=1)  # eV

thermal_gyroradius = np.empty(tf_mag - ti_mag)
for i in range(tf_mag - ti_mag):
    thermal_gyroradius[i] = (
        1.02e02
        * np.sqrt(temp_perp[ti_mag + i])
        / np.linalg.norm(B_avg[ti_mag + i, :])
        * 1e-5
    )  # km

print(
    f"thermal ion gyroradius mean = {np.nanmean(thermal_gyroradius, axis=0):1.3g} km"
)  # nanmean ignora los nans

plt.figure()
for xc in [t1, t2, t3, t4]:
    plt.axvline(x=xc, color="black", linewidth=0.5)
plt.plot(t_swia, density, label="density")
plt.axvline(x=t_diezmado[ti_mag], color="red")
plt.axvline(x=t_diezmado[tf_mag], color="red")
plt.axvline(x=t_swia[ti_swia], color="blue")
plt.axvline(x=t_swia[tf_swia], color="blue")
plt.scatter(
    t_swia[ti_swia - int(paso / 2) : tf_swia - int(paso / 2)],
    density_mean,
    color="r",
    label="density avg",
)
plt.ylabel("ion density (cm⁻³)")
plt.xlabel("time (hdec)")
plt.legend()
# plt.show(block=False)
plt.show()
