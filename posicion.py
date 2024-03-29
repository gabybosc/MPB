import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from funciones import UTC_to_hdec, datenum, donde
from importar_datos import importar_mag_1s

# from datetime import datetime

year = 2016
month = "03"
datos = np.loadtxt("outputs/para_sofi.txt")  # year doy hour X Y Z Rsd
datos_bs = np.genfromtxt(
    "outputs/BScrossings_March2016.txt", dtype="str", skip_header=3
)
dias = datos[:, 1] - 60
horas = datos[:, 2]

"""
year doy t_mpb x y z Rsd_mpb x_bs y_bs z_bs Rsd_bs t_bs
"""

dia_bs = datos_bs[:, 0].astype(int)
hora_bs = datos_bs[:, 1]
t_bs = np.array([UTC_to_hdec(hora_bs[i]) for i in range(len(hora_bs))])


mpb = np.transpose([dias, horas])
t_mpb = np.array(
    [np.datetime64(datenum(2016, 3, mpb[i, 0], mpb[i, 1])) for i in range(len(mpb))]
)

pos_mpb = datos[:, 3:6]
Rsd_mpb = datos[:, 6]
pos_bs = datos[:, 7:10]
Rsd_bs = datos[:, 10]
"""
busco la standoff distance de la MPB para graficar eso después
(y espero que Sofi haga lo propio con el BS)
"""


def ajuste_conico(posicion, x0=0.78, e=0.9, L=0.96):
    """me da la standoff distance en X y en Y según Vignes2000"""
    theta = np.linspace(0, np.pi * 3 / 4, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    R = posicion / 3390  # la posicion de la nave en RM
    # ###### Calculo mi propia elipse que pase por el punto.
    r0 = R - np.array([x0, 0, 0])
    theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

    L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))

    Rsd = x0 + L0 / (1 + e)
    # Rtd = np.sin(PHI) * L0 / (1 + e * np.cos(PHI))

    return Rsd  # , Rtd  # r1, L0, X1, Y1, Z1


# Rsd = []
# for i in range(len(pos_mpb)):
#     Rsd.append(ajuste_conico(pos_mpb[i]))
#
# # para guardar nuevas columnas:
# Rsd = np.reshape(Rsd, (datos.shape[0],1))   # adjust dimension of the new array
# result = np.append(datos, Rsd, 1)         # append as last column
# np.savetxt("outputs/para_sofi.txt", result, delimiter=" ", fmt="%s")

"""
el bs, uso los que tengamos mismo dia/hora
"""
pos_bs = []
for i, day in enumerate(dia_bs):
    ti = t_bs[i]
    mag, t, B, posicion = importar_mag_1s(year, month, day, ti, ti + 0.05)
    pos_bs.append(posicion[0, :])

bs = np.transpose([dia_bs, t_bs])
t = np.array(
    [np.datetime64(datenum(2016, 3, bs[i, 0], bs[i, 1])) for i in range(len(bs))]
)

idx = [donde(t, t_mpb[i]) for i in range(len(t_mpb))]
t_diezmado = t[idx]
bs_diezmado = np.array(pos_bs)[idx]

"""
[:4] year
[5:7] month
[8:10] day
[11:13] hour
[13:15] min
sirve para chequear que estemos en los mismos días/horas
"""
for i in range(len(t_diezmado)):
    if not t_diezmado[i].astype(str)[11:13] == t_mpb[i].astype(str)[11:13]:
        print(i, t_diezmado[i].astype(str)[11:19], t_mpb[i].astype(str)[11:19])

Rsd_bs = []
for i in range(len(bs_diezmado)):
    Rsd_bs.append(ajuste_conico(bs_diezmado[i], 0.64, 1.03, 2.04))


t_bs_dt = np.array(
    [
        np.datetime64(datenum(2016, 3, dias[i], t_diezmado[i]))
        for i in range(len(t_diezmado))
    ]
)
# para guardar nuevas columnas:
# Rsd_bs = np.reshape(Rsd_bs, (datos.shape[0],1))   # adjust dimension of the new array
# result = np.append(datos, Rsd_bs, 1)         # append as last column
# np.savetxt("outputs/para_sofi.txt", result, delimiter=" ", fmt="%s")
"""
Plots
"""

plt.clf()  # clear figure
fig = plt.figure(
    1
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.10, right=0.95, hspace=0.005, wspace=0.15
)
plt.xticks(rotation=25)
xfmt = md.DateFormatter("%H:%M:%S")

ax1 = plt.gca()
ax1.xaxis.set_major_formatter(xfmt)
ax1 = plt.subplot2grid((1, 1), (0, 0))
ax1.plot(t, Rsd_bs, ".")
ax1.plot(t_mpb, Rsd_mpb, ".")


plt.figure()
plt.scatter(Rsd_bs, Rsd_mpb)
plt.xlabel("Rsd BS")
plt.ylabel("Rsd MPB")
plt.show()
