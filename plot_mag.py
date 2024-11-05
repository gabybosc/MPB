import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag_1s
import datetime as dt

"""
Plotea sólo los datos de MAG de baja resolución
"""

np.set_printoptions(precision=4)

year = 2014
for doy in range(352, 366):
    print(doy)
    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)

    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    mag, t, B, posicion = importar_mag_1s(year, month, day, 0.1, 23.9)

    B_norm = np.linalg.norm(B, axis=1)

    fig = plt.figure(
        2, figsize=(8, 30)
    )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
    fig.subplots_adjust(
        top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
    )
    plt.xticks(rotation=25)

    ax1 = plt.subplot2grid((4, 1), (0, 0))
    ax1.plot(t, B_norm, c="k")
    ax1.set_ylabel("|B| (nT)")
    ax1.set_title(f"MAVEN MAG {year}-{month}-{day}")
    # ax1.set_ylim([3, 8])
    # ax1.set_xlim([4, 4.2])
    ax1.grid()

    ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
    ax2.plot(t, B[:, 0], c="C0", label="Bx MSO")
    ax2.set_ylabel("Bx (nT)")

    ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
    ax3.set_ylabel("By (nT)")
    ax3.plot(t, B[:, 1], c="C1", label="By MSO")

    ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)
    ax4.set_ylabel("Bz (nT)")
    ax4.plot(t, B[:, 2], c="C2", label="Bz MSO")
    ax4.set_xlabel("Tiempo (hdec)")

    for ax in [ax2, ax3, ax4]:
        ax.set_ylim([-10, 10])
        ax.grid()

    plt.show()
