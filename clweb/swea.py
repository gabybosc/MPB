import numpy as np
import os
import sys
import datetime as dt
from matplotlib.colors import LogNorm
import cdflib as cdf
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mpl
from pathlib import Path
from importar_datos import (
    find_path,
    tiempo_limite,
)
from cycler import cycler

sys.path.append("..")
from funciones import (
    find_nearest,
    datenum,
    donde,
    UTC_to_hdec,
    unix_to_decimal,
)


# ######### SWEA


def importar_swea_clweb(year, month, day, ti, tf):
    path = f"../../../datos/clweb/{year}-{month}-{day}/"

    swea = np.loadtxt(path + "SWEA.asc")

    """
    Los datos de SWEA se ven de la siguiente manera: para cada tiempo hay muchos
    valores de JE y de energía. Por eso sólo recorto el tiempo y dejo el JE entero
    y en tal caso en el script principal lo recorto.
    """

    t_swea, idx = np.unique(
        swea[:, 3] + swea[:, 4] / 60 + swea[:, 5] / 3600, return_index=True
    )  # hdec
    # al tomar "unique" estoy borrando los t que se repiten para las diferentes energias

    energias = [50 + i * 50 for i in range(3)]

    inicio = donde(t_swea, ti)
    fin = donde(t_swea, tf)  # este corte no sirve para JE ni para energy!!
    energy = swea[:, 7]
    JE_total = swea[:, -1]
    t_cut = t_swea[inicio:fin]

    JE_cut = np.zeros((len(t_cut), len(energias)))
    for i, energia in enumerate(energias):
        index = np.where(energy == find_nearest(energy, energia))[
            0
        ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
        JE = JE_total[index]
        JE_cut[:, i] = JE[inicio:fin]

    return swea, t_cut, JE_cut


def importar_swea_pds(year, month, day, ti, tf):
    date_orbit = dt.date(int(year), int(month), int(day))
    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")

    path = "../../../datos/SWEA/"
    # chequea que el archivo no está vacío
    if os.path.exists(path + f"mvn_swe_l2_svyspec_{year}{month}{day}.cdf"):
        if (
            Path(path + f"mvn_swe_l2_svyspec_{year}{month}{day}.cdf").stat().st_size
            > 1000
        ):
            swea = cdf.CDF(path + f"mvn_swe_l2_svyspec_{year}{month}{day}.cdf")

            flux_all = swea.varget("diff_en_fluxes")
            energia = swea.varget("energy")
            t_unix = swea.varget("time_unix")

            t = unix_to_decimal(t_unix)

            inicio = donde(t, ti)
            fin = donde(t, tf)

            t_cut = t[inicio:fin]

            flux = flux_all[inicio:fin]
            flux_plot = np.transpose(flux)[::-1]

    else:
        swea, t_cut, energia, flux_plot = 0, 0, 0, 0

    return swea, t_cut, energia, flux_plot


year, month, day, doy = 2016, "03", 16, 76
ti, tf = 17.85, 18.4


mpl.rcParams["axes.prop_cycle"] = cycler(
    "color", ["#5F021F", "#336699", "#C70039", "#00270F"]
)

swea_pds, t_pds, energia_pds, flux_plot_pds = importar_swea_pds(
    year, month, day, ti, tf
)
swea_cl, t_cl, JE_cl = importar_swea_clweb(year, month, day, ti, tf)

JE_pds = np.zeros((len(t_pds), 3))
energias = [50 + i * 50 for i in range(3)]

for i, energia in enumerate(energias):
    j = donde(energia_pds, energia)
    JE_pds[:, i] = flux_plot_pds[j]
plt.figure()
# plt.semilogy(t_cl, JE_cl)
plt.semilogy(t_pds, JE_pds)
plt.show()

plt.figure()
plt.imshow(
    flux_plot_pds,
    aspect="auto",
    origin="lower",
    cmap="inferno",
    extent=(t_pds[0], t_pds[-1], energia_pds[-1], energia_pds[0]),
    norm=LogNorm(vmin=1e4, vmax=1e9),
)
plt.show()
