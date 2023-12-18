from importar_datos import (
    importar_tiempos,
    importar_posiciones,
    importar_swia,
    importar_mag_1s,
)
from funciones import (
    long_inercial_iones,
    UTC_to_hdec,
    ancho_mpb,
    diezmar,
    giroradio,
    giroradio_termico,
)
from funciones_metodos import ajuste_conico
import matplotlib.pyplot as plt

"""
Quiero hacer análsis rápidos de caso para todos los que ya encontramos con sofi
Puedo por lo pronto calcular el c/wpi y comparar con las longs.
"""

path = "../../datos/bs_mpb/"
dates, times_mpb, times_bs = importar_tiempos(path)
p_mpb, p_min, p_max = importar_posiciones(path)

cwpi = []
ancho = []
gr = []
gth = []
fechas = []

for d in range(len(dates)):
    date = dates[d]
    t_mpb = [UTC_to_hdec(times_mpb[d, i]) for i in range(3)]
    t_bs = [UTC_to_hdec(times_bs[d, i]) for i in range(3)]
    pos_mpb = p_mpb[d, :]

    year, month, day = date[:4], date[5:7], date[8:]

    if t_mpb[0] > t_bs[0]:
        t_ms_i = t_bs[2]
        t_ms_f = t_mpb[0]
    else:
        t_ms_i = t_mpb[2]
        t_ms_f = t_bs[0]

    swia, t_swia, density, temperature, vel_mso = importar_swia(
        year, month, day, t_ms_i, t_ms_f
    )
    mag, t_mag, B, posicion = importar_mag_1s(year, month, day, t_ms_i, t_ms_f)

    if type(density) is not int:
        ion_length, ion_min, ion_max = long_inercial_iones(density)
        cwpi.append(ion_length)
        # para obtener el ancho de la MPB necesito la normal
        # puedo hacer el MVA en cada caso o hacer un fit rápido

        X1, Y1, Z1, L0, norm_vignes = ajuste_conico(pos_mpb)
        x14, x23 = ancho_mpb(
            t_mpb[0], t_mpb[1], t_mpb[1], t_mpb[2], norm_vignes, vel=2.5
        )
        ancho.append(x14)
        fechas.append(date)

        idx = diezmar(t_mag, t_swia)
        B_diezmado = B[idx]
        fg, rg, rgmin, rgmax = giroradio(B_diezmado, vel_mso)
        gr.append(rg)

        gth.append(giroradio_termico(B_diezmado, temperature))

plt.figure()
plt.plot(cwpi, ancho, ".")
plt.figure()
plt.plot(gr, ancho, ".")
plt.figure()
plt.plot(gth, ancho, ".")
plt.xlim(left=0, right=400)
plt.show()
