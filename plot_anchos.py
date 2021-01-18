import matplotlib.pyplot as plt
from numpy import array

"""
plotea la long inercial, giroradio y el ancho de la MPB para cada cruce
plotea la js, fuerza sobre volumen y SZA.
"""

tiempos = ["24-11-17", "10-10-15", "12-10-15", "05-04-16", "31-03-16", "16-03-16"]
hmin = [115, 39, 18.5, 44.4, 38.5, 82.2]
hmax = [447, 97.1, 73, 175, 122, 174]
long_inercial = [120, 159, 133, 130, 101, 98]
giroradio_termico = [201, 470, 160, 266, 210, 145]
larmor = [45.9, 203, 61.5, 168, 75.9, 68.4]
larmor_normal = [49.8, 319, 68, 20.2, 28.9, 69.1]


js = [
    [1.232, -7.469, -7.503],
    [10.158, 36.428, 10.138],
    [-3.657, 15.663, -9.491],
    [-1.725, -4.477, -15.406],
    [6.242, 8.683, 11.105],
    [-2.859, -19.188, -12.640],
]
normales = [
    [0.981, -0.032, 0.193],
    [0.956, -0.286, 0.070],
    [0.956, 0.048, -0.290],
    [0.815, -0.575, 0.076],
    [0.871, -0.476, -0.117],
    [0.920, -0.302, 0.251],
]
js_norm = [10.7, 39.200, 18.700, 16.100, 15.400, 23.200]

jv = [92.8, 403.224, 255.832, 363.000, 401.000, 282.000]

F = [2.39e-15, 4.37e-14, 3.38e-14, 9.21e-15, 1.34e-14, 1.20e-14]

SZA = [12.10, 17.60, -21.60, 27.40, 23.00, 25.60]

plt.figure(0)
plt.plot(tiempos, hmin, ".", color="C0")
plt.plot(tiempos, hmax, ".", color="C0")
for i in range(len(tiempos)):
    plt.vlines(
        x=tiempos[i], ymin=hmin[i], ymax=hmax[i], zorder=0, linewidth=2, color="C0"
    )
plt.scatter(tiempos, long_inercial, zorder=1, label="long_inercial", color="C1")
plt.scatter(tiempos, giroradio_termico, zorder=1, label="giroradio térmico", color="C2")
plt.scatter(tiempos, larmor, zorder=1, label="Larmor", color="C3")
plt.scatter(
    tiempos, larmor_normal, zorder=1, label="Larmor proyección normal", color="C4"
)
plt.legend()
plt.xlabel("Fecha (dd-mm-aa)")
plt.ylabel("Longitud (km)")
plt.title(
    "Comparación de la longitud inercial y el giroradio de los protones \n con el ancho de la MPB"
)

plt.xlabel("Fecha (dd-mm-aa)")
plt.ylabel("Corriente sup (mA/m2)")


fig = plt.figure(
    1, figsize=(20, 10)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)

ax1 = plt.subplot2grid((3, 1), (0, 0))
ax1.scatter(tiempos, js_norm, marker="s", color="C0")
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("Corriente sup (mA/m2)")

ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax2.scatter(tiempos, array(F) * 1e15, marker="o", color="C1")
ax2.set_ylabel("Fuerza /vol")

ax2 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
ax2.scatter(tiempos, SZA, marker="d", color="C2")
ax2.set_ylabel("SZA")

ax2.set_xlabel("Fecha (dd-mm-aa)")

plt.show()
