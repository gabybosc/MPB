import matplotlib.pyplot as plt
from numpy import array

tiempos = ["24-11-17", "10-10-15", "12-10-15", "05-04-16", "31-03-16", "16-03-16"]
js = [[1.232,-7.469,-7.503],[10.158,36.428,	10.138], [-3.657,15.663,-9.491],
        [-1.725,-4.477,-15.406], [6.242,8.683,11.105], [-2.859,-19.188,-12.640]]
normales = [[0.981,-0.032,0.193], [0.956,-0.286,0.070], [0.956,0.048,-0.290],
        [0.815,-0.575,0.076], [0.871,-0.476,-0.117], [0.920,-0.302,0.251]]
js_norm = [10.7, 39.200,
18.700,
16.100,
15.400,
23.200]

jv = [92.8, 403.224,
255.832,
363.000,
401.000,
282.000]

F = [2.39E-15, 4.37E-14,
3.38E-14,
9.21E-15,
1.34E-14,
1.20E-14]

SZA = [12.10,
17.60,
-21.60,
27.40,
23.00,
25.60]

# plt.scatter(tiempos, jv, marker='o', color="C1", label='corriente en vol')
# plt.scatter(tiempos, F, marker='d', color="C2", label='Fuerza por unidad de vol')

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
ax1.scatter(tiempos, js_norm, marker='s', color="C0")
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("Corriente sup (mA/m2)")

ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax2.scatter(tiempos, array(F)*1E15, marker='o', color="C1")
ax2.set_ylabel("Fuerza /vol")

ax2 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
ax2.scatter(tiempos, SZA, marker='d', color="C2")
ax2.set_ylabel("SZA")

ax2.set_xlabel("Fecha (dd-mm-aa)")

plt.show()
