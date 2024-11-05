<<<<<<< HEAD
import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag
import sys

sys.path.append("..")

from funciones import fechas

"""
Plotea sólo B para todo el día
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = 0.01, 24

# path = f"../../../datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la laptop

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

Bnorm = np.linalg.norm(B, axis=1)

# plt.plot(t, Bnorm)
# plt.ylabel("|B| (nT)")
# plt.xlabel("Time (hdec)")

fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)

ax1 = plt.subplot2grid((2, 1), (1, 0))
ax1.plot(t, B[:, 0], label="Bx MSO")
ax1.plot(t, B[:, 1], label="By MSO")
ax1.plot(t, B[:, 2], label="Bz MSO")
ax1.set_ylabel("B components (nT)")

ax2 = plt.subplot2grid((2, 1), (0, 0), sharex=ax1)
plt.plot(t, Bnorm)
plt.ylabel("|B| (nT)")
ax2.set_ylim([0, 10])
ax2.set_title(f"MAVEN MAG SWEA SWIA {year}-{month}-{day}")
plt.show()

# 2015 09 22
=======
import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag
import sys

sys.path.append("..")

from funciones import fechas

"""
Plotea sólo B para todo el día
"""

np.set_printoptions(precision=4)

year, month, day, doy = fechas()
ti, tf = 0.01, 24

# path = f"../../../datos/clweb/{year}-{month}-{day}/"  # path a los datos desde la laptop

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

Bnorm = np.linalg.norm(B, axis=1)

# plt.plot(t, Bnorm)
# plt.ylabel("|B| (nT)")
# plt.xlabel("Time (hdec)")

fig = plt.figure(
    2, figsize=(8, 30)
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.0, wspace=0.15
)
plt.xticks(rotation=25)

ax1 = plt.subplot2grid((2, 1), (1, 0))
ax1.plot(t, B[:, 0], label="Bx MSO")
ax1.plot(t, B[:, 1], label="By MSO")
ax1.plot(t, B[:, 2], label="Bz MSO")
ax1.set_ylabel("B components (nT)")

ax2 = plt.subplot2grid((2, 1), (0, 0), sharex=ax1)
plt.plot(t, Bnorm)
plt.ylabel("|B| (nT)")
ax2.set_ylim([0, 10])
ax2.set_title(f"MAVEN MAG SWEA SWIA {year}-{month}-{day}")
plt.show()

# 2015 09 22
>>>>>>> 89a36a6263086cb4f9fe8802bd1da03877917733
