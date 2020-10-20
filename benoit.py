import numpy as np
import matplotlib.pyplot as plt
from funciones import fechas

"""
Plotea el campo de MAVEN vs el campo cortical predicho por el modelo de
Benoit Langlais. Necesito sus datos para poder correrlo, por ahora lo tengo
para las seis órbitas del grl2020.
"""

year, month, day, doy = fechas()
data = np.genfromtxt(
    f"../../datos/benoit/{year}{doy}.syn",
    invalid_raise=False,
    missing_values="",
    usemask=False,
    filling_values=0.0,
)
cain = np.loadtxt(f"../../datos/clweb/{year}-{month}-{day}/cain.asc")

B_cain = cain[:, -1]
B_simulation = data[:, 4:7]
B_measured = data[:, 7:10]

B_s_norm = np.linalg.norm(B_simulation, axis=1)
B_m_norm = np.linalg.norm(B_measured, axis=1)

t_cain = cain[:, 3] + cain[:, 4] / 60 + cain[:, 5] / 3600
t_dec = data[:, 0] - int(data[0, 0])
t = t_dec * 24

plt.plot(t, B_m_norm, label="MAVEN")
plt.plot(t, B_s_norm, label="Langlais")
plt.plot(t_cain, B_cain, label="Cain")
plt.xlabel("t (hdec)")
plt.ylabel("B")
plt.title(f"Cain model vs Langlais model vs MAVEN\n{day}-{month}-{year}")
plt.xlim(5, 6)
plt.ylim(-10, 100)
plt.legend()
plt.show()
