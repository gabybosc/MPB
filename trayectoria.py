import numpy as np
import matplotlib.pyplot as plt
from importar_datos import importar_mag
from funciones import fechas, tiempos, donde, hdec_to_UTC

year, month, day, doy = fechas()
ti, tf = tiempos()

mag, t, B, posicion = importar_mag(year, month, day, ti, tf)

periapsis = donde(np.linalg.norm(posicion, axis=1), min(np.linalg.norm(posicion, axis=1)))
t_peri = t[periapsis]

print(f"el periapsis se alcanza a las {t_peri} (hdec) = {hdec_to_UTC(t_peri)} UTC")

plt.plot(t, posicion)

plt.figure()
plt.plot(t, np.linalg.norm(posicion, axis=1))
plt.plot(t, posicion[:, -1])
plt.show()

# year, month, day, hour, minute, second, millisecond, x,y,z
M = len(mag)
date = [np.ones(M) * int(year), np.ones(M) * int(month), np.ones(M) * int(day)]
date = np.transpose(np.array(date))
horario = mag[:, 2:6]
RM = posicion / 3390

tray = np.concatenate((date, horario, RM), axis=1)
np.savetxt('outputs/trajectory.txt', tray)
