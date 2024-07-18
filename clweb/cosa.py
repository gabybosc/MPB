import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import welch

sys.path.append("..")
from funciones import fechas, t_clweb

"""
Hace un filtro butterworth para quitar el ruido de la señal que es de aproximadamente 180 ms.
Primero usa buttord para encontrar el orden.
Es un filtro digital con lo cual pide que las frecuencias estén normalizadas
respecto de la frec de Nyquist. En este caso Nyquist es 32Hz/2 = 16Hz.
Como las frecuencias están normalizadas, da lo mismo usar f o w.
N es el orden del filtro, en general voy a querer que N sea cercano a 10.

Comentado está otro filtro más fuerte.
"""

# year, month, day, doy = fechas()
year, month, day, doy = 2016, "03", 16, 76
hora = 18  # input("hora\n")

path = f"../../../datos/clweb/{year}-{month}-{day}/{hora}/"  # path a los datos desde la compu
mag = np.loadtxt(path + "MAG.asc")

t = t_clweb(mag)

M = np.size(t)  # el numero de datos

# el campo
B = np.zeros((M, 3))
B = mag[:, 6:9]

Bnorm = mag[:, -1]

# Input the provided data


# Calculate the power spectral density using Welch's method
fs = 1 / np.mean(np.diff(t))  # Sampling frequency
frequencies, power_spectral_density = welch(B[:, 0], fs, nperseg=8)

# Plot the power spectral density
plt.figure(figsize=(10, 6))
plt.loglog(frequencies, power_spectral_density)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power Spectral Density")
plt.title("Power Spectrum of the Magnetic Field")
plt.grid(True)
plt.show()
