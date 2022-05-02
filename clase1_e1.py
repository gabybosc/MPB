import numpy as np
import matplotlib.pyplot as plt

"""
Ejercicio 2 guía 1
"""

# como es numérico, uso las variables adimensionalizadas T/T0, x/L

alpha = 1


x = np.linspace(0, 1)  # adimensional

plt.figure()

for t in [0, 0.1, 0.25]:
    T = 1 - alpha * np.exp(-x) * np.sin(2 * np.pi * t)
    plt.plot(x, T, label=rf"t/$\tau$ = {t}")
plt.ylabel("T/T0")
plt.xlabel("x/L")
plt.legend()
plt.show()


"""
Ejercicio 3 guía 1
"""
#  streamlines d
X = np.linspace(-1, 1, 100)
Y = np.linspace(-1, 1, 100)

x, y = np.meshgrid(X, Y)
r = np.sqrt(x ** 2 + y ** 2)
theta = np.arctan(y / x)

Q = 1
u = 1
t = 1

l = 1 / np.tan(theta) + 2 * np.pi * r / Q * u * t / np.sin(theta)

plt.contour(x, y, l)
plt.show()

"""
Ejemplo 2D
"""
x = np.linspace(-10, 10, 100)
x0 = 1
y0 = 1
ux = 1
uy = 1


# ezcast pro
plt.figure()
for x0 in [1, 2, 3]:
    y = y0 + ux / uy * (x - x0)
    plt.plot(x, y, label=f"x0 = {x0}")
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.show()

"""
Ejercicio 4 guía 1
"""

x = np.linspace(0.1, 10, 100)
x0 = 1
y0 = 1
yt = 2
xt = 2
a = 1
b = 1
c = 1
t = 1

# una línea de cada tipo para diferentes a
for a in [1, 2]:
    tray = y0 + c / b * ((x / x0) ** (b / a) - 1)  # trayectoria
    stream = y0 + c / a * (1 + b * t) * np.log(x / x0)  # streamline, a t fijo
    traza = yt + c / b * (b * t + 1) * (1 - (xt / x) ** (b / a))  # traza, a t fijo

    plt.figure()
    plt.plot(x, tray, c="C0", label="trayectoria")
    plt.plot(x, stream, c="C1", label="streamline")
    plt.plot(x, traza, c="C2", label="traza")
    plt.title(f"visualizaciones para a = {a}b")
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.ylim((-5, 5))
plt.show()


#  visualizaciones por separado
plt.figure()
for x0 in range(30):
    tray = y0 + c / b * ((x / x0) ** (b / a) - 1)  # trayectoria
    plt.plot(x, tray, c="C0")
plt.title("trayectoria para diferentes x0")
plt.xlabel("x")
plt.ylabel("y")
plt.ylim((0, 5))
plt.show()

plt.figure()
for x0 in range(30):
    stream = y0 + c / a * (1 + b * t) * np.log(x / x0)  # streamline, a t fijo
    plt.plot(x, stream, c="C0")
plt.title("streamlines para diferentes x0")
plt.xlabel("x")
plt.ylabel("y")
plt.ylim((-5, 5))
plt.show()

plt.figure()
for xt in range(30):
    traza = yt + c / b * (b * t + 1) * (1 - (xt / x) ** (b / a))  # traza, a t fijo
    plt.plot(x, traza, c="C0")
plt.title("traza para diferentes xt")
plt.xlabel("x")
plt.ylabel("y")
plt.ylim((-5, 5))
plt.show()


#  superpuestas
plt.figure()
for x0 in range(30):
    tray = y0 + c / b * ((x / x0) ** (b / a) - 1)  # trayectoria
    stream = y0 + c / a * (1 + b * t) * np.log(x / x0)  # streamline, a t fijo
    plt.plot(x, tray, c="C0")
    plt.plot(x, stream, c="C1")
for xt in range(30):
    traza = yt + c / b * (b * t + 1) * (1 - (xt / x) ** (b / a))  # traza, a t fijo
    plt.plot(x, traza, c="C2")
plt.title("trayectoria, streamline y traza")
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.ylim((-5, 5))
plt.show()
