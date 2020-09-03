"""
Hace RK4 para integrar dy/dx = f(x, y, params)
y = vector desconocido en N-dim (va a ser el vector [by, bz, u])
x = posicion en la direccion entre x0 y xf
params = parametros (vector en k-dim)
y0 = array de N-dim de valores iniciales
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def diff(x, Y, eps, v00, mach):
    dy0 = (eps ** 2 - Y[2]) * Y[1] / (eps * Y[2])
    dy1 = -(v00 + (eps ** 2 - Y[2]) * Y[0]) / (eps * Y[2])
    dy2 = eps * v00 * Y[1] / (Y[2] - 1 / (mach ** 2 * Y[2]))
    return np.array([dy0, dy1, dy2])


# parámetros
eps = 0.05
v0 = 0.03
b0 = 15
v00 = v0 - b0 * eps ** 2
mach = 0.65
params = [eps, v00, mach]

# Integrates dY/dx=F(x,Y,par) using RK4
# Y: N-dim vector of unknowns (N=3)
# x: position in Mars-Sun direction x0 and xf
# params: k-dim array of parameters (k=3)
# Y0: N-dim array of initial values
# res=RK4(Y,Dydx,x,dx,Derivs)

# Resolver para MPB 1D estacionaria, ecs de 2 fluidos
x0 = 0
xf = 5  # va desde el bow shock hacia Marte
Nx = 29000
N = 3
Y = np.zeros(N)

# valores iniciales
Y[0] = b0  # by(0)
Y[1] = 0  # bz(0)
Y[2] = 1  # u(0)
dx = (xf - x0) / Nx
x = x0

sol = solve_ivp(diff, (x0, xf), Y, args=(eps, v00, mach), max_step=dx)
YY = sol.y
XX = sol.t

# presiones
p = YY[2] + 1 / (mach ** 2 * YY[2]) + eps ** 2 * (YY[0] ** 2 + YY[1] ** 2) / 2
# cte = 1 + 1 / mach ** 2 + eps ** 2 * b0 ** 2 / 2

pmag = eps ** 2 * (YY[0] ** 2 + YY[1] ** 2) / 2

plt.plot(pmag / p)
plt.show()

print(f"p_mag/p_total inicial = {pmag[0] / p[0]}")
print(f"p_mag/p_total final = {pmag[-1] / p[-1]}")
