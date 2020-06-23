"""
Hace RK4 para integrar dy/dx = f(x, y, params)
y = vector desconocido en N-dim (va a ser el vector [u, by, bz])
x = posicion en la direccion entre x0 y xf
params = parametros (vector en k-dim)
y0 = array de N-dim de valores iniciales

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def diff(x, Y, eps=0.1, vy0=1, vz0=1):
    dy0 = eps * (vy0 * Y[2] - vz0 * Y[1]) / Y[0]
    dy1 = vz0 / (eps * Y[0]) + (eps / Y[0] - 1.0 / eps) * Y[2]
    dy2 = -vy0 / (eps * Y[0]) - (eps / Y[0] - 1.0 / eps) * Y[1]
    return np.array([dy0, dy1, dy2])


# parámetros
epsilon = 0.1
vy_0 = 1
vz_0 = 1
params = [epsilon, vy_0, vz_0]

# Integrates dY/dx=F(x,Y,par) using RK4
# Y: N-dim vector of unknowns (N=3)
# x: position in Mars-Sun direction x0 and xf
# params: k-dim array of parameters (k=3)
# Y0: N-dim array of initial values
# res=RK4(Y,Dydx,x,dx,Derivs)

# Resolver para MPB 1D estacionaria, ecs de 2 fluidos
x0 = 5
xf = 0  # va de marte al sol
Nx = 29000
N = 3
Y = np.zeros(N)

# valores iniciales
Y[0] = 1  # u(0)
Y[1] = 0  # by(0)
Y[2] = 0  # bz(0)
dx = (xf - x0) / Nx
x = x0

dydx = diff(x, Y)
sol = solve_ivp(diff, (x0, xf), Y)
YY = sol.y
XX = sol.t

plt.plot(YY[1, :], YY[2, :])
plt.xlabel("by")
plt.ylabel("bz")
plt.title("by vs bz, RK4 integration")
plt.show()
