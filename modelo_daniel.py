"""
Hace RK4 para integrar dy/dx = f(x, y, params)
y = vector desconocido en N-dim
x = posicion en la direccion entre x0 y xf
params = parametros (vector en k-dim)
y0 = array de N-dim de valores iniciales

"""

import numpy as np

# Resolver para MPB 1D estacionaria, ecs de 2 fluidos
x0 = 5
xf = 0  # va de marte al sol
Nx = 29000
N = 3
y = np.zeros(N)

# parámetros
epsilon = 0.1
uy_0 = 1
uz_0 = 1
params = [epsilon, uy_0, uz_0]

# valores iniciales
y[0] = 1
y[1] = 0
y[2] = 0

# resolver la ec dif
dy_0 = epsilon * (uy_0 * y[2] - uz_0 * y[1]) / y[0]  # creo que y es el campo B


# Integrar con RK4
dx = (xf - x0) / Nx
YY = np.zeros((N, Nx))
YY[:, 0] = y
XX = np.zeros(Nx)
XX[0] = x0
x = x0

# empzamos el loop
for i in range(Nx):
    x = x + dx
    XX[i] = x
    dydx = np.diff(x, Y)  # buscar como es la sintaxis de esto
    Y = RK4(Y, dydx, x, dx)  # hacer la funcion RK4
    # if Y(0) lt 0. then goto,out
    # YY(*,ix)=Y ???
