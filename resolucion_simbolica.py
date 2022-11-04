import sympy as sp

k, h, omega, g, r2, r1 = sp.symbols('k h omega g r2 r1')

f = -sp.exp(k*h)*(k-omega**2/g)*(-(r2+r1) * omega**2 + g*k * (r2-r1)) + sp.exp(-k*h) * (r2-r1) * (omega**2-g*k) * (k + omega**2/g)

sp.solve(f, omega)
