import numpy as np
import matplotlib.pyplot as plt

path = "../../datos/simulacion_leonardo/"
corrientes = np.loadtxt(path + "corriente.txt")
campo = np.loadtxt(path + "campo.txt")
posicion = np.loadtxt(path + "pos_mhd.txt")

"""
Usamos un poco de list comprehension para elegir las posiciones menores a 2Rm en
valor absoluto. Creo que en realidad para x me alcanza con x entre [0,2],
y entre [-2,2] y z entre [-2,2].
A su vez, a medida que recorto x,y,z tengo que ir recortando también la corriente.
Si quiero recortar más datos tengo que concatenarlos en val.
Ahora querría agregarle Marte y la posición de la MPB de vignes.
A su vez, para tener un corte en los planos x=0, y=0, z=0 tengo que agregar esa
constraint.
"""

val = np.concatenate((posicion, corrientes), axis=1)

xcut = np.array(
    [val[i, :] for i in range(len(val)) if val[i, 0] <= 2 and val[i, 0] >= 0]
)
ycut = np.array([xcut[i, :] for i in range(len(xcut)) if np.abs(xcut[i, 1]) <= 2])
zcut = np.array([ycut[i, :] for i in range(len(ycut)) if np.abs(ycut[i, 2]) <= 2])
# zcut es el cut final y el que voy a usar el resto del tiempo

zcero = np.array([zcut[i, :] for i in range(len(zcut)) if np.abs(zcut[i, 2]) <= 0.01])

x = zcero[:, 0]
y = zcero[:, 1]
z = zcero[:, 2]
j = np.linalg.norm(zcero[:, 3:6], axis=1)

THETA = np.linspace(0, np.pi * 3 / 4, 100)
PHI = np.linspace(0, 2 * np.pi, 100)

L = 0.96
e = 0.9
x0 = 0.78
r1 = L / (1 + e * np.cos(THETA))

X1 = x0 + r1 * np.cos(THETA)
Y1 = r1 * np.sin(THETA)


plt.figure()
plt.scatter(x[::100], y[::100], c=j[::100], cmap="Reds", alpha=0.5)
plt.plot(X1, Y1, label="Mean MPB")
plt.plot(X1, -Y1, c="C0")
plt.clim(0, 0.1)
plt.xlim(left=0)
plt.colorbar(label="corriente (va de 0 a 0.5)")
plt.xlabel("x (RM)")
plt.ylabel("y (RM)")
plt.legend()
plt.title("Simulacion corrientes plano z=0")

# plt.figure()
# plt.scatter(x[::100], z[::100], c=j[::100], cmap="Blues", alpha=0.5)
# plt.plot(X1, Y1, c="C1", label="Mean MPB")
# plt.plot(X1, -Y1, c="C1")
# plt.clim(0, 0.1)
# plt.xlim(left=0)
# plt.colorbar(label="corriente (va de 0 a 0.5)")
# plt.xlabel("x (RM)")
# plt.ylabel("z (RM)")
# plt.legend()
# plt.title("Simulacion corrientes plano (x,z)")
#
# plt.figure()
# plt.scatter(x[::100], y[::100], c=j[::100], cmap="Greens", alpha=0.5)
# plt.plot(X1, Y1, label="Mean MPB")
# plt.plot(X1, -Y1, c="C0")
# plt.clim(0, 0.1)
# plt.xlim(left=0)
# plt.colorbar(label="corriente (va de 0 a 0.5)")
# plt.xlabel("y (RM)")
# plt.ylabel("z (RM)")
# plt.legend()
# plt.title("Simulacion corrientes plano (y,z)")

plt.show()
# adding the Contour lines with labels
# cset = contour(Z,arange(-1,1.5,0.2),linewidths=2,cmap=cm.Set2)
# clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
# colorbar(im) # adding the colobar on the right
# # latex fashion title
# title('$z=(1-x^2+y^3) e^{-(x^2+y^2)/2}$')
