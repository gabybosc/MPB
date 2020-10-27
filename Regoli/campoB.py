import numpy as np
import matplotlib.pyplot as plt

path = "../../../datos/simulacion_leonardo/"
campo = np.loadtxt(path + "campo.txt")
posicion = np.loadtxt(path + "pos_mhd.txt")

"""
Quiero plotear el campo y ver que me aparezca una MPB y un BS.
Primero, voy a buscar las posiciones y=0, z=0 y con eso plotear para esos índices
|B| vs x.
"""

"""
Por como están hechos los datos, cuando una de las posiciones es cero, todas las
variables se anulan. Voy a borrar esos datos entonces.
"""

ceros = [i for i in range(len(posicion)) if posicion[i, 1] == 0]

pos_cut = np.delete(posicion, ceros, axis=0)
B_cut = np.delete(campo, ceros, axis=0)

val = np.concatenate((pos_cut, B_cut), axis=1)

zcero = np.array([val[i, :] for i in range(len(val)) if np.abs(val[i, 2]) <= 0.1])
ycero = np.array([zcero[i, :] for i in range(len(zcero)) if np.abs(zcero[i, 1]) <= 0.1])
zona_interes = np.array(
    [ycero[i, :] for i in range(len(ycero)) if 0 < ycero[i, 0] <= 5]
)
# zcut es el cut final y el que voy a usar el resto del tiempo

x = zona_interes[:, 0]
B = np.linalg.norm(zona_interes[:, 3:6], axis=1)

# quiero reordenar los x de forma creciente

xB = [tuple(np.vstack((x, B))[:, i]) for i in range(len(x))]
xB_sort = sorted(xB, key=lambda f: f[0])

x_val = [x[0] for x in xB_sort]
y_val = [x[1] for x in xB_sort]


plt.scatter(x, B)
plt.plot(x_val, y_val, c="C1")
plt.show()

"""
Para hacer un remesh voy a: poner x,y,z + B y ordenar por x creciente, y creciente,
z creciente (en ese orden).
"""


#
# THETA = np.linspace(0, np.pi * 3 / 4, 100)
# PHI = np.linspace(0, 2 * np.pi, 100)
#
# L = 0.96
# e = 0.9
# x0 = 0.78
# r1 = L / (1 + e * np.cos(THETA))
#
# X1 = x0 + r1 * np.cos(THETA)
# Y1 = r1 * np.sin(THETA)
#
#
# plt.figure()
# plt.scatter(x[::100], y[::100], c=j[::100], cmap="Reds", alpha=0.5)
# plt.plot(X1, Y1, label="Mean MPB")
# plt.plot(X1, -Y1, c="C0")
# plt.clim(0, 0.1)
# plt.xlim(left=0)
# plt.colorbar(label="corriente (va de 0 a 0.5)")
# plt.xlabel("x (RM)")
# plt.ylabel("y (RM)")
# plt.legend()
# plt.title("Simulacion corrientes plano z=0")

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
