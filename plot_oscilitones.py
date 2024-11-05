import numpy as np
import matplotlib.pyplot as plt
from funciones import t_clweb, find_nearest, giroradio, donde, Mij
from funciones_plot import hodograma


# Código de importación y otras funciones


def importar():
    path = "../../datos/oscilitones/"

    swia = np.loadtxt(path + "swifa.asc")
    tswia = t_clweb(swia)
    vel_swia = swia[:, 6:9]
    vel_norm_swia = swia[:, -1]

    swea = np.loadtxt(path + "swea.asc")
    tswea = t_clweb(swea)
    vel_swea = swea[:, 6:9]
    vel_norm_swea = swea[:, -1]

    mag = np.loadtxt(path + "mag.asc")
    tmag = t_clweb(mag)
    B = mag[:, 6:9]

    return tmag, tswia, tswea, B, vel_norm_swia, vel_norm_swea, vel_swia, vel_swea


tmag, tswia, tswea, B, vel_norm_swia, vel_norm_swea, vel_swia, vel_swea = importar()

B_norm = np.linalg.norm(B, axis=1)

# primero diezmamos
paso = int(len(B) / len(vel_swia))
B_diezmado = B[::paso]
t_diezmado = tmag[::paso]

B_cut = B_diezmado[donde(t_diezmado, 3.75) : donde(t_diezmado, 4.15)]
velocidad = vel_swia[donde(tswia, 3.75) : donde(tswia, 4.15)]

B_medio = np.mean(B_cut, axis=0)
B_medio_normalizado = B_medio / np.linalg.norm(B_medio)

# proyecto la velocidad sobre B
dot = np.dot(velocidad, B_medio_normalizado)
N = np.zeros((len(dot), len(B_medio)))

for i in range(len(N)):
    N[i, :] = dot[i] * B_medio_normalizado

v_perp = np.mean(velocidad - N, axis=0)

gf, rg, rg_min, rg_max = giroradio(B_cut, velocidad)
print(f"La girofrecuencia de iones es {gf:1.3g} rad s⁻¹, {gf / 2 / np.pi:1.3g} s⁻¹")
print(
    f"El radio de Larmor de iones SW es {rg:1.3g} km (entre {rg_min:1.3g} y {rg_max:1.3g})"
)

#################
# ahora empieza el MVA con los datos que elegí

t_2 = donde(tmag, 3.931)
t_3 = donde(tmag, 3.96)
B_cut = B[t_2 : t_3 + 1, :]

M_ij = Mij(B_cut)

# ahora quiero los autovectores y autovalores
[lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

# Los ordeno de mayor a menor
idx = lamb.argsort()[::-1]
lamb = lamb[idx]
x = x[:, idx]
# ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
x1 = x[:, 0]
x2 = x[:, 1]
x3 = x[:, 2]

av = np.concatenate([x1, x2, x3])

if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
    x3 = -x3
if any(np.cross(x1, x2) - x3) > 0.01:
    print("Cambio el signo de x1 para que los av formen terna derecha")
    x1 = -x1

# las proyecciones
B1 = np.dot(B_cut, x1)
B2 = np.dot(B_cut, x2)
B3 = np.dot(B_cut, x3)

# Proyecciones del vector medio en los autovectores principales
B_medio_1 = np.dot(B_medio, x1)
B_medio_2 = np.dot(B_medio, x2)

# Plot hodograma con vector medio
plt.figure()
plt.plot(B2, B1, zorder=1)
plt.scatter(B2[0], B1[0], s=50, zorder=2, marker="o", color="r", label="inicio")
plt.scatter(B2[-1], B1[-1], s=50, zorder=2, marker="x", color="r", label="fin")
plt.scatter(
    B_medio_2,
    B_medio_1,
    s=200,
    zorder=3,
    marker=r"$\otimes$",
    color="k",
    label="B medio",
)
plt.xlabel(f"B2 (nT)", fontsize=16)
plt.ylabel(f"B1 (nT)", fontsize=16)
# plt.grid()
plt.tick_params(axis="both", which="major", labelsize=14)
plt.legend(fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

B_medio_vectorial = np.mean(B_cut, axis=0)
B_norm_medio = np.linalg.norm(B_medio_vectorial)
print(B_medio_vectorial / B_norm_medio)

# Calcular la amplitud de las ondas como la desviación estándar
amplitud_ondas = np.std(B_cut)

# Calcular el cociente entre B medio y la amplitud de las ondas
cociente = amplitud_ondas / np.linalg.norm(B_medio)
print(f"Cociente entre B medio y la amplitud de las ondas: {cociente:1.3g}")

VA = 1e-9 * B_norm_medio / np.sqrt(4e-5 * np.pi * 1.67e-27 * 7) * 1e-5
MA = np.linalg.norm(velocidad, axis=1) / VA
