import numpy as np
import matplotlib.pyplot as plt
from funciones import t_clweb, find_nearest, giroradio, donde, Mij
from funciones_plot import hodograma

"""
Plotea sólo los datos de MAG de alta resolución
"""

np.set_printoptions(precision=4)


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


def importar_swea():
    path = "../../datos/oscilitones/"

    swea = np.loadtxt(path + "swea_JE.asc")

    t_swea, idx = np.unique(
        swea[:, 3] + swea[:, 4] / 60 + swea[:, 5] / 3600, return_index=True
    )

    energias = [100 + i * 100 for i in range(5)]

    energy = swea[:, 7]
    JE_total = swea[:, -1]
    JE_cut = np.zeros((len(t_swea), len(energias)))

    for i, energia in enumerate(energias):
        index = np.where(energy == find_nearest(energy, energia))[
            0
        ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
        JE = JE_total[index]
        JE_cut[:, i] = JE

    return swea, t_swea, JE_cut


swea, t_swea, JE = importar_swea()
tmag, tswia, tswea, B, vel_norm_swia, vel_norm_swea, vel_swia, vel_swea = importar()

B_norm = np.linalg.norm(B, axis=1)

fig = plt.figure(
    1, constrained_layout=True
)
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
plt.title("MAVEN MAG SWIA SWEA 24 ENE 2015")

ax1 = plt.subplot2grid((4, 1), (0, 0))
plt.plot(tmag, B_norm, linewidth=1, c="k")
ax1.set_ylabel("|B| (nT)")
ax1.set_xlim([3.7, 4.2])

ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
ax2.plot(tmag, B, linewidth=1, label=["Bx", "By", "Bz"])
ax2.set_ylabel("B (nT)")

ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
ax3.plot(tswia, vel_swia[:, 0] + 300, linewidth=1, label=["Vx"])
ax3.plot(tswia, vel_swia[:, 1:3], linewidth=1, label=["Vy", "Vz"])
ax3.set_ylabel("Vel protones\n SW (km/s)")

ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)
ax4.plot(tswea, vel_swea, linewidth=1, label=["Vx", "Vy", "Vz"])
ax4.set_ylabel("Vel electrones\n SW (km/s)")
ax4.set_xlabel("Tiempo (hdec)")
ax4.grid()

for ax in [ax1, ax2, ax3]:
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.grid()
for ax in [ax4, ax2, ax3]:
    ax.legend(loc="lower left")

plt.show()

fig = plt.figure(
    1, constrained_layout=True
)
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)

ax1 = plt.subplot2grid((3, 1), (0, 0))
plt.plot(tmag, B_norm, linewidth=1, c="k")
ax1.set_title("MAVEN MAG SWIA 24 ENE 2015")
ax1.set_ylabel("|B| (nT)")
ax1.set_xlim([3.74, 4.15])
ax1.set_ylim([0, 10])

ax2 = plt.subplot2grid((3, 1), (1, 0), sharex=ax1)
ax2.plot(tmag, B, linewidth=1, label=["Bx", "By", "Bz"])
ax2.set_ylabel("B (nT)")
ax2.set_ylim([-5, 5])

ax3 = plt.subplot2grid((3, 1), (2, 0), sharex=ax1)
ax4 = ax3.twinx()
ax4.plot(tswia, vel_swia[:, 0], linewidth=1)
ax3.plot(tswia, vel_swia[:, 1], c="C1", linewidth=1)
ax3.plot(tswia, vel_swia[:, 2], c="C2", linewidth=1)
ax3.set_ylabel(r"$v_y$, $v_z$ SW (km/s)")
ax4.set_ylabel(r"$v_x$ SW (km/s)")
ax3.set_xlabel("Tiempo (hdec)")
ax3.set_ylim([0, 50])
ax4.set_ylim([-330, -280])
ax3.grid()


for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.grid()
for ax in [ax2]:
    ax.legend(loc="lower left")

plt.show()

# primero diezmamos
paso = int(len(B) / len(vel_swia))
B_diezmado = B[::paso]
t_diezmado = tmag[::paso]

B_cut = B_diezmado[donde(t_diezmado, 3.75):donde(t_diezmado, 4.15)]
velocidad = vel_swia[donde(tswia, 3.75):donde(tswia, 4.15)]

B_medio = np.mean(B_cut, axis=0)
B_medio_normalizado = B_medio / np.linalg.norm(B_medio)

# proyecto la velocidad sobre B
dot = np.dot(velocidad, B_medio_normalizado)
N = np.zeros((len(dot), len(B_medio)))

for i in range(len(N)):
    N[i, :] = dot[i] * B_medio_normalizado

v_perp = np.mean(velocidad - N, axis=0)

gf, rg, rg_min, rg_max = giroradio(B_cut, velocidad)
print(f"La girofrecuencia de iones es {gf:1.3g} rad s⁻¹, {gf/2/np.pi:1.3g} s⁻¹")
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

plt.figure()
plt.plot(B2, B1, zorder=1)
plt.scatter(B2[0], B1[0], s=50, zorder=2, marker="o", color="r", label="start")
plt.scatter(B2[-1], B1[-1], s=50, zorder=2, marker="x", color="r", label="end")
plt.xlabel(f"B2 (nT)", fontsize=16)
plt.ylabel(f"B1 (nT)", fontsize=16)
plt.grid()
plt.tick_params(axis="both", which="major", labelsize=14)
# plt.suptitle("MAVEN MAG MVA 10-10-2015", fontsize=18)
plt.legend(fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.show()
B_medio_vectorial = np.mean(B_cut, axis=0)
B_norm_medio = np.linalg.norm(B_medio_vectorial)
print(B_medio_vectorial/B_norm_medio)