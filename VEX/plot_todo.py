import matplotlib.pyplot as plt
import numpy as np
import sys
from importar_datos import importar_MAG_pds, importar_ELS_clweb

sys.path.append("..")
from funciones import find_nearest, Bpara_Bperp


energias = [35, 50, 75, 100]

year = 2011
month = "04"
day = 30
doy = 120
ti = 2
tf = 3

t, B, pos = importar_MAG_pds(year, doy, ti, tf)
Bnorm = np.linalg.norm(B, axis=1)
Bpara, Bperp, tpara = Bpara_Bperp(B[::32], t[::32], ti, tf)

t_els, ELS = importar_ELS_clweb(year, month, day, ti, tf)
energy = ELS[:, 7]
JE = ELS[:, -1]

fig = plt.figure(
    1, constrained_layout=True
)  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
fig.subplots_adjust(
    top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
)
# xfmt = md.DateFormatter("%H:%M")
ax1 = plt.gca()

plt.title("Spacebar when ready to click:")
ax1 = plt.subplot2grid((4, 1), (0, 0))
ax1.plot(t, Bnorm, linewidth=0.5)
ax1.set_ylabel("|B| (nT)")
ax1.grid()

ax2 = plt.subplot2grid((4, 1), (1, 0), sharex=ax1)
ax2.plot(t, B[:, 0], label="Bx VSO", linewidth=0.5)
ax2.plot(t, B[:, 1], label="By VSO", linewidth=0.5)
ax2.plot(t, B[:, 2], label="Bz VSO", linewidth=0.5)
ax2.set_ylabel("B components (nT)")
ax2.legend(loc="upper left")
ax2.grid()

ax3 = plt.subplot2grid((4, 1), (2, 0), sharex=ax1)
ax3.plot(tpara, Bpara, linewidth=0.5, label="B ||")
ax3.plot(tpara, Bperp, linewidth=0.5, label="B perp")
ax3.set_ylabel("variación de Bpara perp")
ax3.set_xlabel("Tiempo (hdec)")
ax3.set_ylim([-0.1, 5])
ax3.legend(loc="upper left")
ax3.grid()

ax4 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1)
for energia in energias:
    index = np.where(energy == find_nearest(energy, energia))[
        0
    ]  # no cambiarlo a donde()! Me tiene que dar un array, no un escalar.
    plt.semilogy(
        t_els[index],
        JE[index],
        label=f"{energia} eV",
        linewidth=0.5,
    )
ax4.set_ylabel("Diff energy flux \n of the SW e- \n (cm⁻² sr⁻¹ s⁻¹)")
ax4.legend(loc="center right")
ax4.grid()

plt.show()
