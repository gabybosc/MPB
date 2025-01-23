import numpy as np

from _importar_datos import (
    importar_MAG,
    importar_t1t2t3t4,
    importar_tMVA,
    importar_fila,
    importar_BupBdown,
)
from _update_parametros import update_varios

from _fit import (
    calcular_normal,
    hallar_phi,
    rotacion,
)
from time import sleep
import matplotlib.pyplot as plt
import sys

sys.path.append("..")
from funciones import (
    donde,
    fechas,
    Mij,
    autovectores,
    SZA,
    UTC_to_hdec,
    angulo,
    ancho_mpb,
    corrientes,
    day_to_doy,
    hdec_to_UTC,
    find_nearest_final,
    find_nearest_inicial,
)
from funciones_metodos import bootstrap
from funciones_plot import onpick1

np.set_printoptions(precision=4)

"""
Calcula n mva, n fit, x14, x23, long in, rg, j para poner en una tabla.
"""

"""
calcula las normales de bootstrap y del fit, nada más.
"""

lista = np.genfromtxt("lista-fechas.txt", dtype="str", skip_header=0, skip_footer=12)

for i in range(len(lista)):
    sleep(30)
    year = lista[i, 2]
    month = lista[i, 1]
    day = lista[i, 0]
    doy = day_to_doy(year, month, day)[1]
    print(year, month, day, doy)

    # year, month, day, doy = fechas()
    # t1, t2, t3, t4 = 4.5201, 4.5249, 4.5273, 4.5427

    t1, t2, t3, t4 = importar_t1t2t3t4(year, month, day)

    ti, tf = importar_tMVA(year, month, day)
    Bup, Bdown = importar_BupBdown(year, month, day)

    t, B, posicion, cl, tpos = importar_MAG(year, doy, ti - 1, tf + 1)
    Bnorm = np.linalg.norm(B, axis=1)

    inicio_MVA = donde(t, ti)
    fin_MVA = donde(t, tf)
    B_cut = B[inicio_MVA:fin_MVA]

    n_mva, ang, delta_B3, out, out_phi = bootstrap(1000, B_cut)

    # print("la normal del bootstrap es ", n_mva)

    pos_MPB = int(0.5 * (donde(tpos, ti) + donde(tpos, tf)))

    R = posicion[pos_MPB, :]
    R_2d = np.array([R[0], np.sqrt(R[1] ** 2 + R[2] ** 2)])

    # sza_rad = SZA(posicion, pos_MPB) / 180 * np.pi
    # normal_2d = calcular_normal(sza_rad)
    #
    # """
    # A partir de la normal 2D, la puedo rotar y encontrar la normal 3D
    # Para eso, necesito hallar el ángulo phi primero
    # """
    #
    # phi = hallar_phi(R)[2]
    # n_fit = rotacion(phi, normal_2d)
    #
    # # print("la normal del fit es ", n_fit)
    #
    # angulo_mva = np.arccos(np.clip(np.dot(n_mva, n_fit), -1.0, 1.0))
    #
    # # print(
    # #     f"El ángulo entre las normales 3D de MVA y del fit es {angulo_mva * 180 / np.pi:.3g}º"
    # # )
    #
    # B_mpb = B[donde(t, t1) : donde(t, t4)]
    # B_mva = B[donde(t, ti) : donde(t, tf)]  # el mva es entre t2 y tmva
    #
    # B_medio_mpb = np.mean(B_mpb, axis=0)
    #
    # B3_mva = np.dot(B_mpb, n_mva)
    # B3_fit = np.dot(B_mpb, n_fit)
    #
    # def velocidad(posicion_cut, tpos):
    #     M = len(posicion_cut)
    #     v_punto = np.zeros((M - 1, 3))
    #     deltat = np.zeros(M - 1)
    #     if np.mean(posicion_cut) < 10:  # es decir, pos está en RV en vez de km
    #         posicion_cut = posicion_cut * 6050
    #     for i in range(len(v_punto)):
    #         deltat[i] = (tpos[i + 1] - tpos[i]) * 3600  # delta t en segundos
    #         v_punto[i] = (posicion_cut[i + 1, :] - posicion_cut[i, :]) / deltat[i]
    #         # en km/s
    #     # la velocidad promedio
    #     v_media = np.mean(v_punto, axis=0)
    #     return v_media
    #
    # posi = donde(tpos, t1)
    # posf = donde(tpos, t4)
    #
    # # vel_VEX = (posicion[posf] - posicion[posi]) * 6050 / ((tpos[posf] - tpos[posi]) * 3600)  # para el 30 de abril que tiene un solo punto
    # vel_VEX = velocidad(posicion[posi:posf], tpos[posi:posf])
    #
    # theta_B_mva = angulo(np.mean(B_mpb, axis=0), n_mva) * 180 / np.pi
    # theta_B_fit = angulo(np.mean(B_mpb, axis=0), n_fit) * 180 / np.pi
    #
    # theta_v_mva = angulo(vel_VEX, n_mva) * 180 / np.pi
    # theta_v_fit = angulo(vel_VEX, n_fit) * 180 / np.pi
    #
    # x14_mva, x23_mva = ancho_mpb(t1, t2, t3, t4, n_mva, np.linalg.norm(vel_VEX))
    # x14_fit, x23_fit = ancho_mpb(t1, t2, t3, t4, n_fit, np.linalg.norm(vel_VEX))
    #
    # print(x14_mva, x23_mva)
    #
    # js_mva, jv_mva = corrientes(n_mva, Bup, Bdown, x23_mva)
    # js_fit, jv_fit = corrientes(n_fit, Bup, Bdown, x23_fit)
    #
    # fuerza_mva = np.cross(
    #     jv_mva * 1e-9, B[donde(t, t4 + 0.015), :] * 1e-9
    # )  # N/m^3 #en t4
    # # print(f"la fuerza para el día {year}-{month}-{day} es {fuerza_mva}")
    #
    # """
    # Calcula el giroradio y la long inercial de iones. El giroradio lo calcula usando
    # v_perp al campo B y también usando la velocidad proyectada en la normal.
    # """

    # #########CONSTANTES
    mp = 1.67e-27  # masa del proton en kg
    kB = 1.38e-23  # cte de Boltzmann en J/K
    q_e = 1.602e-19  # carga del electrón en C

    IMA = np.genfromtxt(
        f"dens_{year}{month}{day}-futaana.txt", dtype=str, skip_header=1
    )

    t_IMA = np.array(
        [(UTC_to_hdec(IMA[i, 0]) + UTC_to_hdec(IMA[i, 1])) / 2 for i in range(len(IMA))]
    )

    dens = [float(IMA[i, 2]) for i in range(len(IMA))]
    vel_protones = np.array(
        [
            [float(IMA[i, 3]), float(IMA[i, 4]), float(IMA[i, 5])]
            for i in range(len(IMA))
        ]
    )

    # IMA = np.loadtxt(
    #     "densidad_20140426-clweb.txt", skiprows=1
    # )  # estos son los del clweb
    # t_IMA = np.array(IMA[:, 0] + IMA[:, 1] / 60 + IMA[:, 2] / 3600)  # hdec
    # dens = IMA[:, -1]

    happy = False
    while not happy:
        val = []
        while len(val) < 2:
            plt.clf()  # clear figure
            fig = plt.figure(
                1, constrained_layout=True
            )  # Lo bueno de esta forma es que puedo hacer que solo algunos compartan eje
            fig.subplots_adjust(
                top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.005, wspace=0.15
            )
            # xfmt = md.DateFormatter("%H:%M")

            plt.title("Spacebar when ready to click:")
            ax1 = plt.subplot2grid((2, 1), (0, 0))
            ax1.plot(t, Bnorm, linewidth=1)
            ax1.set_ylabel("|B| (nT)")
            ax1.set_title(f"VEX MAG {year}-{month}-{day}")
            ax1.grid()

            ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1)
            ax2.scatter(t_IMA, dens)
            ax2.set_ylabel("B components (nT)")
            ax2.grid()
            for xc in [t1, t2, t3, t4]:
                for ax in ax1, ax2:
                    ax.axvline(x=xc, color="k", linewidth=1)

            fig.canvas.mpl_connect("pick_event", onpick1)

            zoom_ok = False
            print("\nSpacebar when ready to click:\n")
            while not zoom_ok:
                zoom_ok = plt.waitforbuttonpress(-1)
            print("Click to select magnetosheath: ")
            val = np.asarray(plt.ginput(2))[:, 0]
            print("Selected values: ", val)
            outs = sorted(val)

        print("Happy? Keyboard click for yes, mouse click for no.")
        happy = plt.waitforbuttonpress()

    plt.show()

    """
    giroradio: rg = mp * v_perp / (q_e * B)  en la región upstream
    where vperp is the component of the velocity perpendicular to the
    direction of the magnetic field and B is the strength of the magnetic field
    """

    # selecciono la región upstream usando un único punto y calculo en ese lugar.
    idx1 = np.where(t_IMA == find_nearest_inicial(t_IMA, outs[0]))[0][0]
    idx2 = np.where(t_IMA == find_nearest_final(t_IMA, outs[1]))[0][0]
    density_mean = np.mean(dens[idx1:idx2])
    vel_upstream = np.mean(vel_protones[idx1:idx2], axis=0)
    print(
        f"La densidad de protones en la MS es {density_mean} y la velocidad media es {vel_upstream}"
    )
    B_punto = (
            np.mean(B[donde(t, t_IMA[idx1]): donde(t, t_IMA[idx2])], axis=0) * 1e-9
    )  # para que me de en km

    dot_product = np.dot(vel_upstream, B_punto)

    # Calcular el módulo del campo magnético
    B_norm_squared = np.dot(B_punto, B_punto)

    # Calcular la componente paralela de la velocidad
    v_parallel = (dot_product / B_norm_squared) * B_punto

    # Calcular la componente perpendicular de la velocidad
    v_perp = vel_upstream - v_parallel

    # el giroradio entonces:
    gf = (q_e * np.linalg.norm(B_punto)) / mp
    rg = mp * np.linalg.norm(v_perp) / (q_e * np.linalg.norm(B_punto))

    print(f"El radio de Larmor es {rg:1.3g} km ")

    """
    Longitud inercial
    """
    ion_length = 2.28e07 / np.sqrt(density_mean) * 1e-5  # km
    print(f"La longitud inercial de iones es {ion_length:1.3g} km")

    """
    actualizamos gdocs
    """
    nr, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_Ajuste = importar_fila(
        year, month, day
    )

    # hoja = hoja_parametros
    # hoja.update_acell(f"D{nr}", f"{SZA(posicion, pos_MPB):.3g}")
    # hoja.update_acell(f"S{nr}", f"{np.linalg.norm(vel_VEX):.3g}")
    #
    # cell_vel = hoja.range(f"P{nr}:R{nr}")
    #
    # update_varios(hoja, cell_vel, vel_VEX)
    #
    hoja = hoja_MVA

    # hoja.update_acell(f"T{nr}", f"{round(np.mean(B3_mva), 2)}")
    # hoja.update_acell(
    #     f"V{nr}",
    #     f"{abs(round(np.mean(B3_mva) / np.linalg.norm(np.mean(B_mpb, axis=0))))}",
    # )
    # hoja.update_acell(f"W{nr}", f"{theta_v_mva:.3g}")
    # hoja.update_acell(f"X{nr}", f"{theta_B_mva:.3g}")
    # hoja.update_acell(f"Y{nr}", f"{np.linalg.norm(x23_mva):.3g}")
    # hoja.update_acell(f"Z{nr}", f"{np.linalg.norm(x14_mva):.3g}")
    hoja.update_acell(f"AA{nr}", f"{ion_length:.3g}")  # long inercial
    hoja.update_acell(f"AC{nr}", f"{rg:.3g}")  # giroradio

    # hoja.update_acell(f"AG{nr}", f"{np.linalg.norm(js_mva):.3g}")
    # hoja.update_acell(f"AK{nr}", f"{np.linalg.norm(jv_mva):.3g}")
    # hoja.update_acell(f"AL{nr}", f"{np.linalg.norm(fuerza_mva):.3g}")
    #
    # cell_Js = hoja.range(f"AD{nr}:AF{nr}")
    # cell_Jv = hoja.range(f"AH{nr}:AJ{nr}")
    # update_varios(hoja, cell_Js, js_mva)
    # update_varios(hoja, cell_Jv, jv_mva)
    #
    # hoja = hoja_Ajuste
    # hoja.update_acell(f"J{nr}", f"{angulo_mva * 180 / np.pi:.3g}")
    # hoja.update_acell(f"M{nr}", f"{theta_v_fit:.3g}")
    # hoja.update_acell(f"N{nr}", f"{theta_B_fit:.3g}")
    # hoja.update_acell(f"O{nr}", f"{x23_fit:.3g}")
    # hoja.update_acell(f"P{nr}", f"{x14_fit:.3g}")
    #
    # hoja.update_acell(f"V{nr}", f"{np.linalg.norm(js_fit):.3g}")
    # hoja.update_acell(f"Z{nr}", f"{np.linalg.norm(jv_fit):.3g}")
    # # hoja.update_acell(f'AA{nr}', f'{:.3g}')  # long inercial
    # # hoja.update_acell(f'AC{nr}', f'{:.3g}')  # giroradio
    #
    # cell_Js = hoja.range(f"S{nr}:U{nr}")
    # cell_Jv = hoja.range(f"W{nr}:Y{nr}")
    # # cell_norm = hoja.range(f"G{nr}:I{nr}")
    # update_varios(hoja, cell_Js, js_fit)
    # update_varios(hoja, cell_Jv, jv_fit)
    # update_varios(hoja, cell_norm, n_fit)
