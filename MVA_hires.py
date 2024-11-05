import numpy as np
from funciones import error, donde, SZA, Mij
from funciones_metodos import ajuste_conico, plot_bootstrap, bootstrap
from funciones_plot import hodograma
from importar_datos import importar_fila

"""

Para datos de mag de alta resolución. Para que tarde menos en cargar, skipea las primeras t1*32*3600 rows,
lo malo de eso es que a veces no alcanza con esos datos.

Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o doy-año)
y los tiempos entre los que voy a realizar el MVA. Eventualmente podría simplemente
encontrar todos los cruces que quiero y decirle que lea directamente de algún lugar eso (es lo que hace MVA_automatico)
Antes de correr este programa hay que haber usado plot_seleccionar_puntos, para tener los cuatro tiempos elegidos.
Toma los datos de alta resolución y les aplica un filtro pasabajos con ventana Butterworth con frecuencia de corte de 0.1 Hz de orden 3.
A los datos filtrados les aplica el MVA.
Nos devuelve los lambdas, el cociente de lambdas, el omega, los autovectores y la normal.
Nos da el valor medio de B, de la altitud y el SZA.
Devuelve el ancho de la MPB y la corriente que pasa. Calcula también la fuerza de lorentz y el campo de hall.
Grafica el hodograma, el ajuste de vignes, y la comparación de las normales obtenidas por los distintos métodos.
"""
np.set_printoptions(precision=4)


def MVA(year, month, day, doy, ti_MVA, tf_MVA, t, B, posicion):
    fila, hoja_parametros, hoja_MVA, hoja_Bootstrap, hoja_fit = importar_fila(
        year, month, day, doy, ti_MVA
    )
    t1 = float(hoja_parametros.cell(fila, 6).value)
    t2 = float(hoja_parametros.cell(fila, 7).value)
    t3 = float(hoja_parametros.cell(fila, 8).value)
    t4 = float(hoja_parametros.cell(fila, 9).value)

    inicio = donde(t, ti_MVA)
    fin = donde(t, tf_MVA)

    B_cut = B[inicio:fin, :]

    # ahora empieza el MVA con los datos que elegí
    posicion_cut = posicion[inicio : fin + 1, :]
    altitud = np.mean(np.linalg.norm(posicion_cut, axis=1) - 3390)

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

    # el B medio
    B_medio_vectorial = np.mean(B_cut, axis=0)
    SZAngle = SZA(posicion_cut, 0)
    if any(posicion_cut[:, 2]) < 0:
        SZAngle = -SZAngle

    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    hodograma(B1, B2, B3)
    # el error
    phi, delta_B3 = error(lamb, B_cut, x)

    ###############
    # ###fit
    orbita = posicion[donde(t, t1 - 1) : donde(t, t4 + 1)] / 3390  # radios marcianos

    index = donde(t, (t2 + t3) / 2)
    x0 = 0.78
    e = 0.9
    L0 = 0.96
    normal_fit, X1, Y1, Z1, R, L0 = ajuste_conico(
        posicion, index, orbita, x3, x0, e, L0
    )

    B3_fit = np.dot(B_cut, normal_fit)

    ###############
    # #Bootstrap

    N_boot = 1000
    normal_boot, phi_boot, delta_B3_boot, out, out_phi = bootstrap(N_boot, B_cut)

    muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

    B3_boot = np.dot(B_cut, normal_boot)

    #######
    # Errores
    if phi[2, 1] > phi[2, 0]:
        error_normal = phi[2, 1] * 57.2958
    else:
        error_normal = phi[2, 0] * 57.2958
        # quiero ver si el error más grande es phi31 o phi32

    if sigma31 > sigma32:
        error_boot = sigma31
    else:
        error_boot = sigma32

    print("Fin del MVA")
    #############
    # ahora guardo todo en una spreadsheet
    # ######updateo la hoja de los parámetros
    # hoja_parametros.update_acell(f"D{fila}", f"{SZAngle:.3g}")
    # hoja_parametros.update_acell(f"E{fila}", f"{int(altitud)}")
    # hoja_parametros.update_acell(f"O{fila}", f"{round(B_norm_medio,2)}")
    #
    # cell_B = hoja_parametros.range(f"L{fila}:N{fila}")
    # for i, cell in enumerate(cell_B):
    #     cell.value = round(B_medio_vectorial[i], 2)
    # hoja_parametros.update_cells(cell_B)
    #
    # # if (
    # #     type(B_filtrado) != int
    # # ):  # si no es un int, en particular 0, agrega los datos a la lista
    # #     hoja_parametros.update_acell(f"J{fila}", f"{orden_filtro}")
    # #     hoja_parametros.update_acell(f"K{fila}", f"{frec_filtro}")
    # # else:
    # #     hoja_parametros.update_acell(f"J{fila}", "Sin filtrar")
    #
    # # #######update la hoja de MVA
    # hoja_MVA.update_acell(f"D{fila}", f"{ti_MVA}")
    # hoja_MVA.update_acell(f"E{fila}", f"{tf_MVA}")
    # hoja_MVA.update_acell(f"I{fila}", f"{lamb[1]/lamb[2]:.3g}")
    #
    # hoja_MVA.update_acell(f"S{fila}", f"{error_normal:.3g}")
    # hoja_MVA.update_acell(f"T{fila}", f"{round(np.mean(B3),2)}")
    # hoja_MVA.update_acell(f"U{fila}", f"{round(delta_B3,2)}")
    # hoja_MVA.update_acell(f"V{fila}", f"{abs(round(np.mean(B3)/B_norm_medio,2))}")
    #
    # cell_lambda = hoja_MVA.range(f"F{fila}:H{fila}")
    # for i, cell in enumerate(cell_lambda):
    #     cell.value = round(lamb[i], 2)
    # hoja_MVA.update_cells(cell_lambda)
    #
    # cell_av = hoja_MVA.range(f"J{fila}:R{fila}")
    # for i, cell in enumerate(cell_av):
    #     cell.value = round(av[i], 3)
    # hoja_MVA.update_cells(cell_av)
    #
    # # #######update la hoja de bootstrap
    # hoja_Bootstrap.update_acell(f"D{fila}", f"{len(B)}")
    # hoja_Bootstrap.update_acell(f"E{fila}", f"{N_boot}")
    #
    # cell_normal = hoja_Bootstrap.range(f"F{fila}:H{fila}")
    # for i, cell in enumerate(cell_normal):
    #     cell.value = round(normal_boot[i], 3)
    # hoja_Bootstrap.update_cells(cell_normal)
    #
    # hoja_Bootstrap.update_acell(f"I{fila}", f"{error_boot:.3g}")
    # hoja_Bootstrap.update_acell(f"J{fila}", f"{round(np.mean(B3_boot),2)}")
    # hoja_Bootstrap.update_acell(f"K{fila}", f"{round(sigmaB,2)}")
    # hoja_Bootstrap.update_acell(
    #     f"L{fila}", f"{abs(round(np.mean(B3_boot)/B_norm_medio,2))}"
    # )
    #
    # # #######update la hoja del ajuste
    # hoja_fit.update_acell(f"D{fila}", f"{L0}")
    # hoja_fit.update_acell(f"E{fila}", f"{e}")
    # hoja_fit.update_acell(f"F{fila}", f"{x0}")
    #
    # cell_normal = hoja_fit.range(f"G{fila}:I{fila}")
    # for i, cell in enumerate(cell_normal):
    #     cell.value = round(normal_fit[i], 3)
    # hoja_fit.update_cells(cell_normal)
    #
    # hoja_fit.update_acell(f"K{fila}", f"{round(np.mean(B3_fit),2)}")
    # hoja_fit.update_acell(f"L{fila}", f"{abs(round(np.mean(B3_fit)/B_norm_medio,2))}")
    #
    # print("escribe la spreadsheet")
    return (
        x3,
        normal_boot,
        normal_fit,
        inicio,
        fin,
        B_cut,
        t1,
        t2,
        t3,
        t4,
        B_medio_vectorial,
        fila,
        hoja_parametros,
        hoja_MVA,
        hoja_Bootstrap,
        hoja_fit,
    )


def ajuste():
    return None
