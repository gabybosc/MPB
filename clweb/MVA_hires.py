import numpy as np
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import sys
from importar_datos import importar_mag, importar_fila

sys.path.append("..")
from funciones import (
    angulo,
    error,
    find_nearest,
    Mij,
    autovectores,
    donde,
)
from funciones_metodos import bootstrap, plot_bootstrap, ajuste_conico
from funciones_plot import hodograma

"""

Para datos de mag del clweb.

Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o doy-año)
y los tiempos entre los que voy a realizar el MVA.
Eventualmente podría simplemente encontrar todos los cruces que quiero
y decirle que lea directamente de algún lugar eso. (es lo que hace MVA_automatico)
Antes de correr este programa hay que haber usado plot_seleccionar_puntos, para tener los cuatro tiempos elegidos.
Toma los datos de alta resolución y les aplica un filtro pasabajos
con ventana Butterworth con frecuencia de corte de 0.1 Hz de orden 3.
A los datos filtrados les aplica el MVA.
Nos devuelve los lambdas, el cociente de lambdas, el omega, los autovectores y la normal.
Nos da el valor medio de B, de la altitud y el SZA.
Devuelve el ancho de la MPB y la corriente que pasa. Calcula también la fuerza de lorentz y el campo de hall.
Grafica el hodograma, el ajuste de vignes, y la comparación de las normales obtenidas por los distintos métodos.
"""
np.set_printoptions(precision=4)


def acceso_spreadsheet():
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive.file",
        "https://www.googleapis.com/auth/drive",
    ]

    creds = ServiceAccountCredentials.from_json_keyfile_name("../mpb_api.json", scope)

    client = gspread.authorize(creds)

    hoja_parametros = client.open("MPB").worksheet("Parametros")
    hoja_mva = client.open("MPB").worksheet("MVA")
    hoja_boot = client.open("MPB").worksheet("Bootstrap")
    hoja_fit = client.open("MPB").worksheet("Ajuste")

    fecha_sheet = hoja_mva.col_values(1)
    hora_sheet = hoja_mva.col_values(2)

    print("Acceso a la spreadsheet")
    return hoja_parametros, hoja_mva, hoja_boot, hoja_fit, fecha_sheet, hora_sheet


def poner_fecha(hoja, nr, fecha, hora):
    hoja.update_acell(f"A{nr}", f"{fecha}")
    hoja.update_acell(f"B{nr}", f"{hora}")


def MVA(year, month, day, ti_MVA, tf_MVA):
    mag, t, B, posicion = importar_mag(year, month, day, ti_MVA, tf_MVA)

    M = len(t)
    Bnorm = np.linalg.norm(B, axis=1)

    n_p = int(M / 2)

    M_ij = Mij(B)

    # ahora quiero los autovectores y autovalores
    avec, lamb = autovectores(M_ij)

    # las proyecciones
    B1 = np.dot(B, avec[0])
    B2 = np.dot(B, avec[1])
    B3 = np.dot(B, avec[2])

    # el B medio
    B_medio_vectorial = np.mean(B, axis=0)
    altitud = np.linalg.norm(posicion, axis=1) - 3390  # km
    altitud_media = np.mean(altitud)

    SZA = angulo(posicion[n_p, :], [1, 0, 0]) * 180 / np.pi

    if posicion[n_p, 2] < 0:
        SZA = -SZA

    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    hodograma(B1, B2, B3)

    # el error
    phi, delta_B3 = error(lamb, B, avec[2])
    if phi[2, 1] > phi[2, 0]:
        error_normal = phi[2, 1] * 180 / np.pi
    else:
        error_normal = phi[2, 0] * 180 / np.pi

    print("Fin del MVA.")

    # ######update la hoja de los parámetros

    nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(
        year, month, day, int(ti_MVA)
    )

    hoja_parametros.update_acell(f"D{nr}", f"{SZA:.3g}")
    hoja_parametros.update_acell(f"E{nr}", f"{int(altitud_media)}")
    hoja_parametros.update_acell(f"O{nr}", f"{B_norm_medio}")

    cell_B = hoja_parametros.range(f"L{nr}:N{nr}")
    for i, cell in enumerate(cell_B):
        cell.value = B_medio_vectorial[i]
    hoja_parametros.update_cells(cell_B)

    # #######update la hoja de MVA
    hoja_mva.update_acell(f"D{nr}", f"{ti_MVA}")
    hoja_mva.update_acell(f"E{nr}", f"{tf_MVA}")
    hoja_mva.update_acell(f"I{nr}", f"{lamb[1] / lamb[2]:.3g}")

    hoja_mva.update_acell(f"S{nr}", f"{error_normal:.3g}")
    hoja_mva.update_acell(f"T{nr}", f"{np.mean(B3)}")
    hoja_mva.update_acell(f"U{nr}", f"{delta_B3}")
    hoja_mva.update_acell(f"V{nr}", f"{abs(np.mean(B3) / B_norm_medio)}")

    cell_lambda = hoja_mva.range(f"F{nr}:H{nr}")
    for i, cell in enumerate(cell_lambda):
        cell.value = lamb[i]
    hoja_mva.update_cells(cell_lambda)

    # cell_av = hoja_mva.range(f"J{nr}:R{nr}")
    # for i, cell in enumerate(cell_av):
    #     cell.value = avec[i]
    hoja_mva.update_acell(f"P{nr}", f"{avec[2]}")
    print("Escribió la spreadsheet del MVA.")
    return avec[2], B, t, posicion, nr


def ajuste(year, month, day, doy, ti_MVA, tf_MVA, nr):
    datos_tiempo = np.loadtxt("../outputs/t1t2t3t4.txt")
    idx_d = np.where(int(doy) == datos_tiempo[:, 1].astype(int))[0]
    idx_h = np.where(int(ti_MVA) == datos_tiempo[:, 2].astype(int))[0]
    idx = np.intersect1d(idx_d, idx_h)[0]
    t1 = datos_tiempo[idx, 2]
    t2 = datos_tiempo[idx, 3]
    t3 = datos_tiempo[idx, 4]
    t4 = datos_tiempo[idx, 5]
    tiempos = [t1, t2, t3, t4]

    mag, t, B, posicion = importar_mag(year, month, day, ti_MVA, tf_MVA)

    t_nave = find_nearest(t, (t2 + t3) / 2)
    # el tiempo en el medio de la hoja de corriente
    index = donde(t, t_nave)
    x0 = 0.78
    e = 0.9
    X1, Y1, Z1, L0, normal_ajuste = ajuste_conico(posicion[index])

    B3_fit = np.dot(B, normal_ajuste)
    print("Fin del ajuste. ")

    nr, hoja_parametros, hoja_mva, hoja_boot, hoja_fit = importar_fila(
        year, month, day, int(ti_MVA)
    )

    cell_times = hoja_parametros.range(f"F{nr}:I{nr}")
    for i, cell in enumerate(cell_times):
        cell.value = tiempos[i]
        hoja_parametros.update_cells(cell_times)

    hoja_fit.update_acell(f"D{nr}", f"{L0}")
    hoja_fit.update_acell(f"E{nr}", f"{e}")
    hoja_fit.update_acell(f"F{nr}", f"{x0}")

    cell_normal = hoja_fit.range(f"G{nr}:I{nr}")
    for i, cell in enumerate(cell_normal):
        cell.value = normal_ajuste[i]
    hoja_fit.update_cells(cell_normal)

    hoja_fit.update_acell(f"K{nr}", f"{np.mean(B3_fit)}")
    print("Escribió la spreadsheet del ajuste.")

    return normal_ajuste, t1, t2, t3, t4


def bootstrap_completo(B, M, nr, N=1000):
    normal_boot, phi, delta_B3, out, out_phi = bootstrap(N, B)

    muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

    B3_boot = np.dot(B, normal_boot)

    B_medio_vectorial = np.mean(B, axis=0)
    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    if sigma31 > sigma32:
        error_boot = sigma31
    else:
        error_boot = sigma32

    print("Fin del bootstrap. ")

    (
        hoja_parametros,
        hoja_mva,
        hoja_boot,
        hoja_fit,
        fecha_sheet,
        hora_sheet,
    ) = acceso_spreadsheet()

    hoja_boot.update_acell(f"D{nr}", f"{M}")
    hoja_boot.update_acell(f"E{nr}", f"{N}")

    cell_normal = hoja_boot.range(f"F{nr}:H{nr}")
    for i, cell in enumerate(cell_normal):
        cell.value = normal_boot[i]
    hoja_boot.update_cells(cell_normal)

    hoja_boot.update_acell(f"I{nr}", f"{error_boot:.3g}")
    hoja_boot.update_acell(f"J{nr}", f"{np.mean(B3_boot)}")
    hoja_boot.update_acell(f"K{nr}", f"{sigmaB}")
    hoja_boot.update_acell(f"L{nr}", f"{abs(np.mean(B3_boot) / B_norm_medio)}")
    print("Escribió la spreadsheet del bootstrap.")

    return normal_boot
