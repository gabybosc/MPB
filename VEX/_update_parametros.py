import numpy as np

# import gspread
# from oauth2client.service_account import ServiceAccountCredentials
# from importar_datos import importar_MAG_pds, importar_ELS_clweb, importar_fila

np.set_printoptions(precision=4)


def hoja_MVA_update(hoja, nr, lamb, av, error_normal, B3, delta_B3, B_norm_medio):
    hoja.update_acell(f"I{nr}", f"{lamb[1] / lamb[2]:.3g}")

    hoja.update_acell(f"S{nr}", f"{error_normal:.3g}")
    hoja.update_acell(f"T{nr}", f"{round(np.mean(B3), 2)}")
    hoja.update_acell(f"U{nr}", f"{round(delta_B3, 2)}")
    hoja.update_acell(f"V{nr}", f"{abs(round(np.mean(B3) / B_norm_medio, 2))}")

    cell_lambda = hoja.range(f"F{nr}:H{nr}")
    for i, cell in enumerate(cell_lambda):
        cell.value = round(lamb[i], 2)
    hoja.update_cells(cell_lambda)

    cell_av = hoja.range(f"J{nr}:R{nr}")
    for i, cell in enumerate(cell_av):
        cell.value = round(av[i], 3)
    hoja.update_cells(cell_av)


def hoja_MVA_analisis(hoja, nr, ti, tf, x14, x23, J_s, J_v, fuerza):
    hoja.update_acell(f"D{nr}", f"{ti}")
    hoja.update_acell(f"E{nr}", f"{tf}")
    # hoja.update_acell(f"W{nr}", f"{angulo_v_mva * 180/np.pi:.3g}")
    # hoja.update_acell(f"X{nr}", f"{angulo_B_mva * 180/np.pi:.3g}")
    hoja.update_acell(f"Y{nr}", f"{np.linalg.norm(x23):.3g}")
    hoja.update_acell(f"Z{nr}", f"{np.linalg.norm(x14):.3g}")
    # hoja.update_acell(f'AA{nr}', f'{:.3g}')
    # hoja.update_acell(f'AB{nr}', f'{:.3g}')

    hoja.update_acell(f"AF{nr}", f"{np.linalg.norm(J_s) * 1E-6:.3g}")
    hoja.update_acell(f"AJ{nr}", f"{np.linalg.norm(J_v):.3g}")
    hoja.update_acell(f"AK{nr}", f"{np.linalg.norm(fuerza):.3g}")

    cell_Js = hoja.range(f"AC{nr}:AE{nr}")
    for i, cell in enumerate(cell_Js):
        cell.value = round(J_s[i] * 1e-6, 3)
    hoja.update_cells(cell_Js)

    cell_Jv = hoja.range(f"AG{nr}:AI{nr}")
    for i, cell in enumerate(cell_Jv):
        cell.value = round(J_v[i], 3)
    hoja.update_cells(cell_Jv)


def hoja_param(hoja, nr, sza, v_media, omega, B_upstream, B_downstream):
    hoja.update_acell(f"D{nr}", f"{sza:.3g}")
    hoja.update_acell(f"S{nr}", f"{np.linalg.norm(v_media):.3g}")
    hoja.update_acell(f"Z{nr}", f"{omega * 180 / np.pi:.3g}")
    cell_vel = hoja.range(f"P{nr}:R{nr}")
    for i, cell in enumerate(cell_vel):
        cell.value = round(v_media[i], 2)
    hoja.update_cells(cell_vel)

    cell_Bup = hoja.range(f"T{nr}:V{nr}")
    for i, cell in enumerate(cell_Bup):
        cell.value = round(B_upstream[i], 2)
    hoja.update_cells(cell_Bup)

    cell_Bdown = hoja.range(f"W{nr}:Y{nr}")
    for i, cell in enumerate(cell_Bdown):
        cell.value = round(B_downstream[i], 2)
    hoja.update_cells(cell_Bdown)


def hoja_t1t2t3t4(hoja, nr, t1, t2, t3, t4):
    hoja.update_acell(f"F{nr}", f"{t1}")
    hoja.update_acell(f"G{nr}", f"{t2}")
    hoja.update_acell(f"H{nr}", f"{t3}")
    hoja.update_acell(f"I{nr}", f"{t4}")


def hoja_bootstrap_p2(
        hoja_boot,
        nr,
        angulo_v,
        angulo_B,
        x_23,
        x_14,
        J_s,
        J_v,
        fuerza,
        E_Hall,
):
    hoja_boot.update_acell(f"M{nr}", f"{angulo_v * 180 / np.pi:.3g}")
    hoja_boot.update_acell(f"N{nr}", f"{angulo_B * 180 / np.pi:.3g}")
    hoja_boot.update_acell(f"O{nr}", f"{np.linalg.norm(x_23):.3g}")
    hoja_boot.update_acell(f"P{nr}", f"{np.linalg.norm(x_14):.3g}")

    hoja_boot.update_acell(f"V{nr}", f"{np.linalg.norm(J_s) * 1E-6:.3g}")
    hoja_boot.update_acell(f"Z{nr}", f"{np.linalg.norm(J_v):.3g}")
    hoja_boot.update_acell(f"AA{nr}", f"{np.linalg.norm(fuerza):.3g}")
    hoja_boot.update_acell(f"AE{nr}", f"{np.linalg.norm(E_Hall) * 1E3:.3g}")

    cell_Js = hoja_boot.range(f"S{nr}:U{nr}")
    for i, cell in enumerate(cell_Js):
        cell.value = round(J_s[i], 3)
    hoja_boot.update_cells(cell_Js)

    cell_Jv = hoja_boot.range(f"W{nr}:Y{nr}")
    for i, cell in enumerate(cell_Jv):
        cell.value = round(J_v[i], 3)
    hoja_boot.update_cells(cell_Jv)

    cell_EH = hoja_boot.range(f"AB{nr}:AD{nr}")
    for i, cell in enumerate(cell_EH):
        cell.value = round(E_Hall[i] * 1e3, 3)
    hoja_boot.update_cells(cell_EH)


def hoja_bootstrap_p1(hoja_boot, nr, M, N, normal, error, B3, B_norm_medio, sigmaB):
    hoja_boot.update_acell(f"D{nr}", f"{M}")
    hoja_boot.update_acell(f"E{nr}", f"{N}")

    cell_normal = hoja_boot.range(f"F{nr}:H{nr}")
    for i, cell in enumerate(cell_normal):
        cell.value = round(normal[i], 3)
    hoja_boot.update_cells(cell_normal)

    hoja_boot.update_acell(f"I{nr}", f"{error:.3g}")
    hoja_boot.update_acell(f"J{nr}", f"{round(np.mean(B3), 2)}")
    hoja_boot.update_acell(f"K{nr}", f"{round(sigmaB, 2)}")
    hoja_boot.update_acell(f"L{nr}", f"{abs(round(np.mean(B3) / B_norm_medio, 2))}")
