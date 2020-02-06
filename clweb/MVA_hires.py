import numpy as np
from funciones import error, find_nearest, find_nearest_final, find_nearest_inicial, Mij, next_available_row
from funciones_MVA import ajuste_conico, plot_bootstrap, bootstrap
from funciones_plot import hodograma
import gspread
import os
from oauth2client.service_account import ServiceAccountCredentials


"""

Para datos de mag del clweb.

Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o doy-año) y los tiempos entre los que voy a realizar el MVA.
Eventualmente podría simplemente encontrar todos los cruces que quiero y decirle que lea directamente de algún lugar eso. (es lo que hace MVA_automatico)
Antes de correr este programa hay que haber usado plot_seleccionar_puntos, para tener los cuatro tiempos elegidos.
Toma los datos de alta resolución y les aplica un filtro pasabajos con ventana Butterworth con frecuencia de corte de 0.1 Hz de orden 3.
A los datos filtrados les aplica el MVA.
Nos devuelve los lambdas, el cociente de lambdas, el omega, los autovectores y la normal.
Nos da el valor medio de B, de la altitud y el SZA.
Devuelve el ancho de la MPB y la corriente que pasa. Calcula también la fuerza de lorentz y el campo de hall.
Grafica el hodograma, el ajuste de vignes, y la comparación de las normales obtenidas por los distintos métodos.
"""

"""
USAR PDB PARA DEBUGGEAR

"""
np.set_printoptions(precision=4)


def MVA(year, month, day, doy, ti_MVA, tf_MVA):
    date_entry = f'{year}-{month}-{day}'
    datos_tiempo = np.loadtxt('../outputs/t1t2t3t4.txt')
    idx_d = np.where(int(doy) == datos_tiempo[:,1].astype(int))[0]
    idx_h = np.where(int(ti_MVA) == datos_tiempo[:,2].astype(int))[0]
    idx = np.intersect1d(idx_d, idx_h)[0]
    t1 = datos_tiempo[idx,2]
    t2 = datos_tiempo[idx,3]
    t3 = datos_tiempo[idx,4]
    t4 = datos_tiempo[idx,5]
    tiempos = np.array([t1,t2,t3,t4])

    path = f'../../../datos/clweb/{year}-{month}-{day}/'  # path a los datos desde la laptop
    if os.path.isfile(path + 'mag_filtrado.txt'):
        mag = np.loadtxt(path + 'mag_filtrado.txt', skiprows=2)
        B = mag[:, :3]
        Bnorm = mag[:,-1]
        mag = np.loadtxt(path + 'MAG.asc')

    else:
        mag = np.loadtxt(path + 'MAG.asc')
        B = mag[:, 6:9]
        Bnorm = np.linalg.norm(B, axis=1)

    hh = mag[:,3]
    mm = mag[:,4]
    ss = mag[:,5]

    t = hh + mm/60 + ss/3600  # hdec

    M = np.size(t)  # el numero de datos

    # tengo que asegurarme de que no haya agujeros en mis datos

    # el campo
    # B = np.zeros((M, 3))
    # for i in range(6,9):
    #     B[:,i-6] = mag[:, i]
    #
    # Bnorm = mag[:,-1]

    # la posición(x,y,z)
    posicion = np.zeros((M, 3))
    for i in range(9,12):
        posicion[:,i-9] = mag[:, i]

    # la matriz diaria:
    MD = np.zeros((M, 9))
    MD[:, 0] = t
    for i in range(1,4):
        MD[:, i] = B[:,i-1]
    MD[:,4] = Bnorm
    for i in range(5,8):
        MD[:,i] = posicion[:,i-5]/3390  # en radios marcianos
    MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390  # altitud en km

    inicio = np.where(t == find_nearest_inicial(t, ti_MVA))[0][0]
    fin = np.where(t == find_nearest_final(t, tf_MVA))[0][0]

    B_cut = B[inicio:fin,:]

    # ahora empieza el MVA con los datos que elegí
    MD_cut = MD[inicio:fin+1, :]
    M_cut = np.size(MD_cut[:,0])
    posicion_cut = posicion[inicio:fin+1,:]
    n_p = int(len(posicion_cut)/2)

    M_ij = Mij(B_cut)

    # ahora quiero los autovectores y autovalores
    [lamb, x] = np.linalg.eigh(M_ij)  # uso eigh porque es simetrica

    # Los ordeno de mayor a menor
    idx = lamb.argsort()[::-1]
    lamb = lamb[idx]
    x = x[:,idx]

    # ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
    x1 = x[:,0]
    x2 = x[:,1]
    x3 = x[:,2]

    av = np.concatenate([x1,x2,x3])

    if x3[0] < 0:  # si la normal aputna para adentro me la da vuelta
        x3 = - x3
    if any(np.cross(x1,x2) - x3) > 0.01:
        print('Cambio el signo de x1 para que los av formen terna derecha')
        x1 = -x1

    # las proyecciones
    B1 = np.dot(B_cut, x1)
    B2 = np.dot(B_cut, x2)
    B3 = np.dot(B_cut, x3)

    # el B medio
    B_medio_vectorial = np.mean(B_cut, axis=0)
    altitud = np.mean(MD_cut[:,8])
    SZA = np.arccos(np.clip(np.dot(posicion_cut[n_p,:]/np.linalg.norm(posicion_cut[n_p,:]), [1,0,0]), -1.0, 1.0)) * 180/np.pi

    if posicion_cut[n_p, 2] < 0:
        SZA = - SZA

    B_norm_medio = np.linalg.norm(B_medio_vectorial)

    hodograma(B1, B2, B3, date_entry)

    # el error
    phi, delta_B3 = error(lamb, B_cut, M_cut, x)

    ###############
    # ###fit
    orbita = posicion[np.where(t == find_nearest_inicial(t, t1-1))[0][0]:np.where(t == find_nearest_final(t, t4+1))[0][0], :] / 3390  # radios marcianos

    t_nave = find_nearest(t,(t2+t3)/2)  # el tiempo en el medio de la hoja de corriente
    index = np.where(t == t_nave)[0][0]
    x0 = 0.78
    e = 0.9
    normal_fit, X1, Y1, Z1, R, L0 = ajuste_conico(posicion, index, orbita, date_entry, x3)

    B3_fit = np.dot(B_cut, normal_fit)

    ###############
    # Bootstrap

    N_boot = 1000
    normal_boot, phi_boot, delta_B3_boot, out, out_phi = bootstrap(N_boot, B_cut, M_cut)

    muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

    B3_boot = np.dot(B_cut, normal_boot)

    #######
    # Errores
    if phi[2,1] > phi[2,0]:
        error_normal = phi[2,1]*57.2958
    else:
        error_normal = phi[2,0]*57.2958
        # quiero ver si el error más grande es phi31 o phi32

    if sigma31 > sigma32:
        error_boot = sigma31
    else:
        error_boot = sigma32

    print('Fin del MVA')
    #############
    # ahora guardo todo en una spreadsheet
    scope = ["https://spreadsheets.google.com/feeds","https://www.googleapis.com/auth/spreadsheets","https://www.googleapis.com/auth/drive.file","https://www.googleapis.com/auth/drive"]

    creds = ServiceAccountCredentials.from_json_keyfile_name("../mpb_api.json", scope)

    client = gspread.authorize(creds)

    hoja_parametros = client.open("MPB").worksheet('Parametros')
    hoja_mva = client.open("MPB").worksheet('MVA')
    hoja_boot = client.open("MPB").worksheet('Bootstrap')
    hoja_fit = client.open("MPB").worksheet('Ajuste')

    fecha_sheet = hoja_mva.col_values(1)
    hora_sheet = hoja_mva.col_values(2)

    print('Acceso a la spreadsheet')
    # si ya está esa fecha en la spreadsheet, la sobreescribe. Si no, usa la siguiente línea vacía
    # pdb.set_trace()
    if date_entry in fecha_sheet and str(int(t1)) in hora_sheet:
        listaA = [a for a, fechas in enumerate(fecha_sheet) if fechas == date_entry]
        listaB = [b for b, horas in enumerate(hora_sheet) if horas == str(int(t1))]
        idx = list(set(listaA).intersection(listaB))[0]
        nr = idx+1
    else:
        nr = next_available_row(hoja_mva)

    # ######updatºla hoja de los parámetros
    hoja_parametros.update_acell(f'A{nr}', f'{date_entry}')
    hoja_parametros.update_acell(f'B{nr}', f'{int(t1)}')
    hoja_parametros.update_acell(f'D{nr}', f'{SZA:.3g}')
    hoja_parametros.update_acell(f'E{nr}', f'{int(altitud)}')
    hoja_parametros.update_acell(f'O{nr}', f'{round(B_norm_medio,2)}')

    cell_times = hoja_parametros.range(f'F{nr}:I{nr}')
    for i,cell in enumerate(cell_times):
        cell.value = tiempos[i]
    hoja_parametros.update_cells(cell_times)

    cell_B = hoja_parametros.range(f'L{nr}:N{nr}')
    for i,cell in enumerate(cell_B):
        cell.value = round(B_medio_vectorial[i],2)
    hoja_parametros.update_cells(cell_B)

    # if type(B_filtrado) != int:  # si no es un int, en particular 0, agrega los datos a la lista
    #     hoja_parametros.update_acell(f'J{nr}', f'{orden_filtro}')
    #     hoja_parametros.update_acell(f'K{nr}', f'{frec_filtro}')
    # else:
    hoja_parametros.update_acell(f'J{nr}', 'Sin filtrar')

    # #######update la hoja de MVA
    hoja_mva.update_acell(f'A{nr}', f'{date_entry}')
    hoja_mva.update_acell(f'B{nr}', f'{int(t1)}')
    hoja_mva.update_acell(f'D{nr}', f'{ti_MVA}')
    hoja_mva.update_acell(f'E{nr}', f'{tf_MVA}')
    hoja_mva.update_acell(f'I{nr}', f'{lamb[1]/lamb[2]:.3g}')

    hoja_mva.update_acell(f'S{nr}', f'{error_normal:.3g}')
    hoja_mva.update_acell(f'T{nr}', f'{round(np.mean(B3),2)}')
    hoja_mva.update_acell(f'U{nr}', f'{round(delta_B3,2)}')
    hoja_mva.update_acell(f'V{nr}', f'{abs(round(np.mean(B3)/B_norm_medio,2))}')

    cell_lambda = hoja_mva.range(f'F{nr}:H{nr}')
    for i,cell in enumerate(cell_lambda):
        cell.value = round(lamb[i],2)
    hoja_mva.update_cells(cell_lambda)

    cell_av = hoja_mva.range(f'J{nr}:R{nr}')
    for i,cell in enumerate(cell_av):
        cell.value = round(av[i],3)
    hoja_mva.update_cells(cell_av)

    # #######update la hoja de bootstrap
    hoja_boot.update_acell(f'A{nr}', f'{date_entry}')
    hoja_boot.update_acell(f'B{nr}', f'{int(t1)}')
    hoja_boot.update_acell(f'D{nr}', f'{M_cut}')
    hoja_boot.update_acell(f'E{nr}', f'{N_boot}')

    cell_normal = hoja_boot.range(f'F{nr}:H{nr}')
    for i,cell in enumerate(cell_normal):
        cell.value = round(normal_boot[i],3)
    hoja_boot.update_cells(cell_normal)

    hoja_boot.update_acell(f'I{nr}', f'{error_boot:.3g}')
    hoja_boot.update_acell(f'J{nr}', f'{round(np.mean(B3_boot),2)}')
    hoja_boot.update_acell(f'K{nr}', f'{round(sigmaB,2)}')
    hoja_boot.update_acell(f'L{nr}', f'{abs(round(np.mean(B3_boot)/B_norm_medio,2))}')

    # #######update la hoja del ajuste
    hoja_fit.update_acell(f'A{nr}', f'{date_entry}')
    hoja_fit.update_acell(f'B{nr}', f'{int(t1)}')
    hoja_fit.update_acell(f'D{nr}', f'{L0}')
    hoja_fit.update_acell(f'E{nr}', f'{e}')
    hoja_fit.update_acell(f'F{nr}', f'{x0}')

    cell_normal = hoja_fit.range(f'G{nr}:I{nr}')
    for i,cell in enumerate(cell_normal):
        cell.value = round(normal_fit[i],3)
    hoja_fit.update_cells(cell_normal)

    hoja_fit.update_acell(f'K{nr}', f'{round(np.mean(B3_fit),2)}')
    hoja_fit.update_acell(f'L{nr}', f'{abs(round(np.mean(B3_fit)/B_norm_medio,2))}')

    print('escribe la spreadsheet')
    return x3, normal_boot, normal_fit, t, B, posicion, inicio, fin, B_cut, t1, t2, t3, t4, B_medio_vectorial, nr
