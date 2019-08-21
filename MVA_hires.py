import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
import scipy.signal as signal
import datetime as dt
from scipy.stats import norm
from funciones import error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, next_available_row, datenum, unix_to_decimal
from funciones_MVA import ajuste_conico, plot_velocidades, plot_FLorentz, plot_bootstrap, bootstrap
from funciones_plot import hodograma, set_axes_equal
import gspread
from oauth2client.service_account import ServiceAccountCredentials


"""
DEBUGGEAR MIRANDO EL LOWRES

Para datos de mag de alta resolución. Para que tarde menos en cargar, skipea las primeras t1*32*3600 rows, lo malo de eso es que a veces no alcanza con esos datos.

Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o doy-año) y los cuatro tiempos t1 t2 t3 t4.
Eventualmente podría simplemente encontrar todos los cruces que quiero y decirle que lea directamente de algún lugar eso. (es lo que hace MVA_automatico)
Antes de correr este programa hay que haber usado plot_seleccionar_puntos, para tener los cuatro tiempos elegidos.
Toma los datos de alta resolución y les aplica un filtro pasabajos con ventana Butterworth con frecuencia de corte de 0.1 Hz de orden 3.
A los datos filtrados les aplica el MVA.
Nos devuelve los lambdas, el cociente de lambdas, el omega, los autovectores y la normal.
Nos da el valor medio de B, de la altitud y el SZA.
Devuelve el ancho de la MPB y la corriente que pasa. Calcula también la fuerza de lorentz y el campo de hall.
Grafica el hodograma, el ajuste de vignes, y la comparación de las normales obtenidas por los distintos métodos.
"""


np.set_printoptions(precision=4)

# #si tengo la fecha en dia-mes-año
date_entry = '2016-076'#input('Enter a date in YYYY-DDD format \n')
year, doy = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

t1 = 18.2167
t2 = 18.2204
t3 = 18.235
t4 = 18.2476

# t1 = float(input("t1 = "))
# t2 = float(input("t2 = "))
# while t2 < t1:
#     print('t2 no puede ser menor a t1. \n')
#     t2 = float(input('t2 = '))
# t3 = float(input("t3 = "))
# while t3 < t2:
#     print('t3 no puede ser menor a t2. \n')
#     t3 = float(input('t3 = '))
# t4 = float(input("t4 = "))
# while t4 < t3:
#     print('t4 no puede ser menor a t3. \n')
#     t4 = float(input('t4 = '))



# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop
n = int(t1*32*3600 - 1000)
mag = np.genfromtxt(path + f'MAG_hires/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skip_header=n, skip_footer=1000) #skipea las filas anteriores a la hora que quiero medir para hacer más rápido
lpw = cdf.CDF(path + f'LPW/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf')
#para ver las varaibles del cdf:
# lpw.cdf_info()


dia = mag[:,1]
t = mag[:,6]  #el dia decimal
t = (t - dia) * 24 #para que me de sobre la cantidad de horas

M = np.size(t) #el numero de datos

#tengo que asegurarme de que no haya agujeros en mis datos
for i in range(M-1):
    if t[i+1] - t[i] > 24 * 1.5e-5: #1.1e-5 es 1s, le doy un poco más
        print('salto en la linea {} de {} segundos'.format(i+144, (t[i+1] - t[i]) / (24*1.1e-5)))

#el campo
B = np.zeros((M, 3))
for i in range(7,10):
    B[:,i-7] = mag[:, i]

#la posición(x,y,z)
posicion = np.zeros((M, 3))
for i in range(11,14):
    posicion[:,i-11] = mag[:, i]

#la matriz diaria:
MD = np.zeros((M, 9))
MD[:, 0] = t
for i in range(1,4):
    MD[:, i] = B[:,i-1]
MD[:,4] = np.linalg.norm(B, axis=1)#la norma de B
for i in range(5,8):
    MD[:,i] = posicion[:,i-5]/3390 #en radios marcianos
MD[:, 8] = np.linalg.norm(posicion, axis=1) - 3390 #altitud en km


inicio = np.where(t == find_nearest_inicial(t, t2))[0][0]
fin = np.where(t == find_nearest_final(t, t3))[0][0]

t_1 = np.where(t == find_nearest(t,t1))[0][0]
t_2 = np.where(t == find_nearest(t,t2))[0][0]
t_3 = np.where(t == find_nearest(t,t3))[0][0]
t_4 = np.where(t == find_nearest(t,t4))[0][0]
tiempos = [t1,t2,t3,t4]

#################
#Filtramos
b,a = signal.butter(3,0.01,btype='lowpass')
Bx_filtrado = signal.filtfilt(b, a, B[:,0])
By_filtrado = signal.filtfilt(b, a, B[:,1])
Bz_filtrado = signal.filtfilt(b, a, B[:,2])
B_filtrado = np.linalg.norm(np.array([Bx_filtrado, By_filtrado, Bz_filtrado]), axis=0) #es el axis 0 porque no está traspuesta


#ahora empieza el MVA con los datos que elegí
MD_cut = MD[inicio : fin+1, :]
M_cut = np.size(MD_cut[:,0])
B_cut = np.transpose(np.array([Bx_filtrado[inicio:fin+1], By_filtrado[inicio:fin+1], Bz_filtrado[inicio:fin+1]]))
posicion_cut = posicion[inicio:fin+1,:]

M_ij = Mij(B_cut)

#ahora quiero los autovectores y autovalores
[lamb, x] = np.linalg.eigh(M_ij) #uso eigh porque es simetrica

#Los ordeno de mayor a menor
idx = lamb.argsort()[::-1]
lamb = lamb[idx]
x = x[:,idx]
#ojo que me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
x1 = x[:,0]
x2 = x[:,1]
x3 = x[:,2]

av = np.concatenate([x1,x2,x3])

if x3[0] < 0: #si la normal aputna para adentro me la da vuelta
    x3 = - x3
if any(np.cross(x1,x2) - x3) > 0.01:
    print('Cambio el signo de x1 para que los av formen terna derecha')
    x1 = -x1

#las proyecciones
B1 = np.dot(B_cut, x1)
B2 = np.dot(B_cut, x2)
B3 = np.dot(B_cut, x3)


#el B medio
B_medio_vectorial = np.mean(B_cut, axis=0)
altitud = np.mean(MD_cut[:,8])
SZA = np.arccos(np.clip(np.dot(posicion_cut[0,:]/np.linalg.norm(posicion_cut[0,:]), [1,0,0]), -1.0, 1.0))* 180/np.pi

B_norm_medio = np.linalg.norm(B_medio_vectorial)

# hodograma(B1, B2, B3, 'nT', 'MAVEN MAG MVA')

#el error
phi, delta_B3 = error(lamb, B_cut, M_cut, x)

###############
####fit
orbita = posicion[np.where(t == find_nearest_inicial(t, t1-1))[0][0] : np.where(t == find_nearest_final(t, t4+1))[0][0], :] / 3390 #radios marcianos

t_nave = find_nearest(t,(t2+t3)/2) #el tiempo en el medio de la hoja de corriente
index = np.where(t == t_nave)[0][0]
x0 = 0.78
e = 0.9
normal_fit, X1, Y1, Z1, R, L0 = ajuste_conico(posicion, index, orbita, x3)

B3_fit = np.dot(B_cut, normal_fit)

###############
##Bootstrap

N_boot = 1000
normal_boot, phi_boot, delta_B3_boot, out, out_phi = bootstrap(N_boot, B_cut, M_cut)

muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

B3_boot = np.dot(B_cut, normal_boot)

#######
#Errores
if phi[2,1] > phi[2,0]:
    error_normal = phi[2,1]*57.2958
else:
    error_normal = phi[2,0]*57.2958
    #quiero ver si el error más grande es phi31 o phi32

if sigma31 > sigma32:
    error_boot = sigma31
else:
    error_boot = sigma32

#############
#ahora guardo todo en una spreadsheet
scope = ["https://spreadsheets.google.com/feeds","https://www.googleapis.com/auth/spreadsheets","https://www.googleapis.com/auth/drive.file","https://www.googleapis.com/auth/drive"]

creds = ServiceAccountCredentials.from_json_keyfile_name("mpb_api.json", scope)

client = gspread.authorize(creds)

hoja_parametros = client.open("MPB").worksheet('Parametros')
hoja_mva = client.open("MPB").worksheet('MVA')
hoja_boot = client.open("MPB").worksheet('Bootstrap')
hoja_fit = client.open("MPB").worksheet('Ajuste')

fecha_sheet = hoja_mva.col_values(1)
hora_sheet = hoja_mva.col_values(2)

#si ya está esa fecha en la spreadsheet, la sobreescribe. Si no, usa la siguiente línea vacía
if date_entry in fecha_sheet and str(int(t1)) in hora_sheet:
    listaA = [a for a, fechas in enumerate(fecha_sheet) if fechas == date_entry]
    listaB = [b for b, horas in enumerate(hora_sheet) if horas == str(int(t1))]
    idx = list(set(listaA).intersection(listaB))[0]
    nr = idx+1
else:
    nr = next_available_row(hoja_mva)

#######updatºla hoja de los parámetros
hoja_parametros.update_acell(f'A{nr}', f'{date_entry}')
hoja_parametros.update_acell(f'B{nr}', f'{int(t1)}')
hoja_parametros.update_acell(f'D{nr}', f'{SZA:.3g}')
hoja_parametros.update_acell(f'E{nr}', f'{int(altitud)}')
hoja_parametros.update_acell(f'J{nr}', f'{round(B_norm_medio,2)}')

cell_B = hoja_parametros.range(f'G{nr}:I{nr}')
for i,cell in enumerate(cell_B):
    cell.value = round(B_medio_vectorial[i],2)
hoja_parametros.update_cells(cell_B)


########update la hoja de MVA
hoja_mva.update_acell(f'A{nr}', f'{date_entry}')
hoja_mva.update_acell(f'B{nr}', f'{int(t1)}')
hoja_mva.update_acell(f'K{nr}', f'{lamb[1]/lamb[2]:.3g}')
hoja_mva.update_acell(f'U{nr}', f'{error_normal:.3g}')
hoja_mva.update_acell(f'V{nr}', f'{round(np.mean(B3),2)}')
hoja_mva.update_acell(f'W{nr}', f'{round(delta_B3,2)}')
hoja_mva.update_acell(f'X{nr}', f'{abs(round(np.mean(B3)/B_norm_medio,2))}')



cell_times = hoja_mva.range(f'D{nr}:G{nr}')
for i,cell in enumerate(cell_times):
    cell.value = tiempos[i]
hoja_mva.update_cells(cell_times)

cell_lambda = hoja_mva.range(f'H{nr}:J{nr}')
for i,cell in enumerate(cell_lambda):
    cell.value = round(lamb[i],2)
hoja_mva.update_cells(cell_lambda)

cell_av = hoja_mva.range(f'L{nr}:T{nr}')
for i,cell in enumerate(cell_av):
    cell.value = round(av[i],3)
hoja_mva.update_cells(cell_av)


########update la hoja de bootstrap
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

########update la hoja del ajuste
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
