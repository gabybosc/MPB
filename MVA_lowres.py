import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
from scipy.stats import norm
import datetime as dt
from funciones import error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, datenum, unix_to_decimal
from funciones_MVA import ajuste_conico, plot_velocidades, plot_FLorentz, plot_bootstrap
from funciones_plot import hodograma, set_axes_equal
import gspread
from oauth2client.service_account import ServiceAccountCredentials

"""
Para datos de MAg de baja resolución

Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o dia_del_año-año) y los cuatro tiempos t1 t2 t3 t4.
Eventualmente podría simplemente encontrar todos los cruces que quiero y decirle que lea directamente de algún lugar eso.
Antes de correr este programa hay que haber usado plot_seleccionar_puntos, para tener los cuatro tiempos elegidos.
Nos devuelve los lambdas, el cociente de lambdas, el omega, los autovectores y la normal para el MVA, el ajuste y el bootstrap.
Nos da el valor medio de B, de la altitud y el SZA. Grafica el hodograma.

Guarda los datos en una spreadsheet de google
"""


np.set_printoptions(precision=4)

# #si tengo la fecha en dia-mes-año
# date_entry = input('Enter a date in YYYY-MM-DD format \n')
# year, month, day = map(int, date_entry.split('-'))
# date_orbit = dt.date(year, month, day)

#si tengo la fecha en dia del año
# date_entry = input('Enter a date in YYYY-DDD format \n')
date_entry = '2016-076'
year, doy = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop
mag = np.loadtxt(path + 'MAG_1s/2016/mvn_mag_l2_{0}{3}ss1s_{0}{1}{2}_v01_r01.sts'.format(year, month, day, doy), skiprows=148) #datos MAG 1s (para plotear no quiero los datos pesados)
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

#si quiero elegir entre ciertas horas:

# t1 = float(input("t1 = "))
# t2 = float(input("t2 = "))
# t3 = float(input("t3 = "))
# t4 = float(input("t4 = "))
t1 = 18.2167
t2 = 18.2204
t3 = 18.235
t4 = 18.2476

inicio = np.where(t == find_nearest_inicial(t, t2))[0][0]
fin = np.where(t == find_nearest_final(t, t3))[0][0]

t_1 = np.where(t == find_nearest(t,t1))[0][0]
t_2 = np.where(t == find_nearest(t,t2))[0][0]
t_3 = np.where(t == find_nearest(t,t3))[0][0]
t_4 = np.where(t == find_nearest(t,t4))[0][0]
tiempos = [t1,t2,t3,t4]

#################

#ahora empieza el MVA con los datos que elegí
MD_cut = MD[inicio : fin+1, :]
M_cut = np.size(MD_cut[:,0])
B_cut = B[inicio:fin+1, :]
posicion_cut = posicion[inicio:fin+1,:]

M_ij = Mij(B_cut)
# np.savetxt('outs/matriz_%d'%dia[0], M_ij) #guarda la matriz de cada dia

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

#lambda2/lambda3
print('lambda1 = {0:1.3g} \nlambda2 = {1:1.3g} \nlambda3 = {2:1.3g}'.format(lamb[0], lamb[1], lamb[2]))
print('lambda2/lambda3 = {0:1.3g}'.format(lamb[1]/lamb[2]))
#las proyecciones
B1 = np.dot(B_cut, x1)
B2 = np.dot(B_cut, x2)
B3 = np.dot(B_cut, x3)


#el B medio
B_medio_vectorial = np.mean(B_cut, axis=0)
altitud = np.mean(MD_cut[:,8])
SZA = np.arccos(np.clip(np.dot(posicion_cut[0,:]/np.linalg.norm(posicion_cut[0,:]), [1,0,0]), -1.0, 1.0))* 180/np.pi

print('El valor medio de la altitud = {0:1.3g} km'.format(np.mean(MD_cut[:,8])))
print(f'El SZA es {SZA:1.3g}')
B_norm_medio = np.mean(np.linalg.norm(B_cut))
print('El valor medio de |B| es = {0:1.3g} nT'.format(B_norm_medio))

# hodograma(B1, B2, B3, 'nT', 'MAVEN MAG MVA ')

#el error
phi, delta_B3 = error(lamb, B_cut, M_cut, x)

###############
####fit
orbita = posicion[np.where(t == find_nearest_inicial(t, t1-1))[0][0] : np.where(t == find_nearest_final(t, t4+1))[0][0], :] / 3390 #radios marcianos

t_nave = find_nearest(t,(t2+t3)/2) #el tiempo en el medio de la hoja de corriente
index = np.where(t == t_nave)[0][0]
x0 = 0.78
e = 0.9
norm_vignes, X1, Y1, Z1, R = ajuste_conico(posicion, index, orbita, x3)



#######
out = np.zeros(1000)
out_phi = np.zeros((1000, 2))
normal_ran = np.zeros((1000, 3))
for a in range(1000):
    index = np.random.choice(B_cut.shape[0], M_cut, replace=True) #elije M índices de B, puede repetir (replace = True)
    B_random = B_cut[index,:] #me da un B hecho a partir de estos índices random

    Mij_random = Mij(B_random)

    [lamb_ran, x_ran] = np.linalg.eigh(Mij_random) #uso eigh porque es simetrica
    idx_ran = lamb_ran.argsort()[::-1]
    lamb_ran = lamb_ran[idx_ran]
    x_ran = x_ran[:,idx_ran]
    #ojo que a veces me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
    x3_ran = x_ran[:,2]
    if x3_ran[0] < 0:
        x3_ran = - x3_ran
    normal_ran[a, :] = x3_ran

    B3_ran = np.dot(B_random, x3_ran)
    phi, delta_B3 = error(lamb_ran, B_random, M_cut, x_ran)
    out[a] = np.mean(B3_ran)
    out_phi[a,0] = phi[2,0]
    out_phi[a,1] = phi[2,1]

normal_boot = np.mean(normal_ran, axis=0)
normal_boot = normal_boot/np.linalg.norm(normal_boot)

muB, sigmaB, mu31, sigma31, mu32, sigma32 = plot_bootstrap(out, out_phi)

print(f'La normal del MVA es = {x3}')
print(f'La normal de bootstrap es {normal_boot}')
print(f'La normal del ajuste es {norm_vignes}')

print('El valor medio de B a lo largo de la normal del MVA es B3 = {0:1.3g} nT'.format(np.mean(B3)))
print('El valor medio de B a lo largo de la normal del bootstrap es = {0:1.3g} nT'.format(np.mean(np.dot(B_cut, normal_boot))))
print('El valor medio de B a lo largo de la normal del ajuste es = {0:1.3g} nT'.format(np.mean(np.dot(B_cut, norm_vignes))))

print('|B3|/B_medio = ', np.abs(np.mean(B3)/np.mean(B_cut)))
print('|B3_boot|/B_medio = ', np.abs(np.mean(np.dot(B_cut, normal_boot))/np.mean(B_cut)))
print('|B3_ajuste|/B_medio = ', np.abs(np.mean(np.dot(B_cut, norm_vignes))/np.mean(B_cut)))

print('Matriz de incerteza angular (grados): \n{}'.format(phi  *  180 / np.pi))
print('<B3> = {0:1.3g} +- {1:1.3g} nT'.format(np.mean(B3),delta_B3))


if phi[2,1] > phi[2,0]:
    print('El error phi32 = {0:1.3g}º es mayor a phi31 = {1:1.3g}º'.format(phi[2,1]*57.2958, phi[2,0]*180 / np.pi))
else:
    print('El error phi31 = {0:1.3g}º es mayor a phi32 = {1:1.3g}º'.format(phi[2,0]*57.2958, phi[2,1]*180 / np.pi))
    #quiero ver si el error más grande es phi31 o phi32

print('mu_boot = {0:1.3g} nT, std_boot={1:1.3g} nT'.format(muB, sigmaB))
print('mean_phi31 = {0:1.3g}º, std_phi31={1:1.3g}º'.format(mu31, sigma31))
print('mean_phi32 = {0:1.3g}º, std_phi32={1:1.3g}º'.format(mu32, sigma32))

#############
#ahora guardo todo en una spreadsheet
scope = ["https://spreadsheets.google.com/feeds","https://www.googleapis.com/auth/spreadsheets","https://www.googleapis.com/auth/drive.file","https://www.googleapis.com/auth/drive"]

creds = ServiceAccountCredentials.from_json_keyfile_name("google_api_MPB.json", scope)

client = gspread.authorize(creds)

sheet = client.open("MPB").sheet1

def next_available_row(sheet):
    str_list = list(filter(None, sheet.col_values(1)))
    return str(len(str_list)+1)

nr = next_available_row(sheet)

cell_times = sheet.range(f'D{nr}:G{nr}')
cell_lambda = sheet.range(f'H{nr}:J{nr}')
cell_av = sheet.range(f'K{nr}:S{nr}')

for i,cell in enumerate(cell_times):
    cell.value = tiempos[i]

for i,cell in enumerate(cell_lambda):
    cell.value = lamb[i]

for i,cell in enumerate(cell_av):
    cell.value = av[i]


sheet.update_acell(f'A{nr}', f'{doy}-{year}')
sheet.update_acell(f'B{nr}', f'{int(t1)}')
sheet.update_cells(cell_times)
sheet.update_cells(cell_lambda)
sheet.update_cells(cell_av)
