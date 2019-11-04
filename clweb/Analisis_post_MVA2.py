import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
import scipy.signal as signal
import datetime as dt
import pdb
from matplotlib.mlab import normpdf
from scipy.stats import norm
from funciones import error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, next_available_row, datenum, unix_to_decimal, UTC_to_hdec
from funciones_MVA import ajuste_conico, plot_velocidades, plot_FLorentz, plot_bootstrap, bootstrap
from funciones_plot import hodograma, set_axes_equal
import gspread
from oauth2client.service_account import ServiceAccountCredentials

"""
Hace el MVA
calcula el ángulo entre normales
calcula el ancho de la mpb
calcula la corriente
"""
from MVA_hires import MVA


date_entry = input('Enter a date in YYYY-DDD or YYYY-MM-DD format \n')\

if len(date_entry.split('-')) < 3:
    year, doy = map(int, date_entry.split('-'))
    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date
else:
    year, month, day = map(int, date_entry.split('-'))
    date_orbit = dt.date(year, month, day)

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doy = date_orbit.strftime("%j")

# ti_MVA = float(input("t_incial = "))
# tf_MVA = float(input("t_final = "))
ti_MVA = UTC_to_hdec(input('Tiempo inicial hh:mm:ss\n'))
tf_MVA = UTC_to_hdec(input('Tiempo final hh:mm:ss\n'))
while tf_MVA < ti_MVA:
    print('t final no puede ser menor a t inicial. \n')
    # ti_MVA = float(input('t_inicial = '))
    # tf_MVA = float(input("t_final = "))
    ti_MVA = UTC_to_hdec(input('Tiempo inicial hh:mm:ss\n'))
    tf_MVA = UTC_to_hdec(input('Tiempo final hh:mm:ss\n'))


path = f'../../../datos/clweb/{year}-{month}-{day}/' #path a los datos desde la laptop
mag = np.loadtxt(path + 'mag.asc')
swia = np.loadtxt(path + 'swia_density.asc')
lpw = np.loadtxt(path + 'lpw_density.asc')

x3, normal_boot, normal_fit, t, B, posicion, inicio, fin, B_cut, t1, t2, t3, t4,B_medio_vectorial, nr = MVA(date_entry, ti_MVA, tf_MVA, mag)

#########
#buscamos el ángulo entre las normales
angulo_mva = np.arccos(np.clip(np.dot(normal_fit, x3), -1.0, 1.0)) #el clip hace que si por algun motivo el dot me da >1 (i.e. 1,00002), me lo convierte en 1
angulo_boot = np.arccos(np.clip(np.dot(normal_boot, x3), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!


##############
#Calculo la velocidad de la nave
v_punto = np.zeros((fin-inicio, 3))
norma_v = np.zeros(fin-inicio)
posicion_cut = posicion[inicio : fin+1, :]
t_cut = t[inicio : fin+1] * 3600 #en segundos
for i in range(fin-inicio):
    v_punto[i,:] = (posicion[inicio+1,:] - posicion[inicio]) / (t_cut[i+1]-t_cut[i]) #en km/s
    norma_v[i] = np.linalg.norm(v_punto[i,:])
#veamos que no cambia mucho punto a punto, usemos la norma
diff = max(norma_v)- min(norma_v)
#la velocidad promedio
v_media = np.array([np.mean(v_punto[:,0]), np.mean(v_punto[:,1]), np.mean(v_punto[:,2])])

#ahora quiero ver si la nave atraviesa perpendicularmente a la MPB
v_media_norm = v_media/np.linalg.norm(v_media) #la normalizo por cuestion de  cuentas nomas
angulo_v_fit = np.arccos(np.clip(np.dot(normal_fit, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_v_mva = np.arccos(np.clip(np.dot(x3, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_v_boot = np.arccos(np.clip(np.dot(normal_boot, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!

B_intermedio = B_medio_vectorial/np.linalg.norm(B_medio_vectorial) #la normalizo por cuestion de  cuentas nomas
angulo_B_fit = np.arccos(np.clip(np.dot(normal_fit, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_B_mva = np.arccos(np.clip(np.dot(x3, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_B_boot = np.arccos(np.clip(np.dot(normal_boot, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!

######
##Espesor de la MPB
#ahora veamos v_para
deltat_14 = (t4 - t1) * 3600
deltat_23 = (t3 - t2) * 3600

v_para = np.dot(v_media, normal_fit) * normal_fit
x_14_fit = v_para * deltat_14 #en km# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x_23_fit = v_para * deltat_23

#si ahora proyecto sobre la normal de la MVA
v_para_MVA = np.dot(v_media, x3) * x3
x_14_MVA = v_para_MVA * deltat_14 #en km
x_23_MVA = v_para_MVA * deltat_23

#si ahora proyecto sobre la normal del bootstrap
v_para_boot = np.dot(np.array([-2.487,  0.479,  2.836]), normal_boot) * normal_boot
x_14_boot = v_para_boot * deltat_14 #en km
x_23_boot = v_para_boot * deltat_23

# plot_velocidades(X1, Y1, Z1, R, normal_fit, x3, v_media, v_para, v_para_MVA)

###########
#giroradio


#########
###análisis de corrientes

inicio_up = np.where(t == find_nearest_inicial(t, t1-0.015))[0][0] #las 18:12:00
fin_up = np.where(t == find_nearest_final(t, t1))[0][0] #las 18:13:00
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0) #nT

inicio_down = np.where(t == find_nearest_inicial(t, t4))[0][0] #las 18:14:51
fin_down = np.where(t == find_nearest_final(t, t4+0.015))[0][0] #las 18:15:52
B_downstream = np.mean(B[inicio_down:fin_down,:], axis=0) #nT



omega = np.arccos(np.dot(B_upstream,B_downstream)/(np.linalg.norm(B_upstream)*np.linalg.norm(B_downstream)))

mu = 4* np.pi * 1E-7 #Henry/m

J_s_MVA = np.cross(x3, (B_upstream-B_downstream)) / mu #nA/m
ancho_mpb = np.linalg.norm(x_23_MVA) #considero que es el tiempo corto a lo largo de la normal del MVA
J_v_MVA = J_s_MVA / (1000*ancho_mpb) #nA/m²

J_s_fit = np.cross(normal_fit, (B_upstream-B_downstream)) / mu #nA/m
ancho_mpb_fit = np.linalg.norm(x_23_fit) #considero que es el tiempo corto a lo largo de la normal del ajuste
J_v_fit = J_s_fit / (1000*ancho_mpb) #nA/m²

J_s_boot = np.cross(normal_boot, (B_upstream-B_downstream)) / mu #nA/m
ancho_mpb_boot = np.linalg.norm(x_23_boot) #km
J_v_boot = J_s_boot / (1000*ancho_mpb_boot) #nA/m²

fuerza_mva = np.cross(J_v_MVA * 1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
fuerza_fit = np.cross(J_v_fit *1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
fuerza_boot = np.cross(J_v_boot * 1E-9, B[inicio_down,:]*1E-9) #N/m^3


e_density = lpw.varget('data')[:,3]
t_unix = lpw.varget('time_unix')
t_lpw = unix_to_decimal(t_unix)
ti_lpw = np.where(t_lpw == find_nearest(t_lpw, t2))[0][0]
tf_lpw = np.where(t_lpw == find_nearest(t_lpw, t3))[0][0]
n_e = np.nanmean(e_density[ti_lpw:tf_lpw]) #hace el mean ignorando los nans #cm⁻³
n_e = n_e * 1E6 #m⁻³
if np.isnan(n_e):
    n_e = 1E7
    print('LPW no tiene datos de densidad, asumí n_e = 1E7')
q_e = 1.6E-19 #carga electron #C

E_Hall = np.cross(J_v_MVA * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
E_Hall_fit = np.cross(J_v_fit * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
E_Hall_boot = np.cross(J_v_boot * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m


#
# plot_FLorentz(X1, Y1, Z1, R, J_v, B_upstream, B_downstream, fuerza_mva, x3)
#
# plt.figure()
# plt.plot(t[inicio_up:fin_down], fuerza_mva)
# plt.plot(t[inicio_up:fin_down], fuerza_ajuste)
# plt.legend(['Fx', 'Fy', 'Fz'])
# plt.xlabel('Tiempo')
# plt.ylabel('Fuerza (N/m^3)')
#
plt.show(block=False)


############
#Ahora guardamos todo en la spreadsheet
scope = ["https://spreadsheets.google.com/feeds","https://www.googleapis.com/auth/spreadsheets","https://www.googleapis.com/auth/drive.file","https://www.googleapis.com/auth/drive"]

creds = ServiceAccountCredentials.from_json_keyfile_name("mpb_api.json", scope)

client = gspread.authorize(creds)

hoja_parametros = client.open("MPB").worksheet('Parametros')
hoja_mva = client.open("MPB").worksheet('MVA')
hoja_boot = client.open("MPB").worksheet('Bootstrap')
hoja_fit = client.open("MPB").worksheet('Ajuste')


##########
#Parámetros

hoja_parametros.update_acell(f'Z{nr}', f'{omega*180/np.pi:.3g}')
hoja_parametros.update_acell(f'S{nr}', f'{np.linalg.norm(v_media):.3g}')

cell_vel = hoja_parametros.range(f'P{nr}:R{nr}')
for i,cell in enumerate(cell_vel):
    cell.value = round(v_media[i],2)
hoja_parametros.update_cells(cell_vel)

cell_Bup = hoja_parametros.range(f'T{nr}:V{nr}')
for i,cell in enumerate(cell_Bup):
    cell.value = round(B_upstream[i],2)
hoja_parametros.update_cells(cell_Bup)

cell_Bdown = hoja_parametros.range(f'W{nr}:Y{nr}')
for i,cell in enumerate(cell_Bdown):
    cell.value = round(B_downstream[i],2)
hoja_parametros.update_cells(cell_Bdown)


#La hoja del MVA

hoja_mva.update_acell(f'W{nr}', f'{angulo_v_mva * 180/np.pi:.3g}')
hoja_mva.update_acell(f'X{nr}', f'{angulo_B_mva * 180/np.pi:.3g}')
hoja_mva.update_acell(f'Y{nr}', f'{np.linalg.norm(x_23_MVA):.3g}')
hoja_mva.update_acell(f'Z{nr}', f'{np.linalg.norm(x_14_MVA):.3g}')
# hoja_mva.update_acell(f'AA{nr}', f'{:.3g}')
# hoja_mva.update_acell(f'AB{nr}', f'{:.3g}')

hoja_mva.update_acell(f'AF{nr}', f'{np.linalg.norm(J_s_MVA)*1E-6:.3g}')
hoja_mva.update_acell(f'AJ{nr}', f'{np.linalg.norm(J_v_MVA):.3g}')
hoja_mva.update_acell(f'AK{nr}', f'{np.linalg.norm(fuerza_mva):.3g}')
hoja_mva.update_acell(f'AO{nr}', f'{np.linalg.norm(E_Hall)*1E3:.3g}') #mV/m

cell_Js = hoja_mva.range(f'AC{nr}:AE{nr}')
for i,cell in enumerate(cell_Js):
    cell.value = round(J_s_MVA[i] * 1E-6,3)
hoja_mva.update_cells(cell_Js)

cell_Jv = hoja_mva.range(f'AG{nr}:AI{nr}')
for i,cell in enumerate(cell_Jv):
    cell.value = round(J_v_MVA[i],3)
hoja_mva.update_cells(cell_Jv)

cell_EH = hoja_mva.range(f'AL{nr}:AN{nr}')
for i,cell in enumerate(cell_EH):
    cell.value = round(E_Hall[i]*1E3,3)
hoja_mva.update_cells(cell_EH)

#La hoja del bootstrap

hoja_boot.update_acell(f'M{nr}', f'{angulo_v_boot * 180/np.pi:.3g}')
hoja_boot.update_acell(f'N{nr}', f'{angulo_B_boot * 180/np.pi:.3g}')
hoja_boot.update_acell(f'O{nr}', f'{np.linalg.norm(x_23_boot):.3g}')
hoja_boot.update_acell(f'P{nr}', f'{np.linalg.norm(x_14_boot):.3g}')
# hoja_boot.update_acell(f'Q{nr}', f'{:.3g}')
# hoja_boot.update_acell(f'R{nr}', f'{:.3g}')

hoja_boot.update_acell(f'V{nr}', f'{np.linalg.norm(J_s_boot)*1E-6:.3g}')
hoja_boot.update_acell(f'Z{nr}', f'{np.linalg.norm(J_v_boot):.3g}')
hoja_boot.update_acell(f'AA{nr}', f'{np.linalg.norm(fuerza_boot):.3g}')
hoja_boot.update_acell(f'AE{nr}', f'{np.linalg.norm(E_Hall_boot)*1E3:.3g}')

cell_Js = hoja_boot.range(f'S{nr}:U{nr}')
for i,cell in enumerate(cell_Js):
    cell.value = round(J_s_boot[i] * 1E-6,3)
hoja_boot.update_cells(cell_Js)

cell_Jv = hoja_boot.range(f'W{nr}:Y{nr}')
for i,cell in enumerate(cell_Jv):
    cell.value = round(J_v_boot[i],3)
hoja_boot.update_cells(cell_Jv)

cell_EH = hoja_boot.range(f'AB{nr}:AD{nr}')
for i,cell in enumerate(cell_EH):
    cell.value = round(E_Hall_boot[i]*1E3,3)
hoja_boot.update_cells(cell_EH)


#La hoja del ajuste

hoja_fit.update_acell(f'J{nr}', f'{angulo_mva * 180/np.pi:.3g}')
hoja_fit.update_acell(f'M{nr}', f'{angulo_v_fit * 180/np.pi:.3g}')
hoja_fit.update_acell(f'N{nr}', f'{angulo_B_fit * 180/np.pi:.3g}')
hoja_fit.update_acell(f'O{nr}', f'{np.linalg.norm(x_23_fit):.3g}')
hoja_fit.update_acell(f'P{nr}', f'{np.linalg.norm(x_14_fit):.3g}')
# hoja_fit.update_acell(f'Q{nr}', f'{:.3g}')
# hoja_fit.update_acell(f'R{nr}', f'{:.3g}')

hoja_fit.update_acell(f'V{nr}', f'{np.linalg.norm(J_s_fit)*1E-6:.3g}')
hoja_fit.update_acell(f'Z{nr}', f'{np.linalg.norm(J_v_fit):.3g}')
hoja_fit.update_acell(f'AA{nr}', f'{np.linalg.norm(fuerza_fit):.3g}')
hoja_fit.update_acell(f'AE{nr}', f'{np.linalg.norm(E_Hall_fit)*1E3:.3g}')

cell_Js = hoja_fit.range(f'S{nr}:U{nr}')
for i,cell in enumerate(cell_Js):
    cell.value = round(J_s_fit[i] * 1E-6,3)
hoja_fit.update_cells(cell_Js)

cell_Jv = hoja_fit.range(f'W{nr}:Y{nr}')
for i,cell in enumerate(cell_Jv):
    cell.value = round(J_v_fit[i],3)
hoja_fit.update_cells(cell_Jv)

cell_EH = hoja_fit.range(f'AB{nr}:AD{nr}')
for i,cell in enumerate(cell_EH):
    cell.value = round(E_Hall_fit[i]*1E3,3)
hoja_fit.update_cells(cell_EH)
