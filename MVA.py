from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from scipy.stats import norm
import datetime as dt
import scipy.signal as signal
from funciones import hodograma, error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij, set_axes_equal,datenum

"""
Este script pide como user input una fecha (puede ser o fecha dd-mm-aaaa o dia_del_año-año) y los cuatro tiempos t1 t2 t3 t4.
Eventualmente podría simplemente encontrar todos los cruces que quiero y decirle que lea directamente de algún lugar eso.
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
# date_entry = input('Enter a date in YYYY-MM-DD format \n')
# year, month, day = map(int, date_entry.split('-'))
# date_orbit = dt.date(year, month, day)

#si tengo la fecha en dia del año
date_entry = input('Enter a date in YYYY-DDD format \n')
year, doty = map(int, date_entry.split('-'))
date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doty - 1) #para convertir el doty en date

year = date_orbit.strftime("%Y")
month = date_orbit.strftime("%m")
day = date_orbit.strftime("%d")
doty = date_orbit.strftime("%j")

# path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
path = '../../datos/' #path a los datos desde la laptop
mag = np.loadtxt(path + 'MAG_1s/mvn_mag_l2_{0}{3}ss1s_{0}{1}{2}_v01_r01.sts'.format(year, month, day, doty), skiprows=148) #datos MAG 1s (para plotear no quiero los datos pesados)
lpw = cdf.CDF(path + 'mvn_lpw_l2_lpnt_{0}{1}{2}_v03_r02.cdf'.format(year, month, day, doty))

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

t1 = float(input("t1 = "))
t2 = float(input("t2 = "))
t3 = float(input("t3 = "))
t4 = float(input("t4 = "))

inicio = np.where(t == find_nearest_inicial(t, t2))[0][0]
fin = np.where(t == find_nearest_final(t, t3))[0][0]

t_1 = np.where(t == find_nearest(t,t1))[0][0]
t_2 = np.where(t == find_nearest(t,t2))[0][0]
t_3 = np.where(t == find_nearest(t,t3))[0][0]
t_4 = np.where(t == find_nearest(t,t4))[0][0]


#################
#Filtramos los cantidad_datos
b,a = signal.butter(3,0.1,btype='lowpass')
Bx_filtrado = signal.filtfilt(b, a, B[:,0])
By_filtrado = signal.filtfilt(b, a, B[:,1])
Bz_filtrado = signal.filtfilt(b, a, B[:,2])
B_filtrado = np.linalg.norm(np.array([Bx_filtrado, By_filtrado, Bz_filtrado]), axis=0) #es el axis 0 porque no está traspuesta


#ahora empieza el MVA con los datos que elegí
MD_cut = MD[inicio : fin+1, :]
M_cut = np.size(MD_cut[:,0])
B_cut = np.transpose(np.array([Bx_filtrado[inicio:fin+1], By_filtrado[inicio:fin+1], Bz_filtrado[inicio:fin+1]]))

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
if x3[0] < 0: #si la normal aputna para adentro me la da vuelta
    x3 = - x3

#lambda2/lambda3
print('lambda1 = {0:1.3g} \nlambda2 = {1:1.3g} \nlambda3 = {2:1.3g}'.format(lamb[0], lamb[1], lamb[2]))
print('lambda2/lambda3 = {0:1.3g}'.format(lamb[1]/lamb[2]))
print('La normal es = {}'.format(x3))
#las proyecciones
B1 = np.dot(B_cut, x1)
B2 = np.dot(B_cut, x2)
B3 = np.dot(B_cut, x3)

hodograma(B1, B2, B3, 'nT', 'MAVEN MAG MVA ')

#el error
phi, delta_B3 = error(lamb, B_cut, M_cut, x, dia)
print('Matriz de incerteza angular (grados): \n{}'.format(phi  *  180 / np.pi))
print('<B3> = {0:1.3g} +- {1:1.3g} nT'.format(np.mean(B3),delta_B3))


if phi[2,1] > phi[2,0]:
    print('El error phi32 = {0:1.3g}º es mayor a phi31 = {1:1.3g}º'.format(phi[2,1]*57.2958, phi[2,0]*180 / np.pi))
else:
    print('El error phi31 = {0:1.3g}º es mayor a phi32 = {1:1.3g}º'.format(phi[2,0]*57.2958, phi[2,1]*180 / np.pi))
#quiero ver si el error más grande es phi31 o phi32

#el B medio
B_medio_vectorial = np.mean(B_cut, axis=0)
print('El valor medio de B = {0:1.3g} nT'.format(np.linalg.norm(B_medio_vectorial)))#np.mean(B_cut, 0))) para que me de los tres valores medios
print('|B3|/B_medio = ', np.abs(np.mean(B3)/np.linalg.norm(B_medio_vectorial)))
print('El valor medio de la altitud = {0:1.3g} km'.format(np.mean(MD_cut[:,8])))

# x3 = np.array([ 0.9182, -0.3183,  0.2357])
###############
orbita = posicion[np.where(t == find_nearest_inicial(t, t1-1))[0][0] : np.where(t == find_nearest_final(t, t4+1))[0][0], :] / 3390 #radios marcianos

# usamos vignes et al:
x0 = 0.78
e = 0.9
L = 0.96

#ec conica
theta = np.linspace(0, np.pi *3/4, 100)
phi = np.linspace(0, 2 * np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(THETA))

####dibujamos el punto por el que pasa la nave:
t_nave = find_nearest(t,(t2+t3)/2) #el tiempo en el medio de la hoja de corriente
index = np.where(t == t_nave)[0][0]
R = posicion[index,:] / 3390 #la posicion de la nave en RM


####### Calculo mi propia elipse que pase por el punto.
r0 = R - np.array([x0,0,0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(THETA))


fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.set_xlabel(r'$X_{MSO} (R_m)$')
ax.set_ylabel(r'$Y_{MSO} (R_m)$')
ax.set_zlabel(r'$Z_{MSO} (R_m)$')
ax.set_aspect('equal')
ax.plot(orbita[:,0], orbita[:,1], orbita[:,2], color='green', label='Órbita')
ax.scatter(R[0], R[1], R[2], label='MAVEN', color='k', s=40)#, marker='x')
X1 = x0 + r1 * np.cos(THETA)
Y1 = r1 * np.sin(THETA) * np.cos(PHI)
Z1 = r1 * np.sin(THETA) * np.sin(PHI)
plot = ax.plot_surface(
    X1, Y1, Z1, rstride=4, cstride=4, alpha=0.5, edgecolor='none', cmap=plt.get_cmap('Blues_r'))
asc = L0 / (1 - e**2)  #semieje mayor
bsc = np.sqrt(asc*L0)
csc = e*asc - x0 #donde está centrada. Hay que ver el signo

norm_vignes = np.array([(R[0]+csc)*2 / asc**2, R[1]*2 / (bsc)**2, R[2]*2 / (bsc)**2]) #la normal de vignes
norm_vignes = norm_vignes/np.linalg.norm(norm_vignes) #normalizado
ax.quiver(R[0], R[1], R[2], norm_vignes[0], norm_vignes[1], norm_vignes[2], color='b',length=0.5, label='Normal del ajuste') #asi se plotea un vector
ax.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color='k',length=0.5, label='Normal del MVA')
# normal_boot = np.array([0.9183, 0.3186, 0.2351])
# ax.quiver(R[0], R[1], R[2], normal_boot[0], normal_boot[1], normal_boot[2], color='m',length=0.5, label='Normal del bootstrap')

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)
ax.legend()
set_axes_equal(ax) #para que tenga forma de esfera la esfera
# plt.show()

print('la normal del ajuste es {}'.format(norm_vignes))
print('El valor medio de B a lo largo de la normal del ajuste es = {0:1.3g} nT'.format(np.mean(np.dot(B_cut, norm_vignes))))#np.mean(B_cut, 0))) para que me de los tres valores medios
print('|B3_ajuste|/B_medio = ', np.abs(np.mean(np.dot(B_cut, norm_vignes))/np.mean(B_cut)))

#########
#buscamos el ángulo entre ambas normales
angulo_mva = np.arccos(np.clip(np.dot(norm_vignes, x3), -1.0, 1.0)) #el clip hace que si por algun motivo el dot me da >1 (i.e. 1,00002), me lo convierte en 1
print('El ángulo entre la normal del ajuste y la del MVA = {0:1.3g}º'.format(angulo_mva * 180/np.pi))

angulo_boot = np.arccos(np.clip(np.dot(np.array([0.9162, 0.3068, 0.2577]), x3), -1.0, 1.0)) #el clip hace que si por algun motivo el dot me da >1 (i.e. 1,00002), me lo convierte en 1
print('El ángulo entre la normal del bootstrap y la del MVA = {0:1.3g}º'.format(angulo_boot * 180/np.pi))


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
print('La velocidad minima es = {0:1.3g} km/s, la diferencia entre la minima y la maxima es = {1:1.3g} km/s'.format(min(norma_v), diff))
#la velocidad promedio
v_media = np.array([np.mean(v_punto[:,0]), np.mean(v_punto[:,1]), np.mean(v_punto[:,2])])
print('La velocidad media de la nave a través de la MPB es = {0:1.3g} km/s'.format(np.linalg.norm(v_media)))

#ahora quiero ver si la nave atraviesa perpendicularmente a la MPB
v_media_norm = v_media/np.linalg.norm(v_media) #la normalizo por cuestion de  cuentas nomas
angulo_v_ajuste = np.arccos(np.clip(np.dot(norm_vignes, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_v_mva = np.arccos(np.clip(np.dot(x3, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre la velocidad y la normal del ajuste es = {0:1.3g}º'.format(angulo_v_ajuste * 180/np.pi))
print('El ángulo entre la velocidad y la normal del MVA es = {0:1.3g}º'.format(angulo_v_mva * 180/np.pi))

B_intermedio = B_medio_vectorial/np.linalg.norm(B_medio_vectorial) #la normalizo por cuestion de  cuentas nomas
angulo_B_ajuste = np.arccos(np.clip(np.dot(norm_vignes, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
angulo_B_mva = np.arccos(np.clip(np.dot(x3, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre el campo B y la normal del ajuste es = {0:1.3g}º'.format(angulo_B_ajuste * 180/np.pi))
print('El ángulo entre el campo B y la normal del MVA es = {0:1.3g}º'.format(angulo_B_mva * 180/np.pi))



#ahora veamos v_para
v_para = np.dot(v_media, norm_vignes) * norm_vignes

#deltat = np.loadtxt('tiempos.txt', delimiter=',', usecols = (3,4,5,6))
# deltat = deltat[-1] #si estoy de acuerdo con usar la última seleccion de puntos que hice

deltat_14 = (t4 - t1) * 3600
deltat_23 = (t3 - t2) * 3600

x_14 = v_para * deltat_14 #en km# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x_23 = v_para * deltat_23

print('El ancho de la MPB proyectando sobre la normal de Vignes para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14), np.linalg.norm(x_23)))

#si ahora proyecto sobre la normal de la MVA
v_para_MVA = np.dot(v_media, x3) * x3

x_14_MVA = v_para_MVA * deltat_14 #en km
x_23_MVA = v_para_MVA * deltat_23

print('El ancho de la MPB proyectando sobre la normal del MVA para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14_MVA), np.linalg.norm(x_23_MVA)))

for i in [t_1,t_2,t_3,t_4]:
    sza = np.arccos(np.clip(np.dot(posicion[i]/np.linalg.norm(posicion[i]), [1,0,0]), -1.0, 1.0))
    altitud = MD[i, 8]
    print('El tiempo = {2:1.6g} tiene SZA = {0:1.3g}º y altitud = {1:1.3g} km '.format(sza * 180/np.pi, altitud, t[i]))

fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1, projection='3d')
ax1.set_xlabel('x mso')
ax1.set_ylabel('y mso')
ax1.set_zlabel('z mso')
ax1.set_aspect('equal')
plot = ax1.plot_surface(
    X1, Y1, Z1, rstride=4, cstride=4, cmap=plt.get_cmap('Blues_r'), alpha=0.5)
# ax1.scatter(R[0], R[1], R[2])
ax1.plot(orbita[6000:11000,0], orbita[6000:11000,1], orbita[6000:11000,2])
# ax1.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)
ax1.quiver(R[0], R[1], R[2], norm_vignes[0], norm_vignes[1], norm_vignes[2], color='g',length=0.5, label='fit normal') #asi se plotea un vector
ax1.quiver(R[0], R[1], R[2], v_media[0], v_media[1], v_media[2], color='b',length=0.5, label='velocity') #asi se plotea un vector
ax1.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color='k',length=0.5, label='MVA normal')
ax1.quiver(R[0], R[1], R[2], v_para[0], v_para[1], v_para[2], color='r',length=0.5, label='v parallel') #asi se plotea un vector
ax1.quiver(R[0], R[1], R[2], v_para_MVA[0], v_para_MVA[1], v_para_MVA[2], color='m',length=0.5, label='v parallel MVA') #asi se plotea un vector
ax1.legend()
set_axes_equal(ax1)

#########
inicio_up = np.where(t == find_nearest_inicial(t, t1-0.015))[0][0] #las 18:12:00
fin_up = np.where(t == find_nearest_final(t, t1))[0][0] #las 18:13:00
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0) #nT

inicio_down = np.where(t == find_nearest_inicial(t, t4))[0][0] #las 18:14:51
fin_down = np.where(t == find_nearest_final(t, t4+0.015))[0][0] #las 18:15:52
B_downstream = np.mean(B[inicio_down:fin_down,:], axis=0) #nT

print('B upstream es {} nT y su módulo es {}'.format(B_upstream, np.linalg.norm(B_upstream)))
print('B downstream es {} nT y su módulo es {}'.format(B_downstream, np.linalg.norm(B_downstream)))

omega = np.arccos(np.dot(B_upstream,B_downstream)/(np.linalg.norm(B_upstream)*np.linalg.norm(B_downstream)))
print('El ángulo omega es {0:1.3g}º'.format(omega * 180/np.pi))

mu = 4* np.pi * 1E-7 #Henry/m

J_s = np.cross(x3, (B_upstream-B_downstream)) / mu #nA/m

ancho_mpb = np.linalg.norm(x_23_MVA) #considero que es el tiempo corto a lo largo de la normal del MVA
J_v = J_s / (1000*ancho_mpb) #nA/m²

print('La corriente superficial con la normal del MVA es Js = {} mA/m, |Js| = {} mA/m'.format(J_s *1E-6, np.linalg.norm(J_s)*1E-6))
print('La corriente en volumen con la normal del MVA es Jv = {} nA/m², |Jv| = {} nA/m²'.format(J_v, np.linalg.norm(J_v)))

J_s_ajuste = np.cross(norm_vignes, (B_upstream-B_downstream)) / mu #nA/m

ancho_mpb_ajuste = np.linalg.norm(x_23) #considero que es el tiempo corto a lo largo de la normal del MVA
J_v_ajuste = J_s_ajuste / (1000*ancho_mpb) #nA/m²
print('La corriente superficial con la normal del ajuste es Js = {} mA/m, |Js| = {} mA/m'.format(J_s_ajuste *1E-6, np.linalg.norm(J_s_ajuste)*1E-6))
print('La corriente en volumen con la normal del ajuste es Jv = {} nA/m², |Jv| = {} nA/m²'.format(J_v_ajuste, np.linalg.norm(J_v_ajuste)))


# fuerza_mva = np.cross(J_v * 1E-9, B[inicio_up:fin_down, :] * 1E-9) #A/m² * T = N/m³
# fuerza_ajuste = np.cross(J_v_ajuste * 1E-9, B[inicio_up:fin_down, :] * 1E-9) #A/m² * T = N/m³
fuerza_mva = np.cross(J_v * 1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
fuerza_ajuste = np.cross(J_v_ajuste * 1E-9, B[inicio_down,:]*1E-9) #N/m^3 #en t4
print('La fuerza de lorentz del MVA es {} V/m, su magnitud es {}'.format(fuerza_mva, np.linalg.norm(fuerza_mva)))
print('La fuerza de lorentz del ajuste es {} V/m, su magnitud es {}'.format(fuerza_ajuste, np.linalg.norm(fuerza_ajuste)))


e_density = lpw[:,3]
ti_lpw = np.where(t_lpw == find_nearest(t_lpw, ti))[0][0]
tf_lpw = np.where(t_lpw == find_nearest(t_lpw, tf))[0][0]
n_e = np.nanmean(e_density[ti_lpw:tf_lpw]) #hace el mean ignorando los nans #cm⁻³
n_e = n_e * 1E6 #m⁻³
# n_e = 1E7
q_e = 1.6E-19 #carga electron #C

E_Hall = np.cross(J_v * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
print('El campo de Hall del MVA es {} mV/m, su magnitud es {}'.format(E_Hall*1E3, np.linalg.norm(E_Hall)*1E3))

E_Hall_ajuste = np.cross(J_v_ajuste * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
print('El campo de Hall del ajuste es {} mV/m, su magnitud es {}'.format(E_Hall_ajuste*1E3, np.linalg.norm(E_Hall_ajuste)*1E3))


fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1, projection='3d')
ax2.set_xlabel(r'$X_{MSO} (R_M)$')
ax2.set_xlim(left=2, right=0)
ax2.set_ylabel(r'$Y_{MSO} (R_M)$')
ax2.set_zlabel(r'$Z_{MSO} (R_M)$')
ax2.set_aspect('equal')
ax2.scatter(R[0], R[1], R[2])
plot = ax2.plot_surface(
    X1, Y1, Z1, rstride=4, cstride=4, alpha=0.5, cmap=plt.get_cmap('Blues_r'))
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
ax2.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)
ax2.quiver(R[0], R[1], R[2], J_v[0], J_v[1], J_v[2], color='r',length=1E-3, label='Corriente en volumen')
ax2.quiver(R[0], R[1], R[2], B_upstream[0], B_upstream[1], B_upstream[2], color='b',length=1E-2, label='B upstream')
ax2.quiver(R[0], R[1], R[2], B_downstream[0], B_downstream[1], B_downstream[2], color='g',length=1E-2, label='B downstream')
ax2.quiver(R[0], R[1], R[2], fuerza_mva[0], fuerza_mva[1],fuerza_mva[2], color='m',length=2E13, label='Fuerza MVA')
ax2.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color='k',length=0.5, label='Normal del MVA', linewidths=0.5)
set_axes_equal(ax2)
ax2.legend(loc='upper right',bbox_to_anchor=(1.1, 1.05))

# plt.figure()
# plt.plot(t[inicio_up:fin_down], fuerza_mva)
# plt.plot(t[inicio_up:fin_down], fuerza_ajuste)
# plt.legend(['Fx', 'Fy', 'Fz'])
# plt.xlabel('Tiempo')
# plt.ylabel('Fuerza (N/m^3)')
plt.show(block=False) #para que siga andando aunque no lo haya cerrado

# #bootstrap error
# out = np.zeros(1000)
# out_phi = np.zeros((1000, 2))
# normal_ran = np.zeros((1000, 3))
# for a in range(1000):
#     index = np.random.choice(B_cut.shape[0], M_cut, replace=True) #elije M índices de B, puede repetir (replace = True)
#     B_random = B_cut[index,:] #me da un B hecho a partir de estos índices random
#
#     Mij_random = Mij(B_random)
#
#     [lamb_ran, x_ran] = np.linalg.eigh(Mij_random) #uso eigh porque es simetrica
#     idx_ran = lamb_ran.argsort()[::-1]
#     lamb_ran = lamb_ran[idx_ran]
#     x_ran = x_ran[:,idx_ran]
#     #ojo que a veces me da las columnas en vez de las filas como autovectores: el av x1 = x[:,0]
#     x3_ran = x_ran[:,2]
#     if x3_ran[0] < 0:
#         x3_ran = - x3_ran
#     normal_ran[a, :] = x3_ran
#
#     B3_ran = np.dot(B_random, x3_ran)
#     phi, delta_B3 = error(lamb_ran, B_random, M_cut, x_ran, dia)
#     out[a] = np.mean(B3_ran)
#     out_phi[a,0] = phi[2,0]
#     out_phi[a,1] = phi[2,1]
#
# normal_boot = np.linalg.norm(normal_ran, axis=0)/np.linalg.norm(normal_ran)
#
# plt.figure()
# plt.subplot(311)
# n, bins, patches = plt.hist(out, 50, normed = True, alpha=0.5)
# (muB, sigmaB) = norm.fit(out)
# y = normpdf(bins, muB, sigmaB)
# plt.plot(bins, y)
# plt.xlabel(r'$\langle B_3 \rangle$ (nT)')
#
# plt.subplot(312)
# n, bins, patches = plt.hist(out_phi[:,0]*57.2958, 50, normed=1, alpha=0.5)
# (mu31, sigma31) = norm.fit(out_phi[:,0]*57.2958)
# y = normpdf(bins, mu31, sigma31)
# plt.plot(bins, y)
# plt.xlabel(r'$\Delta \phi_{31}$ (º)')
#
# plt.subplot(313)
# n, bins, patches = plt.hist(out_phi[:,1]*57.2958, 50, normed=1, alpha=0.5)
# (mu32, sigma32) = norm.fit(out_phi[:,1]*57.2958)
# y = normpdf(bins, mu32, sigma32)
# plt.plot(bins, y)
# plt.xlabel(r'$\Delta \phi_{32}$ (º)')
# plt.tight_layout()
#
#
# plt.show(block=False) #para que siga andando aunque no lo haya cerrado
#
# print('mean_B = {0:1.3g} nT, std_B={1:1.3g} nT'.format(muB, sigmaB))
# print('mean_phi31 = {0:1.3g}º, std_phi31={1:1.3g}º'.format(mu31, sigma31))
# print('mean_phi32 = {0:1.3g}º, std_phi32={1:1.3g}º'.format(mu32, sigma32))
# print('la normal del bootstrap es {}'.format(normal_boot))
# angulo_v_boot = np.arccos(np.clip(np.dot(normal_boot, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
# print('El ángulo entre la velocidad y la normal del bootstrap es = {0:1.3g}º'.format(angulo_v_boot * 180/np.pi))
#
# angulo_B_boot = np.arccos(np.clip(np.dot(normal_boot, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
# print('El ángulo entre el campo B y la normal del bootstrap es = {0:1.3g}º'.format(angulo_B_boot * 180/np.pi))
#
# angulo_normales = np.arccos(np.clip(np.dot(normal_boot, x3), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
# print('El ángulo entre la normal del mva y la normal del bootstrap es = {0:1.3g}º'.format(angulo_normales * 180/np.pi))
#
# #si ahora proyecto sobre la normal de la MVA
# v_para_boot = np.dot(np.array([-2.487,  0.479,  2.836]), normal_boot) * normal_boot
#
# x_14_boot = v_para_boot * deltat_14 #en km
# x_23_boot = v_para_boot * deltat_23
#
# print('El ancho de la MPB proyectando sobre la normal del bootstrap para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14_boot), np.linalg.norm(x_23_boot)))
#
# J_s_boot = np.cross(normal_boot, (B_upstream-B_downstream)) / mu #nA/m
#
# ancho_mpb_boot = 77 #km, esta copiado de lo que da en el otro script
# J_v_boot = J_s_boot / (1000*ancho_mpb_boot) #nA/m²
# print('La corriente superficial con la normal del bootstrap es Js = {} mA/m, |Js| = {} mA/m'.format(J_s_boot *1E-6, np.linalg.norm(J_s_boot)*1E-6))
# print('La corriente en volumen con la normal del bootstrap es Jv = {} nA/m², |Jv| = {} nA/m²'.format(J_v_boot, np.linalg.norm(J_v_boot)))
#
# fuerza_mva_boot = np.cross(J_v_boot * 1E-9, B[inicio_down,:]*1E-9) #N/m^3
# print('La fuerza de lorentz es {} V/m, su magnitud es {}'.format(fuerza_mva_boot, np.linalg.norm(fuerza_mva_boot)))
# E_Hall_boot = np.cross(J_v_boot * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
# print('El campo de Hall es {} mV/m, su magnitud es {}'.format(E_Hall_boot*1E3, np.linalg.norm(E_Hall_boot)*1E3))
