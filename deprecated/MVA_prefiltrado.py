import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from scipy.stats import norm
from datetime import datetime
from funciones import hodograma, error, find_nearest, find_nearest_final, find_nearest_inicial, deltaB, Mij,datenum
from funciones_plot import hodograma, set_axes_equal

np.set_printoptions(precision=4)

path = 'datos/marzo 2016/16/' #path a los datos
# datos = np.loadtxt(path + 'mvn_mag_l2_2016076ss1s_20160316_v01_r01.sts', skiprows=148) #lee todo y me da todo
# path = 'datos/marzo_2016_hires/' #path a los datos
datos = np.loadtxt(path + 'mag_hires_mso_16032016_17_19_filt.txt', skiprows=1) #datos filtrados
# cdf_swia = cdf.CDF(path + 'mvn_swi_l2_onboardsvymom_20160316_v01_r01.cdf')
lpw = np.loadtxt(path + 'mvn_kp_insitu_20160316_v14_r03_orbita18h.csv') #son los datos entre las 18 y las 19h
t_lpw = lpw[:,0] + lpw[:,1]/60 + lpw[:,2]/3600

ti = 18.227
tf = 18.235
dia = 16
t = datos[:,0] #para que me de sobre la cantidad de horas

M = np.size(t) #el numero de datos

#el campo
B = np.zeros((M, 3))
for i in [1,2,3]:
    B[:,i-1] = datos[:, i]

MD = np.zeros((M, 5))
MD[:, 0] = t
for i in range(1,4):
    MD[:, i] = B[:,i-1]
MD[:,4] = datos[:,-1]#la norma de B

#si quiero elegir entre ciertas horas:
t1 = find_nearest_inicial(t, ti)
t2 = find_nearest_final(t, tf)
inicio = np.where(t == t1)[0][0]
fin = np.where(t == t2)[0][0]

t_1 = np.where(t == find_nearest(t,18.2167))[0][0]
t_2 = np.where(t == find_nearest(t,18.2204))[0][0]
t_3 = np.where(t == find_nearest(t,18.235))[0][0]
t_4 = np.where(t == find_nearest(t,18.2476))[0][0]


#################
#ahora empieza el MVA con los datos que elegí
MD_cut = MD[inicio : fin+1, :]
M_cut = np.size(MD_cut[:,0])
B_cut = B[inicio:fin+1, :]

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

#el error
phi, delta_B3 = error(lamb, B_cut, M_cut, x)
# print('Matriz de incerteza angular (grados): \n{}'.format(phi  *  57.2958))
print('<B3> = {0:1.3g} +- {1:1.3g} nT'.format(np.mean(B3),delta_B3))


if phi[2,1] > phi[2,0]:
    print('El error phi32 = {0:1.3g}º es mayor a phi31 = {1:1.3g}º'.format(phi[2,1]*57.2958, phi[2,0]*180/np.pi))
else:
    print('El error phi31 = {0:1.3g}º es mayor a phi32 = {1:1.3g}º'.format(phi[2,0]*57.2958, phi[2,1]*180/np.pi))
#quiero ver si el error más grande es phi31 o phi32

#el B medio
B_medio_vectorial = np.mean(B_cut, axis=0)
print('El valor medio de B = {0:1.3g} nT'.format(np.linalg.norm(B_medio_vectorial)))#np.mean(B_cut, 0))) para que me de los tres valores medios
print('|B3|/B_medio = ', np.abs(np.mean(B3)/np.linalg.norm(B_medio_vectorial)))


deltat_14 = (18.2476 - 18.2167) * 3600
deltat_23 = (18.2204 - 18.235) * 3600

v_para_MVA = np.dot(np.array([-2.487,  0.479,  2.836]), x3) * x3

x_14 = v_para_MVA * deltat_14 #en km# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x_23 = v_para_MVA * deltat_23

print('El ancho de la MPB proyectando sobre la normal para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14), np.linalg.norm(x_23)))


inicio_up = np.where(t == find_nearest_inicial(t, 18.20))[0][0] #las 18:12:00
fin_up = np.where(t == find_nearest_final(t, 18.2167))[0][0] #las 18:13:00
B_upstream = np.mean(B[inicio_up:fin_up, :], axis=0) #nT

inicio_down = np.where(t == find_nearest_inicial(t, 18.2476))[0][0] #las 18:14:51
fin_down = np.where(t == find_nearest_final(t, 18.2645))[0][0] #las 18:15:52
B_downstream = np.mean(B[inicio_down:fin_down,:], axis=0) #nT

print('B upstream es {} nT y su módulo es {}'.format(B_upstream, np.linalg.norm(B_upstream)))
print('B downstream es {} nT y su módulo es {}'.format(B_downstream, np.linalg.norm(B_downstream)))


omega = np.arccos(np.dot(B_upstream,B_downstream)/(np.linalg.norm(B_upstream)*np.linalg.norm(B_downstream)))
print('El ángulo omega es {0:1.3g}º'.format(omega * 180/np.pi))

B_intermedio = B_medio_vectorial/np.linalg.norm(B_medio_vectorial) #la normalizo por cuestion de  cuentas nomas
angulo_B_mva = np.arccos(np.clip(np.dot(x3, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre el campo B y la normal del MVA es = {0:1.3g}º'.format(angulo_B_mva * 180/np.pi))


v_media_norm = np.array([-0.6541,  0.126 ,  0.7459]) #la normalizo por cuestion de  cuentas nomas
angulo_v_mva = np.arccos(np.clip(np.dot(x3, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre la velocidad y la normal del MVA es = {0:1.3g}º'.format(angulo_v_mva * 180/np.pi))

hodograma(B1, B2, B3, 'nT', 'MAVEN MAG MVA 18:13:33 - 18:14:06 UTC')

mu = 4* np.pi * 1E-7 #Henry/m

J_s = np.cross(x3, (B_upstream-B_downstream)) / mu #nA/m

ancho_mpb = 90 #km, esta copiado de lo que da en el otro script
J_v = J_s / (1000*ancho_mpb) #nA/m²
print('La corriente superficial con la normal del MVA es Js = {} mA/m, |Js| = {} mA/m'.format(J_s *1E-6, np.linalg.norm(J_s)*1E-6))
print('La corriente en volumen con la normal del MVA es Jv = {} nA/m², |Jv| = {} nA/m²'.format(J_v, np.linalg.norm(J_v)))

e_density = lpw[:,3]
ti_lpw = np.where(t_lpw == find_nearest(t_lpw, ti))[0][0]
tf_lpw = np.where(t_lpw == find_nearest(t_lpw, tf))[0][0]
n_e = np.nanmean(e_density[ti_lpw:tf_lpw]) #hace el mean ignorando los nans #cm⁻³
n_e = n_e * 1E6 #m⁻³
# n_e = 1E7
q_e = 1.6E-19 #carga electron #C

fuerza_mva = np.cross(J_v * 1E-9, B[inicio_down,:]*1E-9) #N/m^3
print('La fuerza de lorentz es {} V/m, su magnitud es {}'.format(fuerza_mva, np.linalg.norm(fuerza_mva)))
E_Hall = np.cross(J_v * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
print('El campo de Hall es {} mV/m, su magnitud es {}'.format(E_Hall*1E3, np.linalg.norm(E_Hall)*1E3))


plt.show(block=False)


#bootstrap error
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
    phi, delta_B3 = error(lamb_ran, B_random, M_cut, x_ran, dia)
    out[a] = np.mean(B3_ran)
    out_phi[a,0] = phi[2,0]
    out_phi[a,1] = phi[2,1]

normal_boot = np.linalg.norm(normal_ran, axis=0)/np.linalg.norm(normal_ran)

plt.figure()
plt.subplot(311)
n, bins, patches = plt.hist(out, 50, normed = True, alpha=0.5)
(muB, sigmaB) = norm.fit(out)
y = normpdf(bins, muB, sigmaB)
plt.plot(bins, y)
plt.xlabel(r'$\langle B_3 \rangle$ (nT)')

plt.subplot(312)
n, bins, patches = plt.hist(out_phi[:,0]*57.2958, 50, normed=1, alpha=0.5)
(mu31, sigma31) = norm.fit(out_phi[:,0]*57.2958)
y = normpdf(bins, mu31, sigma31)
plt.plot(bins, y)
plt.xlabel(r'$\Delta \phi_{31}$ (º)')

plt.subplot(313)
n, bins, patches = plt.hist(out_phi[:,1]*57.2958, 50, normed=1, alpha=0.5)
(mu32, sigma32) = norm.fit(out_phi[:,1]*57.2958)
y = normpdf(bins, mu32, sigma32)
plt.plot(bins, y)
plt.xlabel(r'$\Delta \phi_{32}$ (º)')
plt.tight_layout()


plt.show(block=False) #para que siga andando aunque no lo haya cerrado

print('mean_B = {0:1.3g} nT, std_B={1:1.3g} nT'.format(muB, sigmaB))
print('mean_phi31 = {0:1.3g}º, std_phi31={1:1.3g}º'.format(mu31, sigma31))
print('mean_phi32 = {0:1.3g}º, std_phi32={1:1.3g}º'.format(mu32, sigma32))
print('la normal del bootstrap es {}'.format(normal_boot))
angulo_v_boot = np.arccos(np.clip(np.dot(normal_boot, v_media_norm), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre la velocidad y la normal del bootstrap es = {0:1.3g}º'.format(angulo_v_boot * 180/np.pi))

angulo_B_boot = np.arccos(np.clip(np.dot(normal_boot, B_intermedio), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre el campo B y la normal del bootstrap es = {0:1.3g}º'.format(angulo_B_boot * 180/np.pi))

angulo_normales = np.arccos(np.clip(np.dot(normal_boot, x3), -1.0, 1.0))  #Es importante que los vectoers estén normalizados!
print('El ángulo entre la normal del mva y la normal del bootstrap es = {0:1.3g}º'.format(angulo_normales * 180/np.pi))

# angulo_normales = np.arccos(np.clip(np.dot(normal_boot, norm_vignes), -1.0, 1.0))
#si ahora proyecto sobre la normal de la MVA
v_para_boot = np.dot(np.array([-2.487,  0.479,  2.836]), normal_boot) * normal_boot

x_14_boot = v_para_boot * deltat_14 #en km
x_23_boot = v_para_boot * deltat_23

print('El ancho de la MPB proyectando sobre la normal del bootstrap para el tiempo largo es = {0:1.3g} km, para el tiempo corto es = {1:1.3g} km'.format(np.linalg.norm(x_14_boot), np.linalg.norm(x_23_boot)))

J_s_boot = np.cross(normal_boot, (B_upstream-B_downstream)) / mu #nA/m

ancho_mpb_boot = 77 #km, esta copiado de lo que da en el otro script
J_v_boot = J_s_boot / (1000*ancho_mpb_boot) #nA/m²
print('La corriente superficial con la normal del bootstrap es Js = {} mA/m, |Js| = {} mA/m'.format(J_s_boot *1E-6, np.linalg.norm(J_s_boot)*1E-6))
print('La corriente en volumen con la normal del bootstrap es Jv = {} nA/m², |Jv| = {} nA/m²'.format(J_v_boot, np.linalg.norm(J_v_boot)))

fuerza_mva_boot = np.cross(J_v_boot * 1E-9, B[inicio_down,:]*1E-9) #N/m^3
print('La fuerza de lorentz de bootstrap es {} V/m, su magnitud es {}'.format(fuerza_mva_boot, np.linalg.norm(fuerza_mva_boot)))
E_Hall_boot = np.cross(J_v_boot * 1E-9, B[inicio_down, :] * 1E-9) / (q_e * n_e) #V/m
print('El campo de Hall del bootstrap es {} mV/m, su magnitud es {}'.format(E_Hall_boot*1E3, np.linalg.norm(E_Hall_boot)*1E3))


##########

###############
# orbita = posicion[np.where(t == find_nearest_inicial(t, 17.6))[0][0] : np.where(t == find_nearest_final(t, 19))[0][0], :] / 3390 #radios marcianos
#
# # usamos vignes et al:
# x0 = 0.78
# e = 0.9
# L = 0.96
#
# #ec conica
# theta = np.linspace(0, np.pi *3/4, 100)
# phi = np.linspace(0, 2 * np.pi, 100)
# THETA, PHI = np.meshgrid(theta, phi)
#
# r = L / (1 + e * np.cos(THETA))
#
# ####dibujamos el punto por el que pasa la nave:
# t_nave = find_nearest(t,18.2303) #las 18:13:49, es el tiempo en el meio de la hoja de corriente
# index = np.where(t == t_nave)[0][0]
# R = posicion[index,:] / 3390 #la posicion de la nave en RM
#
#
# ####### Calculo mi propia elipse que pase por el punto.
# r0 = R - np.array([x0,0,0])
# theta0 = np.arccos(r0[0] / np.linalg.norm(r0))
#
# L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
# r1 = L0 / (1 + e * np.cos(THETA))
#
#
# fig = plt.figure()
# ax = fig.add_subplot(1,1,1, projection='3d')
# ax.set_xlabel(r'$X_{MSO} (R_m)$')
# ax.set_ylabel(r'$Y_{MSO} (R_m)$')
# ax.set_zlabel(r'$Z_{MSO} (R_m)$')
# ax.set_aspect('equal')
# ax.plot(orbita[:,0], orbita[:,1], orbita[:,2], color='green')
# ax.scatter(R[0], R[1], R[2])
# X1 = x0 + r1 * np.cos(THETA)
# Y1 = r1 * np.sin(THETA) * np.cos(PHI)
# Z1 = r1 * np.sin(THETA) * np.sin(PHI)
# plot = ax.plot_surface(
#     X1, Y1, Z1, rstride=4, cstride=4, cmap=plt.get_cmap('Blues_r'), alpha=0.5, edgecolor='none')
# asc = L0 / (1 - e**2)  #semieje mayor
# bsc = np.sqrt(asc*L0)
# csc = e*asc - x0 #donde está centrada. Hay que ver el signo
#
# norm_vignes = np.array([(R[0]+csc)*2 / asc**2, R[1]*2 / (bsc)**2, R[2]*2 / (bsc)**2]) #la normal de vignes
# norm_vignes = norm_vignes/np.linalg.norm(norm_vignes) #normalizado
# ax.quiver(R[0], R[1], R[2], norm_vignes[0], norm_vignes[1], norm_vignes[2], color='b',length=0.5, label='Normal del ajuste') #asi se plotea un vector
# ax.quiver(R[0], R[1], R[2], x3[0], x3[1], x3[2], color='k',length=0.5, label='Normal del MVA')
# ax.quiver(R[0], R[1], R[2], 0.9162, 0.3068, 0.2577, color='m',length=0.5, label='Normal del bootstrap')
#
# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)
# ax.legend()
# set_axes_equal(ax) #para que tenga forma de esfera la esfera
#
# print('la normal del ajuste es {}'.format(norm_vignes))
# print('El valor medio de B a lo largo de la normal del ajuste es = {0:1.3g} nT'.format(np.mean(np.dot(B_cut, norm_vignes))))#np.mean(B_cut, 0))) para que me de los tres valores medios
# print('|B3_ajuste|/B_medio = ', np.abs(np.mean(np.dot(B_cut, norm_vignes))/np.mean(B_cut)))
