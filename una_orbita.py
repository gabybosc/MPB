import numpy as np
from os import listdir
import glob as glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import scipy.signal as signal
from funciones import find_nearest, set_axes_equal
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos en los cuales se cumplen:
SZA < 45º (o incluso <30º)
Altitud entre 300 y 1300 km
Z_MSO > 0 (si después se pueden volver a bajar los datos y hacer Z_pc > 0, mejor)

Se fija dónde es que coincide la posicion de MAVEN con el fit de vignes y mira estas condiciones en ese punto.

tenemos datos desde 10/2014 hasta 02/2018
"""

# Ajuste de Vignes:
x0 = 0.78
e = 0.9
L = 0.96

theta = np.linspace(0, 3*np.pi/4, 1000)
phi = np.linspace(0, np.pi, 1000)
THETA, PHI = np.meshgrid(theta, phi)

r = L / (1 + e * np.cos(theta))

X = x0 + r * np.cos(theta)
Y = r * np.sin(theta) * np.cos(phi)
Z = r * np.sin(theta) * np.sin(phi)

R = np.transpose(np.array([X,Y,Z]))
# mag = np.loadtxt('../../../MAVEN/mag_1s/2016/03/mvn_mag_l2_2016085ss1s_20160325_v01_r01.sts', skiprows=148) #en la compu del iafe
mag = np.loadtxt('../../datos/MAG_1s/mvn_mag_l2_2016092ss1s_20160401_v01_r01.sts', skiprows=148) #en mi compu

B = np.zeros((len(mag[:,0]), 3))
for j in range(7,10):
    B[:,j-7] = mag[:, j]

B_norm = np.linalg.norm(B, axis=1)

posicion = np.zeros((len(mag[:,0]), 3))
for j in range(11,14):
    posicion[:,j-11] = mag[:, j]

orbita = posicion / 3390 #radios marcianos
#quiero encontrar cuándo la órbita cruza a la superficie dada por (X,Y,Z)

idx = np.zeros(len(orbita))
for i in range(len(orbita)):
    idx[i] = (np.abs(orbita[i,:]-R)).argmin()

# b,a = signal.butter(3,0.01,btype='lowpass')
# filtered = signal.filtfilt(b, a, B_norm)
# peaks = signal.find_peaks(filtered, 40)
# min = peaks[0][0]-2000
# max = peaks[0][0]+2000
# rango = posicion[min:max]
#
# find_nearest(rango, np.transpose(R))

# fig = plt.figure()
# ax = fig.add_subplot(1,1,1, projection='3d')
# ax.set_xlabel(r'$X_{MSO} (R_m)$')
# ax.set_ylabel(r'$Y_{MSO} (R_m)$')
# ax.set_zlabel(r'$Z_{MSO} (R_m)$')
# ax.set_aspect('equal')
# # ax.plot(orbita[:,0], orbita[:,1], orbita[:,2], color='C2', label='Órbita')
# ax.plot(orbita[:,0], orbita[:,1], orbita[:,2], color='C2', label='Órbita')
# plot = ax.plot_surface(
#     X, Y, Z, rstride=4, cstride=4, alpha=0.5, edgecolor='none', cmap=plt.get_cmap('Blues_r'))
# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="r", linewidth=0.5)
# ax.legend()
# set_axes_equal(ax) #para que tenga forma de esfera la esfera
# plt.show(block=False)


"""
Clasificación por SZA, es el que menos varía. Si el SZA medio es < 45, probablemente todos los SZA lo sean.
"""
# B_norm = np.linalg.norm(B, axis=1)
#
# b,a = signal.butter(3,0.01,btype='lowpass')
# filtered = signal.filtfilt(b, a, B_norm)
# peaks = signal.find_peaks(filtered, 40)
# mpb = peaks[0]-500
#
# SZA = np.zeros(len(mpb))
#
# for j in range(len(mpb)):
#     SZA[j] = np.arccos(np.clip(np.dot(posicion[mpb[j]]/np.linalg.norm(posicion[mpb[j]]), [1,0,0]), -1.0, 1.0))* 180/np.pi
#
# altitud = np.linalg.norm(posicion[mpb], axis=1) - 3390
#
# plt.figure(0)
# plt.plot(filtered)
#
# plt.figure(1)
# plt.plot(SZA)
# plt.show()
