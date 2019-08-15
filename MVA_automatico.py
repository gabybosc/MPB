import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
import datetime as dt
from funciones import find_nearest, find_nearest_final, find_nearest_inicial, Mij, datenum, unix_to_decimal
from funciones_plot import hodograma

"""
Para datos de Mag de baja resolución.
Agarra los datos de t1t2t3t4 que están en el txt y les hace mva. Es bueno para una primera aproximación.
Devuelve un hodograma y el cociente de lambdas.

"""


np.set_printoptions(precision=4)

# #si tengo la fecha en dia-mes-año
# date_entry = input('Enter a date in YYYY-MM-DD format \n')
# year, month, day = map(int, date_entry.split('-'))
# date_orbit = dt.date(year, month, day)

#si tengo la fecha en dia del año
fechas = np.loadtxt('t1t2t3t4.txt')
for j in range(len(fechas)):
    year = int(fechas[j,0])
    doy = int(fechas[j,1])
    t1 = fechas[j,2]
    t2 = fechas[j,3]
    t3 = fechas[j,4]
    t4 = fechas[j,5]

    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")
    print(f'Fecha {doy}-{year} entre t2 = {t2} y t3 = {t3}')

    # path = '../../../MAVEN/mag_1s/2016/03/' #path a los datos desde la desktop
    path = '../../datos/' #path a los datos desde la laptop
    mag = np.loadtxt(path + f'MAG_1s/subsolares/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', skiprows=148)
    lpw = cdf.CDF(path + f'LPW/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf')
    #para ver las varaibles del cdf:
    # lpw.cdf_info()

    dia = mag[:,1]
    t = mag[:,6]  #el dia decimal
    t = (t - dia) * 24 #para que me de sobre la cantidad de horas

    M = np.size(t) #el numero de datos

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
    #las proyecciones
    B1 = np.dot(B_cut, x1)
    B2 = np.dot(B_cut, x2)
    B3 = np.dot(B_cut, x3)

    hodograma(B1, B2, B3, 'nT', f'MAVEN MAG MVA fecha {doy}-{year}')

plt.show()
