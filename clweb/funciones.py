import numpy as np
import matplotlib.dates as md
import pandas as pd
from numpy.random import rand
import datetime as dt
import calendar

def deltaB(B):
    B_medio = np.mean(B, axis=0)
    abs_deltaB_para = np.abs(np.dot(B - B_medio, B_medio)) / np.linalg.norm(B_medio)**2   #|deltaB_para / B|

    dot = np.dot(B - B_medio, B_medio / np.linalg.norm(B_medio))
    N = np.zeros((len(dot), len(B_medio)))

    for i in range(len(N)):
        N[i,:] = dot[i] * B_medio/np.linalg.norm(B_medio)
        deltaB_perp = (B - B_medio) - N
        #y ahora necesito el valor abs de perp
        abs_deltaB_perp = np.abs(deltaB_perp) / np.linalg.norm(B_medio)

    return(abs_deltaB_para, abs_deltaB_perp)

def Bpara_Bperp(B, t, ti, tf):
    j_inicial = int(np.where(t == find_nearest(t, ti))[0][0])
    j_final =  int(np.where(t == find_nearest(t, tf))[0][0])

    #Lo hago en ventanas de 60s, moviendose de a 1s.
    B_para = np.zeros(j_final - j_inicial)
    B_perp = np.zeros((j_final - j_inicial, 3))
    B_perp_norm = np.zeros(j_final - j_inicial)
    for j in range(j_inicial, j_final):
        Mi = j
        Mf = j + 25
        M_delta = 25
        B_delta = B[Mi:Mf]
        t_delta = t[Mi:Mf]
        deltaB_para, deltaB_perp = deltaB(B_delta)
        B_para[j-j_inicial] = deltaB_para[12]
        B_perp[j-j_inicial, :] = deltaB_perp[12, :]
        B_perp_norm[j-j_inicial] = np.linalg.norm(deltaB_perp[12,:])


    return(B_para, B_perp_norm, j_inicial, j_final)

def error(lamb, B, M, x):
    phi = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i == j:
                phi[i,j] = 0
            else:
                phi[i,j] = np.sqrt(lamb[2] / (M - 1) * (lamb[i] + lamb[j] - lamb[2])/(lamb[i] - lamb[j])**2) #en radianes


    delta_B3 = np.sqrt(lamb[2] / (M-1) + (phi[2,1] * np.dot(np.mean(B,0), x[1]) )**2 + (phi[2,0] * np.dot(np.mean(B,0), x[0]))**2)
    return(phi, delta_B3)


#para encontrar el valor mas cercano
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin() #argmin me da el valor del minimo en el array
    return array[idx]

def find_nearest_inicial(array,value): #el valor más cercano pero más grande
    idx = (np.abs(array-value)).argmin() #el indice de la minima dif entre el array y el valor
    if array[idx] < value: #si es menor que el valor que busco:
        idx = np.abs((array - (value + 1/3600))).argmin() #el indice va a ser la min dif entre el array y el valor+1seg
    return array[idx]

def find_nearest_final(array,value): #el valor más cercano pero más chico
    idx = (np.abs(array-value)).argmin() #el indice de la minima dif entre el array y el valor
    if array[idx] > value: #si es mayor que el valor que busco:
        idx = np.abs((array - (value - 1/3600))).argmin() #el indice va a ser la min dif entre el array y el valor-1seg
    return array[idx]


def Mij(B):
    Mij = np.zeros((3,3))
    for i in range(3): #para las tres coordenadas
        for j in range(3):
            Mij[i,j] = np.mean(B[:,i] * B[:,j]) - np.mean(B[:,i]) * np.mean(B[:,j])
    return Mij

def next_available_row(sheet):
    str_list = list(filter(None, sheet.col_values(1)))
    return str(len(str_list)+1)

def rms(x):
    rms = np.sqrt(np.vdot(x, x)/x.size)
    return rms

def unix_to_decimal(t_unix):
    t = np.zeros(np.size(t_unix))
    for i in range(np.size(t)):
        u = dt.datetime.utcfromtimestamp(int(t_unix[i]))
        t[i] = u.second/3600 + u.minute/60 + u.hour #tiempo decimal
    return(t)
    #t me da el mismo array que tengo en el MVA!

def unix_to_timestamp(t_unix):
    n = len(t_unix)
    timestamps = np.linspace(t_unix[0], t_unix[-1], n)
    dates = [dt.datetime.utcfromtimestamp(ts) for ts in timestamps] #me lo da en UTC
    datenums = md.date2num(dates)
    return(datenums)

def UTC_to_hdec(t_UTC):
    (h, m, s) = t_UTC.split(':')
    t_hdec = int(h) + int(m) / 60 + int(s) / 3600

    return(t_hdec)

def getrem(input):
    "this function yields the value behind the decimal point"
    import numpy as np
    output=abs(input-np.fix(input))
    return output

def datenum(Yr,Mo=1,Da=1,Hr=0,Mi=0,Se=0,Ms=0):
    "this function works as regular datetime.datetime, but allows for float input"

    #correct faulty zero input
    if Mo<1:
        Mo+=1
    if Da<1:
        Da+=1

    #distribute the year fraction over days
    if  getrem(Yr)>0:
        if calendar.isleap(np.floor(Yr)):
            fac=366
        else:
            fac=365
        Da=Da+getrem(Yr)*fac
        Yr=int(Yr)
    #if months exceeds 12, pump to years
    while int(Mo)>12:
        Yr=Yr+1
        Mo=Mo-12
    #distribute fractional months to days
    if getrem(Mo)>0:
        Da=Da+getrem(Mo)*calendar.monthrange(Yr,int(Mo))[1]
        Mo=int(Mo)
    #datetime input for 28 days always works excess is pumped to timedelta
    if Da>28:
        extraDa=Da-28
        Da=28
    else:
        extraDa=0
    # sometimes input is such that you get 0 day or month values, this fixes this anomaly
    if int(Da)==0:
       Da+=1
    if int(Mo)==0:
       Mo+=1

    #datetime calculation
    mytime=dt.datetime(int(Yr),int(Mo),int(Da))+dt.timedelta(days=extraDa+getrem(Da),hours=Hr,minutes=Mi,seconds=Se,microseconds=Ms)
    return mytime

def array_datenums(year, month, day, t):
    year = int(year)
    month = int(month)
    day = int(day)
    timestamps = np.array([np.datetime64(datenum(year, month, day, x)) for x in t]) #datenum es una función mía
    return(timestamps)
