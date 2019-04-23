import numpy as np
from datetime import datetime, date

"""
Este código va a buscar entre todos los archivos que tenemos de MAG los lapsos en los cuales se cumplen:
SZA < 45º (o incluso <30º)
Altitud entre 300 y 1300 km
Z_MSO > 0 (si después se pueden volver a bajar los datos y hacer Z_pc > 0, mejor)

A su vez, va a guardar cada una de estas condiciones y ver cuántas veces por mes se repite

tenemos datos desde 10/2014 hasta 02/2018
"""
#
# date_entry = input('Enter a date in YYYY-MM-DD format \n')
# year, month, day = map(int, date_entry.split('-'))
# date_orbit = date(year, month, day)
#
# year = date_orbit.strftime("%Y")
# month = date_orbit.strftime("%m")
# day = date_orbit.strftime("%d")
# doty = date_orbit.strftime("%j")

# path = '~/../MAVEN/mag_1s/{0}/{1}/'.format(year, month)
# mag = np.loadtxt(path + 'mvn_mag_l2_{0}{3}ss1s_{0}{1}{2}_v01_r01.sts'.format(year, month, day, doty), skiprows=148)

path = '~/../MAVEN/mag_1s/2016/03/'
mag = np.loadtxt(path + 'mvn_mag_l2_2016091ss1s_20160331_v01_r01.sts', skiprows=148)
