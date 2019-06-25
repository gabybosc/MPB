import urllib.request
import shutil
import numpy as np
import datetime as dt


"""
descarga los archivos de mag hi res (pues si llegue hasta acá es porque ya vi los low res), lpw, swea y swia de la lista de fechas del año que le pida
El problema es que no puede asegurarse de que existan los archivos antes porque la página no tira 404. (Si tira 404, acá está la solución https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head)
"""

fechas = np.loadtxt('outputs/fechas_MVA_2016.txt')
for j in range(len(fechas)):
    year = 2016
    doy = int(fechas[j])


    date_orbit = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1) #para convertir el doty en date

    year = date_orbit.strftime("%Y")
    month = date_orbit.strftime("%m")
    day = date_orbit.strftime("%d")
    doy = date_orbit.strftime("%j")

    mag_hires = f'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/highres/{year}/{month}/mvn_mag_l2_{year}{doy}ss_{year}{month}{day}_v01_r01.sts'

    swea = f'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swea.calibrated/data/svy_spec/{year}/{month}/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf'

    lpw = f'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.lpw.derived/data/lp-nt/{year}/{month}/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf'

    swia = f'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.swia.calibrated/data/onboard_svy_mom/{year}/{month}/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf'

    # with urllib.request.urlopen(mag_hires) as response, open(f'../../datos/MAG_hires/mvn_mag_l2_{year}{doy}ss1s_{year}{month}{day}_v01_r01.sts', 'wb') as out_file:
    #     shutil.copyfileobj(response, out_file)
    # print(f'mag dia {doy} listo')

    with urllib.request.urlopen(swea) as response, open(f'../../datos/SWEA/mvn_swe_l2_svyspec_{year}{month}{day}_v04_r01.cdf', 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    print(f'swea dia {doy} listo')

    # with urllib.request.urlopen(swia) as response, open(f'../../datos/SWIA/mvn_swi_l2_onboardsvymom_{year}{month}{day}_v01_r01.cdf', 'wb') as out_file:
    #     shutil.copyfileobj(response, out_file)
    # print(f'swia dia {doy} listo')
    #
    # with urllib.request.urlopen(lpw) as response, open(f'../../datos/LPW/mvn_lpw_l2_lpnt_{year}{month}{day}_v03_r02.cdf', 'wb') as out_file:
    #     shutil.copyfileobj(response, out_file)
    # print(f'lpw dia {doy} listo')

    print(f'voy {int(j/len(fechas)*100)}%')
