"""
descarga los archivos de mag del pds
El problema es que no puede asegurarse de que existan los archivos antes porque
la página no tira 404. Si tira 404, acá está la solución:
https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head
creo que sería mejor que baje el mes entero, si bien ocupa más espacio y tiempo,
ya que no necesita el nombre exacto del archivo.
"""

import urllib.request
import shutil
import numpy as np
from socket import gethostname


year = 2009
doy = np.random.choice(366, 60)  # elije 60 días aleatorios del año para descargar
# doy = range(1, 366)

if gethostname() == "magneto2":
    path = f"../../../../../media/gabybosc/datos/VEX/{year}/"
else:
    path = f"../../../datos/VEX/{year}/"

for dd in doy:
    mag = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/venus_express_cleaned_high_res_mag/data/{year}/fg128HzY{str(year)[-2:]}D{str(dd).zfill(3)}.tab"

    with urllib.request.urlopen(mag) as response, open(
        path + f"VEX_MAG_{year}{str(dd).zfill(3)}.tab",
        "wb",
    ) as out_file:
        shutil.copyfileobj(response, out_file)
    print(f"mag dia {dd} listo")
