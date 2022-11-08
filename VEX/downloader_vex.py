"""
descarga los archivos de mag del pds
El problema es que no puede asegurarse de que existan los archivos antes porque
la página no tira 404. Si tira 404, acá está la solución:
https://stackoverflow.com/questions/20387246/checking-file-exists-before-download-using-head
creo que sería mejor que baje el mes entero, si bien ocupa más espacio y tiempo,
ya que no necesita el nombre exacto del archivo.
A veces el problema es que el archivo se llama "_r02" o algo así y por eso no
sirve.
"""

import urllib.request
import shutil
import numpy as np
from socket import gethostname


year = 2010
doy = np.random.choice(366, 60)  # elije 60 días aleatorios del año para descargar


if gethostname() == "magneto2":
    path = f"../../../../../media/gabybosc/datos/VEX/{year}/"
else:
    path = f"../../../datos/VEX/{year}/"

for dd in doy:
    mag = f"https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/venus_express_cleaned_high_res_mag/data/{year}/fg128HzY{str(year)[-2:]}D{dd}.tab"

    with urllib.request.urlopen(mag) as response, open(
        path + f"VEX_MAG_{year}{dd}.tab", "wb",
    ) as out_file:
        shutil.copyfileobj(response, out_file)
    print(f"mag dia {dd} listo")
