import numpy as np
import csv as csv
import pandas as pd
from funciones import donde

catalogo = np.genfromtxt("outputs/catalogo_grupo2.txt", dtype="str")
fecha = np.genfromtxt("outputs/catalogo_grupo2.txt", dtype="str", usecols=0)
hora = np.genfromtxt("outputs/catalogo_grupo2.txt", dtype="str", usecols=1)
theta = np.genfromtxt("outputs/catalogo_grupo2.txt", dtype="float", usecols=2)
pdyn = np.genfromtxt("outputs/catalogo_grupo2.txt", dtype="float", usecols=3)

"""
El catálogo es muy grande, entonces hace una subselección
"""

idx_para = []
for i in range(len(fecha)):
    ang = theta[i]
    if ang < 20:
        idx_para.append(i)

fecha_para = fecha[idx_para]
hora_para = hora[idx_para]
theta_para = theta[idx_para]
p_para = pdyn[idx_para]

"""
estos fueron los i donde es quasipara, así que están preseleccionados y
los saco del filtro en idx_perp
"""

"""
No quiero pdyn repetidos
"""

p_unique, idx = np.unique(pdyn, return_index=True)
# ojo que esto me los da ordenados de menor a mayor

idx_perp = [i for i in idx if i not in idx_para]

fecha_cut = fecha[idx_perp]
hora_cut = hora[idx_perp]
theta_cut = theta[idx_perp]
p_cut = pdyn[idx_perp]

"""
Sin theta repetidos
"""

th_unique, idx_th = np.unique(theta_cut, return_index=True)

fecha_final = fecha_cut[idx_th]
hora_final = hora_cut[idx_th]
theta_final = theta_cut[idx_th]
p_final = p_cut[idx_th]

"""
Los escribo en un archivo
No olvidarse de agregar los para!!
"""

# los perp
with open("outputs/hoja_grupo2.txt", "a") as file:
    for i in range(len(fecha_final)):
        file.write(f"{fecha_final[i]}\t{hora_final[i]}\t{theta_final[i]}\t{p_final[i]}")
        file.write("\n")

# los para
with open("outputs/hoja_grupo2.txt", "a") as file:
    for i in range(len(fecha_para)):
        file.write(f"{fecha_para[i]}\t{hora_para[i]}\t{theta_para[i]}\t{p_para[i]}")
        file.write("\n")
