import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""Time [UNIX],Time [UTC],Density H+ [/cc],Density He++ (m/q=2) [/cc],Density O+ [/cc],
Density O2+ [/cc],Density CO2+ [/cc],Uncertainty in H+ [%],Uncertainty in He++ [%],
Uncertainty in O+ [%],Uncertainty in O2+ [%],Uncertainty in CO2+ [%],SZA [degrees],
"Altitude [km, IAU]","Pos X [MSO, km]","Pos Y [MSO, km]","Pos Z [MSO, km]",Radius pos [km]"""

path = "../../datos/STATIC/"
datos = pd.read_csv(path + "STATIC_Ni_2016-03-16.txt")

t_hdec = []
for t_UTC in datos["Time [UTC]"]:
    (dia, hora) = t_UTC.split("/")
    (h, m, s) = hora.split(":")
    t_hdec.append(int(h) + int(m) / 60 + int(s) / 3600)

plt.semilogy(t_hdec, datos["Density H+ [/cc]"], ".", label="H+")
plt.semilogy(t_hdec, datos["Density O+ [/cc]"], ".", label="O+")
plt.semilogy(t_hdec, datos["Density O2+ [/cc]"], ".", label="O2+")
plt.semilogy(t_hdec, datos["Density CO2+ [/cc]"], ".", label="CO2+")
plt.ylim(ymin=0.1)
plt.legend()
plt.show()
