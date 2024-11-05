import numpy as np
from generar_npys import generar_npys_limites
import pandas as pd

"""
lee el txt de sofi y lo pasa a npy
lee el catálogo de jacob y agarra los datos que no están en lo de sofi y lo pasa a npy
"""

path = "../../../../datos/bs_mpb/"
g1 = np.genfromtxt(path + "group1_aberrated_withVignesfit.txt", dtype="str")
g2 = np.genfromtxt(
    path + "group2_aberrated_withVignesfit.txt", dtype="str", skip_header=1
)
g3 = np.genfromtxt(
    path + "group3_aberrated_withVignesfit.txt", dtype="str", skip_header=1
)
g4 = np.genfromtxt(
    path + "group4_aberrated_withVignesfit.txt", dtype="str", skip_header=1
)

allgroups = np.concatenate((g1, g2, g3, g4))
# for c in range(len(allgroups.T)):
#     np.save(path + f"{allgroups[0, c]}", allgroups[1:, c])

generar_npys_limites("../../outputs/allgroups/")
rbs_1 = np.load("../../outputs/grupo1/Rsd_bs.npy")
rbs_2 = np.load("../../outputs/grupo2/Rsd_bs.npy")
rbs_3 = np.load("../../outputs/grupo3/Rsd_bs.npy")
rbs_4 = np.load("../../outputs/grupo4/Rsd_bs.npy")
rmpb_1 = np.load("../../outputs/grupo1/Rsd_mpb.npy")
rmpb_2 = np.load("../../outputs/grupo2/Rsd_mpb.npy")
rmpb_3 = np.load("../../outputs/grupo3/Rsd_mpb.npy")
rmpb_4 = np.load("../../outputs/grupo4/Rsd_mpb.npy")

rbs = np.hstack((rbs_1, rbs_2, rbs_3, rbs_4))
rmpb = np.hstack((rmpb_1, rmpb_2, rmpb_3, rmpb_4))

rbslow = np.hstack((rbs_3, rbs_4))
rmpblow = np.hstack((rmpb_3, rmpb_4))

rbsh = np.hstack((rbs_1, rbs_2))
rmpbh = np.hstack((rmpb_1, rmpb_2))

np.save("../../outputs/allgroups/Rsd_bs_high.npy", rbsh)
np.save("../../outputs/allgroups/Rsd_mpb_high.npy", rmpbh)

"""
Ahora jacob
"""

jacob = pd.read_csv("../../outputs/allgroups/Jacob_catalogue_grupos.csv", decimal=",")

date = [fecha[:10] for fecha in jacob["TIME"]]
time = [fecha[11:] for fecha in jacob["TIME"]]
x = jacob["XMSO [km]"]
y = jacob["YMSO [km]"]
z = jacob["ZMSO [km]"]
theta = jacob["thetaBN [deg]"]
beta = jacob["beta_p"]
np.save(path + "beta_p.npy", np.array(beta), allow_pickle=True)
np.save(path + "date.npy", date)
np.save(path + "time.npy", time)
np.save(path + "vup.npy", np.array(jacob["|Vupfine| [km/s]"]))
np.save(path + "wci.npy", np.array(jacob["Wci [rad/s]"]))
np.save(path + "rci.npy", np.array(jacob["rci [km]"]))
np.save(path + "pdyn.npy", np.array(jacob["Pdyn [nPa]"]))
np.save(path + "cwpe.npy", np.array(jacob["c/Wpe [km]"]))
np.save(path + "cwpi.npy", np.array(jacob["c/Wpi [km]"]))
np.save(path + "tau_ci.npy", np.array(jacob["tau_ci [s]"]))
np.save(path + "cone_angle.npy", np.array(jacob["cone angle [º]"]))
np.save(path + "malf_local.npy", np.array(jacob["Malfven_local"]))
np.save(path + "malf_global.npy", np.array(jacob["Malfven_global"]))
np.save(path + "Mfms.npy", np.array(jacob["Mfms"]))
np.save(path + "Ls.npy", np.array(jacob["Ls [deg]"]))
