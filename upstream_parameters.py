import numpy as np

path = "../../datos/"
n_ion = np.loadtxt(path + "n_ion.asc")[:, -1]
T_ion = np.loadtxt(path + "t_ion.asc")[:, -1]
v_ion = np.loadtxt(path + "v_ion.asc")[:, -3:]
T_e = np.loadtxt(path + "T_e.asc")[:, -1]
B = np.loadtxt(path + "B.asc")[:, -3:]

v_ion_norm = np.linalg.norm(v_ion, axis=1)
B_norm = np.linalg.norm(B, axis=1)

print("n ion = ", np.mean(n_ion))
print("T ion = ", np.mean(T_ion))
print("v ion = ", np.mean(v_ion_norm))
print("T electron = ", np.mean(T_e))
print("B mean = ", np.mean(B_norm))
