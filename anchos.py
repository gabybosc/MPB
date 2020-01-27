import matplotlib.pyplot as plt

tiempos = ['24-11-17','10-10-15', '12-10-15', '05-04-16', '31-03-16', '16-03-16']
hmin = [115, 39, 18.5, 44.4, 38.5, 82.2]
hmax = [447, 97.1, 73, 175, 122, 174]
long_inercial = [120, 159, 133, 130, 101, 98]
giroradio_termico = [201, 470, 160, 266, 210, 145]
larmor = [45.9, 203, 61.5, 168, 75.9, 68.4]
larmor_normal = [49.8, 319, 68, 20.2, 28.9, 69.1]

plt.plot(tiempos, hmin, '.', color = 'C0')
plt.plot(tiempos, hmax, '.', color='C0')
for i in range(len(tiempos)):
    plt.vlines(x=tiempos[i], ymin=hmin[i], ymax = hmax[i], zorder=0, linewidth=2, color='C0')
plt.scatter(tiempos, long_inercial, zorder=1,  label='long_inercial', color = 'C1')
plt.scatter(tiempos, giroradio_termico, zorder=1, label = 'giroradio térmico', color='C2')
plt.scatter(tiempos, larmor, zorder=1, label = 'Larmor', color='C3')
plt.scatter(tiempos, larmor_normal, zorder=1, label = 'Larmor proyección normal', color='C4')
plt.legend()
plt.xlabel('Fecha (dd-mm-aa)')
plt.ylabel('Longitud (km)')
plt.title('Comparación de la longitud inercial y el giroradio de los protones \n con el ancho de la MPB')
plt.show()
