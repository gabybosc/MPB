import numpy as np
import matplotlib.pyplot as plt
import cdflib as cdf
from funciones import find_nearest, unix_to_decimal, unix_to_timestamp
from funciones_plot import plot_select
from matplotlib.colors import LogNorm
import datetime as dt
import matplotlib.dates as md
plt.ion()
# np.set_printoptions(precision=4)

path = '../../datos/SWEA/' #path a los datos
cdf_file = cdf.CDF(path + 'mvn_swe_l2_svyspec_20160316_v04_r01.cdf')

flux_all = cdf_file.varget('diff_en_fluxes')
energia = cdf_file.varget('energy')
t_unix = cdf_file.varget('time_unix')

tu = unix_to_decimal(t_unix)
ti = np.where(tu == find_nearest(tu, 17.85))[0][0]
tf = np.where(tu == find_nearest(tu, 18.5))[0][0]
t = tu[ti:tf]
flux = flux_all[ti:tf]

datenums = unix_to_timestamp(t_unix[ti:tf])

log_flux = np.flip(np.log(flux), axis=1)
log_flux[log_flux<-1000] = None# np.min(log_flux[log_flux>-1000])

flux_plot = np.transpose(flux)[::-1]

plt.figure()
plt.subplots_adjust(bottom=0.2)
plt.xticks( rotation=25 )
ax = plt.gca()
xfmt = md.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(xfmt)
plt.imshow(flux_plot, aspect = 'auto',origin = 'lower', extent=(datenums[0], datenums[-1],  energia[-1], energia[0]), cmap='jet', norm=LogNorm(vmin=1E4, vmax=10**(9.5)))
clb = plt.colorbar()
clb.ax.set_title('log')
plt.xlabel('Tiempo')
plt.ylabel('Energía (eV)')
plt.title('SWEA')
