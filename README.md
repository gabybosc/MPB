# MPB

Python files for analysing the MPB

READ CDF:
Con cdflib, para ver todas las variables hacer variable_cdf.cdf_info()

PLOTTING:

plot_tiempos.py: plotea los datos de MAG, SWEA, SWIA, LPW con los tiempos t1t2t3t4 y los guarda como pickle. Para leerlos: usar outputs/leer_figs_MPB.py

plot_seleccionar_MVA.py: plotea los datos MAG, SWEA, SWIA, LPW y me da un cursor interactivo para poder elegir los tiempos t1t2t3t4.

plot_mpb.py: la figura principal de mi tesis. Plotea MAG, SWEA, SWIA y LPW en toda la región y de nuevo en un zoom.

plot_mag.py: plotea MAG en alta resolución.

plot_escalas_lambda.py: plotea los archivos que devuelve escalas_lambda.py

MVA:

MVA_automatico.py: hace el MVA a partir del archivo de tiempos. Devuelve un hodograma y el cociente de lambdas. Está bueno como primera aproximación o si estoy segura de que los tiempos t1t2t3t4 son los que quiero.

MVA_hires.py: Es el script completo para datos de alta resolución. Los filtra.

MVA_lowres.py: Hace todo para datos de baja resolución.


OTROS ANÁLISIS:
bup_bdown.py: grafica la variación de Bup y Bdown en el tiempo.

densidades.py: ?

Ehall_convectivo.py: Calcula el EHall a partir del Ecv.

Ehall_jxB.py: Calcula EHall a partir del rotor del campo B: E ~ rotB x B usando la ecuación con B_para y B_perp

Escalas_lambda.py: Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas
