# MPB

Python files for analysing the MPB

READ CDF:
Con cdflib, para ver todas las variables hacer archivo_cdf.cdf_info()
Para abrir una variable: archivo_cdf.varget('nombre de la variable')


UNPICKLE:
with open('../outputs/figs_MPB/MPB_y2016_d336_t11.pkl', 'rb') as file:
     to_file = pkl.load(file)
plt.show()

PLOTTING:

clweb/normales.py: plotea las normales y el fit

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

giroradio_termico.py: Calcula el giroradio térmico (es decir, a partir de los datos de temperatura) y la long inercial.

giroradio_in/outbound.py: Calcula el giroradio (usando v_perp al campo B y también usando la velocidad proyectada en la normal) y la long inercial.
                          Versión para cruce inbound y outbound.

Ehall_convectivo.py: Calcula el EHall a partir del Ecv.

Ehall_jxB.py: Calcula EHall a partir del rotor del campo B: E ~ rotB x B usando la ecuación con B_para y B_perp

Escalas_lambda.py: Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas


CLWEB:
Analisis_in/outbound / Analisis_sin_spreadsheet: Hacen el análisis completo (todo lo que está en el excel) para órbitas in/outbound. Plotea el hodograma y las normales.

escalas_lambda: Hace un barrido en ventanas para ver cuándo se maximiza el cociente de lambdas

filtro.py: Hace un filtro butterworth para quitar el ruido de la señal del magnetómetro que es de aproximadamente 180 ms.
            Se puede usar para filtrar más cosas, pero esto se hace siempre.

funciones_mva: tiene funciones varias para hacer el mva y el análisis de vignes.
funciones_plot: Tiene varias funciones que uso en los plots en general.
funciones: Funciones tanto que uso en el MVA como funciones más generales.

importar_datos: Tiene las funciones para importar los datos de mag, swea, swia y lpw. La gracia es para ahorrarme siempre las mismas líneas de código en el comienzo de los scripts.

MVA_hires / MVA_sin_spreadsheet: La función que hace el MVA que después llama el Analisis_inbound /etc.

normales: plotea la MPB de vignes y las diferentes normales, la del MVA y del fit, para varios cruces.

planetocentric: agarra los datos pc y calcula la latitud y longitud del cruce para comparar después con el mapa de connerney 2000 de los mapas corticales.

plot_escalas_lambda.py: plotea los archivos que devuelve escalas_lambda.py

plot_mag.py: Plotea los datos de B en coordenadas esféricas. Sirve para descartar campos corticales.

plot_mpb.py: Este script plotea mag, swea, swia y lpw en la región de interés. Es la fig principal del grl.

plot_orbitas: Hace la fig2 del poster de la AGU2019. Básicamente pone en una proyección 2D la órbita en torno a Marte.

plot_seleccionar_MPB.py: plotea los datos MAG, SWEA, SWIA, LPW y me da un cursor interactivo para poder elegir los tiempos t1t2t3t4.

plot_tiempos.py: plotea los datos de MAG, SWEA, SWIA, LPW con los tiempos t1t2t3t4 y los guarda como pickle. Para leerlos: usar outputs/leer_figs_MPB.py
