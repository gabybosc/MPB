# Chuanfei

Si tengo los archivos binarios .out y el IDL, usando write_mpb.pro + un archivo que
ahora está perdido, se obtienen todos los .csv: posicion, densidades, velocidades,
campo B, presion, corrientes, grad_p.

Para juntar todo eso en un único archivo usamos **convertir_en_un_txt.py**, que
me va a devolver un .gz donde están ordenados los resultados. Hay que hacer esto
para y=0, z=0 (x=0 no lo uso pero también se puede). Puedo ya usarlos.

Si necesitamos, en particular si vemos algunas cosas raras, conviene usar
**recortar.py**. Me devuelve el archivo recortado en la zona que elija.
Por default está en modo que recorte para tener sólo valores positivos de z/y
