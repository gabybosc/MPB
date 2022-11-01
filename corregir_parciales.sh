
fdupes -rdN ./  # elimina el archivo repetido (chequea tipo y tamaño de archivo. Puede tener diferente nombre)
find . -empty -type d -delete  # elimina las carpetas vacías

for folder in */; do
  mv "$folder"/* .
done

rm -v !(*1*)  # borra todo lo que no tenga un 1 escrito
