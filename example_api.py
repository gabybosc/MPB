import gspread
from oauth2client.service_account import ServiceAccountCredentials
from pprint import pprint

scope = ["https://spreadsheets.google.com/feeds","https://www.googleapis.com/auth/spreadsheets","https://www.googleapis.com/auth/drive.file","https://www.googleapis.com/auth/drive"]

creds = ServiceAccountCredentials.from_json_keyfile_name("google_api_MPB.json", scope)

client = gspread.authorize(creds)

sheet = client.open("MPB").sheet1

data = sheet.get_all_records()

row = sheet.row_values(1) #selecciona la fila 1

col = sheet.col_values(4) #selecciona la col 4

cell = sheet.cell(2,3).value #selecciona la celda de la fila 2, columna 3

sheet.update_acell('B4', 'hola')

# sheet.insert_row(row, 4) #en la fila 4 pone lo que tengo en la variable row
# sheet.delete_row(4) #ahora borra la fila 4
# sheet.update_cell(4,2, 'Reescrito lince') #sobreescribe la celda

# pprint(cell)
