path = "../../../datos/simulacion_chuanfei/"
file = open(path + "z=0_HallOn.out", 'rb')


byte = file.read(1)

while byte:

    print(byte)

    byte = file.read(1)

file.close()
