import numpy as np

def hallar_minimo(mag, orbita, int una_vuelta, calendario, R):
    cdef int paso
    orbitas = [orbita[:una_vuelta], orbita[una_vuelta:una_vuelta*2], orbita[una_vuelta*2:una_vuelta*3], orbita[una_vuelta*3:una_vuelta*4], orbita[una_vuelta*4:]]
    resta = np.zeros((len(R),3))
    paso = 50

    for indice, l in enumerate(orbitas):
        calendario[i*5+indice,0] = mag[1,1]
        calendario[i*5+indice,1] = indice+1 
        pos = l * 3390
        X_MSO = pos[:, 0]
        Z_MSO = pos[:, 2]
        idx_min = np.zeros(int(una_vuelta/paso))
        max_acercamiento = np.zeros(int(una_vuelta/paso))
        minimo = 0
        for k in range(0,int(una_vuelta)-100, paso):
            if Z_MSO[k] > 0 and X_MSO[k] > 0:
                for m in range(len(R)):
                    resta[m, :] = l[k,:] - R[m,:]
                A = np.linalg.norm(resta, axis=1)
                idx_min[int(k/paso)] = np.argmin(A)
                max_acercamiento[int(k/paso)] = A[int(idx_min[int(k/paso)])]
        if sum(max_acercamiento) == 0: #si es cero, va a fallar todo el script, así que digo que esa órbita es mala y listo
            calendario[i*5+indice,2] = 0
            calendario[i*5+indice, 3] = 0
            calendario[i*5+indice,4] = 0
        else:
            minimo = np.where( max_acercamiento==np.min(max_acercamiento[np.nonzero(max_acercamiento)]))[0][0] #busca el minimo que no sea cero

        idx = minimo * paso

        altitud = np.linalg.norm(pos[int(idx),:]) - 3390

        if altitud < 1300 and altitud > 300:
            calendario[i*5+indice,2] = 1

        SZA = np.arccos(np.clip(np.dot(pos[int(idx)]/np.linalg.norm(pos[int(idx)]), [1,0,0]), -1.0, 1.0))* 180/np.pi
        # print(j, indice, SZA)
        if SZA < 30:
            calendario[i*5+indice,3] = 1
        elif SZA > 60:
            calendario[i*5+indice, 4] = 1
        else:
            calendario[i*5+indice, 5] = 1

        return(calendario)
