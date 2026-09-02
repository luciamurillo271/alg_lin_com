import numpy as np
def error(x,y):
    return abs(np.float64(x) - np.float64(y))
'''Recibe dos numeros x e y, y calcula el error de aproximar x usando y en float64'''

def error_relativo(x,y):
    return abs(error(x,y)) / abs(np.float64(x))
'''Recibe dos numeros x e y, y calcula el error relativo de aproximar x usando y en float64'''
    
def matricesIguales(A, B):
    iguales = True
    f,c = A.shape
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        iguales = False
    else:
        i = 0
        while iguales and i < f:
            for j in range(c):
                iguales = iguales and np.isclose(A[i][j], B[i][j])
            i += 1
    return iguales
'''Devuelve True si ambas matrices son iguales y False en otro caso. Considerar que las matrices pueden tener distintas dimensiones, ademas de distintos valores.'''