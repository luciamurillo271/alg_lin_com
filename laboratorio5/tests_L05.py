def multiplicarMatrices(A,B):
    res = np.zeros((len(A),len(B[0])))
    for f in range(len(A)):
        for c in range(len(B[0])):
            sum = 0
            for fb in range(len(B)):
                sum += A[f][fb]*B[fb][c]
            res[f][c] = sum

    return res

def traspuesta(A):
    f,c = A.shape
    T = np.zeros((c,f))
    i = 0
    while i < c:
        T[i,:] = A[:,i]
        i+=1
    return T

def norma(x,p):
    norm = 0
    if p == "inf":
        for elem in x:
            if abs(elem) > norm:
                norm = abs(elem)
    else:
        suma = 0
        for elem in x:
            suma += abs(elem)**p
            norm = suma**(1/p)
    return norm

def producto_interno(x, y):
    if len(x) != len(y):
        return None
    prod = 0
    for i in range(len(x)):
        prod += x[i]*y[i]
    return prod

def inversa(A):
    inv_A = np.eye(len(A))
    if len(A) != len(A[0]):
        return None
    for c in range(len(A)):
        for f in range(c + 1, len(A)):
            m = A[f][c]/A[c][c]

'''
def Gram_Schmidt(col_A):
    Q = col_A.copy()

    for i in range(len(col_A)):
        for j in range(i):
            Q[i] -= (producto_interno(Q[i],Q[j])/norma(Q[j],2))*Q[j]

    for i in range(len(Q)):
        Q[i] = Q[i]/norma(Q[i],2)

    return Q

def QR_con_GS(A,tol=1e-12,retorna_nops=False):
    if len(A) != len(A[0]):
        return None
    
    columnas_A = traspuesta(A)
    Q = Gram_Schmidt(columnas_A)

    R = np.zeros((len(A),len(A)))

    R[0][0] = norma(columnas_A[0], 2)
    for c in range(len(R)):
        for f in range(len(R)):
            R
'''
    


### Funciones L05-QR

def QR_con_GS(A,tol=1e-12,retorna_nops=False):
    """
    A una matriz de n x n 
    tol la tolerancia con la que se filtran elementos nulos en R
    retorna_nops permite (opcionalmente) retornar el numero de operaciones realizado
    retorna matrices Q y R calculadas con Gram Schmidt (y como tercer argumento opcional, el numero de operaciones).
    Si la matriz A no es de n x n, debe retornar None
    """
    if len(A) != len(A[0]):
        return None

    Q = np.zeros((len(A), len(A)))
    R = np.zeros((len(A), len(A)))
    A_tras = traspuesta(A)
    
    no_descartados = []
    for j in range(len(A)):
        qj = A_tras[j].copy()
        for k in no_descartados:
            R[k][j] = producto_interno(Q[k],qj)
            qj -= R[k][j]*Q[k]
        R[j][j] = norma(qj, 2)   
        if norma(qj, 2) > tol:
            Q[j] = qj/R[j][j]
            no_descartados.append(j)
    
    Q_mono = traspuesta(Q[no_descartados])
    R_mono = R[:len(Q_mono[0])]
    return Q_mono, R_mono
    

def QR_con_HH(A,tol=1e-12,extras=False):
    """
    A una matriz de m x n (m>=n)
    tol la tolerancia con la que se filtran elementos nulos en R
    retorna matrices Q y R calculadas con reflexiones de Householder
    Si la matriz A no cumple m>=n, debe retornar None
    extras : bool, opcional
        Si es True, devuelve informacion extra sobre el proceso de factorizacion.
        Por defecto es False. Esto lo hacemos para poder graficar el proceso.
    Devuelve la factorizacion QR de A usando reflectores de Householder.
    Devuelve: 
        Q, R, extra_info (si extras es True)
        Q, R (si extras es False)
    extra_info es un diccionario con la clave:
        'R_matrices': lista de las matrices R en cada paso
        'Q_matrices': lista de las matrices Q en cada paso
    """
    if len(A) < len(A[0]):
        return None
    
    R = A.copy().astype(float)
    Q = np.eye(len(A))
    for k in range(len(A)):
        x = np.zeros(len(A) - k)
        for f in range(len(A) - k):
            x[f] = R[f + k][k]
        a = -x[0]*norma(x, 2)/abs(x[0])
        u = np.array(x.copy())
        u[0] -= a 
        if norma(u, 2) > tol:
            u = u/norma(u, 2)
            Hk = np.eye(len(A)-k)
            Hk -= 2*multiplicarMatrices(u.reshape((len(u),1)), traspuesta(u.reshape((len(u),1))))
            Hk_mono = np.eye(len(A))
            Hk_mono[k:,k:] = Hk
            R = multiplicarMatrices(Hk_mono, R)
            Q = multiplicarMatrices(Q, traspuesta(Hk_mono))
    return Q, R

def calculaQR(A,metodo='RH',tol=1e-12):
    """
    A una matriz de n x n 
    tol la tolerancia con la que se filtran elementos nulos en R    
    metodo = ['RH','GS'] usa reflectores de Householder (RH) o Gram Schmidt (GS) para realizar la factorizacion
    retorna matrices Q y R calculadas con Gram Schmidt (y como tercer argumento opcional, el numero de operaciones)
    Si el metodo no esta entre las opciones, retorna None
    """
    if metodo == 'RH':
        Q, R = QR_con_HH(A, tol)
    elif metodo == 'GS':
        Q, R = QR_con_GS(A,tol)
    else:
        return None
    return Q, R



# Tests L05-QR:

import numpy as np

# --- Matrices de prueba ---
A2 = np.array([[1., 2.],
               [3., 4.]])

A3 = np.array([[1., 0., 1.],
               [0., 1., 1.],
               [1., 1., 0.]])

A4 = np.array([[2., 0., 1., 3.],
               [0., 1., 4., 1.],
               [1., 0., 2., 0.],
               [3., 1., 0., 2.]])

# --- Funciones auxiliares para los tests ---
def check_QR(Q,R,A,tol=1e-10):
    # Comprueba ortogonalidad y reconstrucción
    assert np.allclose(Q.T @ Q, np.eye(Q.shape[1]), atol=tol)
    assert np.allclose(Q @ R, A, atol=tol)

# --- TESTS PARA QR_by_GS2 ---
Q2,R2 = QR_con_GS(A2)
check_QR(Q2,R2,A2)

Q3,R3 = QR_con_GS(A3)
check_QR(Q3,R3,A3)

Q4,R4 = QR_con_GS(A4)
check_QR(Q4,R4,A4)

# --- TESTS PARA QR_by_HH ---
Q2h,R2h = QR_con_GS(A2)
check_QR(Q2h,R2h,A2)

Q3h,R3h = QR_con_HH(A3)
check_QR(Q3h,R3h,A3)

Q4h,R4h = QR_con_HH(A4)
check_QR(Q4h,R4h,A4)

# --- TESTS PARA calculaQR ---
Q2c,R2c = calculaQR(A2,metodo='RH')
check_QR(Q2c,R2c,A2)

Q3c,R3c = calculaQR(A3,metodo='GS')
check_QR(Q3c,R3c,A3)

Q4c,R4c = calculaQR(A4,metodo='RH')
check_QR(Q4c,R4c,A4) 