import numpy as np
def multiplicarMatrices(A,B):
    res = np.zeros((len(A),len(B[0])))
    for f in range(len(A)):
        for c in range(len(B[0])):
            sum = 0
            for fb in range(len(B)):
                sum += A[f][fb]*B[fb][c]
            res[f][c] = sum

    return res

def multiplicarAx(A,x):
    Ax = np.zeros(len(x))
    for i in range(len(A)):
        sum = 0
        for j in range(len(x)):
            sum += A[i][j]*x[j]
        Ax[i] = sum
    return Ax

def traspuesta(A):
    f,c = A.shape
    T = np.zeros((c,f))
    i = 0
    while i < c:
        T[i,:] = A[:,i]
        i+=1
    return T

def esSimetrica(A,error):
    f,c = A.shape
    if f == c:
        T = traspuesta(A)
        simetrica = True
        i = 0
        while simetrica and i < f:
            for j in range(c):
                simetrica = simetrica and abs(A[i,j] - T[i,j]) < error
            i += 1
        
        return simetrica

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

def f_A(A, v, k):
    w = v.copy()
    for i in range(k):
        w = multiplicarAx(A,w)
        if norma(w, 2) > 0:
            w = w/norma(w, 2)
    return w

def metpot2k(A,tol=1e-15,K=1000):
    v = np.random.randn(len(A))
    v_m = f_A(A,v,2)
    e = producto_interno(v_m, v)
    k = 0
    while abs(e - 1)> tol and k < K:
        v = v_m
        v_m = f_A(A,v,2)
        e = producto_interno(v_m,v)
        k += 1
    Av_m = multiplicarAx(A, v_m)
    lamda = producto_interno(v_m, Av_m)
    error = e - 1
    return v_m, lamda, k

def diagRH(A,tol=1e-15,K=1e5):
    if not esSimetrica(A, tol):
        return None
    v, lamda, error = metpot2k(A, tol, K)
    if lamda<tol:
        return np.zeros((len(A),len(A))), np.zeros((len(A),len(A)))
    I = np.eye(len(A))
    u = I[0].copy()
    u -= v
    H = I - (2/norma(u,2)**2)*multiplicarMatrices(u.reshape((len(u),1)), traspuesta(u.reshape((len(u),1))))
    if len(A) == 2:
        S = H
        D = multiplicarMatrices(H, multiplicarMatrices(A, traspuesta(H)))
    else:
        B = multiplicarMatrices(H, multiplicarMatrices(A, traspuesta(H)))
        A_mon = B[1:, 1:]
        S_mon, D_mon = diagRH(A_mon, tol, K)
        D = np.eye(len(A))
        D[0][0] = lamda
        D[1:, 1:] = D_mon
        S1 = np.eye(len(A))
        S1[1:, 1:] = S_mon
        S = multiplicarMatrices(H, S1)    
    return S, D


def svd_reducida(A,k="max",tol=1e-15):
    m , n = A.shape
    A_t = traspuesta(A)
    if (m<n):
        M = multiplicarMatrices(A,A_t) # M e mxm
    else:
        M = multiplicarMatrices(A_t,A) # M e nxn
    avec, aval = diagRH(M)
    aval_no_nulos = []
    i=0
    while i<len(aval[0]) and aval[i,i]>tol:
        aval_no_nulos.append(aval[i,i]**(1/2))
        i+=1
    r = len(aval_no_nulos)
    if(k=="max"):
        k=r
    p = min(k,r)
    hat_Sig = np.zeros((p,1)) #se podria cortar diag hasta p
    for i in range (p): #eficiencia? pasar de lista a columna
        hat_Sig[i,0]=aval_no_nulos[i]

    B = np.zeros((m,p))
    for i in range(p):
        B[:,i] = avec[:, i]/norma(avec[:, i], 2)
    C = np.zeros((n,p))
    for i in range(p):
        C[:,i] = (multiplicarAx(A_t,B[:,i]))/hat_Sig[i,0]
    if (m<n):
        hat_U = B
        hat_V = C
    else:
        hat_U = C
        hat_V = B
    return hat_U, hat_Sig, hat_V

# Tests L08
    
# Matrices al azar
def genera_matriz_para_test(m,n=2,tam_nucleo=0):
    if tam_nucleo == 0:
        A = np.random.random((m,n))
    else:
        A = np.random.random((m,tam_nucleo))
        A = np.hstack([A,A])
    return(A)

def test_svd_reducida_mn(A,tol=1e-15):
    m,n = A.shape
    hU,hS,hV = svd_reducida(A,tol=tol)
    nU,nS,nVT = np.linalg.svd(A)
    r = len(hS)+1
    assert np.all(np.abs(np.abs(np.diag(hU.T @ nU))-1)<10**r*tol), 'Revisar calculo de hat U en ' + str((m,n))
    assert np.all(np.abs(np.abs(np.diag(nVT @ hV))-1)<10**r*tol), 'Revisar calculo de hat V en ' + str((m,n))
    assert len(hS) == len(nS[np.abs(nS)>tol]), 'Hay cantidades distintas de valores singulares en ' + str((m,n))
    assert np.all(np.abs(hS-nS[np.abs(nS)>tol])<10**r*tol), 'Hay diferencias en los valores singulares en ' + str((m,n))

for m in [2,5,10,20]:
    for n in [2,5,10,20]:
        for _ in range(10):
            A = genera_matriz_para_test(m,n)
            test_svd_reducida_mn(A)


# Matrices con nucleo

m = 12
for tam_nucleo in [2,4,6]:
    for _ in range(10):
        A = genera_matriz_para_test(m,tam_nucleo=tam_nucleo)
        test_svd_reducida_mn(A)

# Tamaños de las reducidas
A = np.random.random((8,6))
for k in [1,3,5]:
    hU,hS,hV = svd_reducida(A,k=k)
    assert hU.shape[0] == A.shape[0], 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[0] == A.shape[1], 'Dimensiones de hV incorrectas(caso a)'
    assert hU.shape[1] <= k, 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[1] <= k, 'Dimensiones de hV incorrectas(caso a)'
    assert len(hS) <= k, 'Tamaño de hS incorrecto'
