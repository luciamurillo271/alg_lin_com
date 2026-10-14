# Test L06-metpot2k, Aval

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
                simetrica = np.isclose(A[i,j], T[i,j],atol=error)
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

def normaExacta(A,p=[1,"inf"]):
    if isinstance(p,int):
        if p == 1:
            A_traspuesta = traspuesta(A)
            norm = normaExacta(A_traspuesta,"inf")
        else:
            norm = None
    elif isinstance(p,str):
        if p == "inf":
            norm = 0
            for fila in A:
                norma_fila = norma(fila,1)
                if norma_fila > norm:
                    norm = norma_fila
        else:
            norm = None
    else:
        if p[0] == 1:
            norm = [normaExacta(A,1),normaExacta(A,"inf")]
        else:
            norm = None
    return norm

def producto_interno(x, y):
    if len(x) != len(y):
        return None
    prod = 0
    for i in range(len(x)):
        prod += x[i]*y[i]
    return prod

#### EJERCICIOS
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

#%%--------------------------------------------------------------
#ejercicio 2
import matplotlib.pyplot as plt

def ej_metpot2k(A,tol=1e-15,K=1000):
    v = np.random.randn(len(A))
    v_m = f_A(A,v,2)
    e = producto_interno(v_m, v)
    k = 0
    est_lamda = []
    while abs(e - 1)> tol and k < K:
        Av_m = multiplicarAx(A, v_m)
        a_k = producto_interno(v_m, Av_m)
        est_lamda.append(a_k)
        v = v_m
        v_m = f_A(A,v,2)
        e = producto_interno(v_m,v)
        k += 1
    Av_m = multiplicarAx(A, v_m)
    lamda = producto_interno(v_m, Av_m)
    error = e - 1
    return v_m, lamda,error, est_lamda

C = np.random.rand(100,100)
A = (1/2)*(C + traspuesta(C))
I = np.eye(100)

B = A + 500*I

v_m, lamda, error, est_lamda = ej_metpot2k(B,1e-8,1000)
for i in range(len(est_lamda)):
    est_lamda[i] = abs(lamda - est_lamda[i])

plt.semilogy(est_lamda)      # el último es 0 y rompe la escala log
plt.xlabel("k")
plt.ylabel(r"$|\lambda - \lambda_k|$")
plt.grid(True, which="both", alpha=0.3)
plt.show()

#%%-----------------------------------------------------------------------------
#### TESTEOS
# Tests metpot2k

S = np.vstack([
    np.array([2,1,0])/np.sqrt(5),
    np.array([-1,2,5])/np.sqrt(30),
    np.array([1,-2,1])/np.sqrt(6)
              ]).T

# Pedimos que pase el 95% de los casos
exitos = 0
for i in range(100):
    D = np.diag(np.random.random(3)+1)*100
    A = S@D@S.T
    v,l,_ = metpot2k(A,1e-15,1e5)
    if np.abs(l - np.max(D))< 1e-8:
        exitos += 1
assert exitos > 95


#Test con HH
exitos = 0
for i in range(100):
    v = np.random.rand(9)
    #v = np.abs(v)
    #v = (-1) * v
    ixv = np.argsort(-np.abs(v))
    D = np.diag(v[ixv])
    I = np.eye(9)
    H = I - 2*np.outer(v.T, v)/(np.linalg.norm(v)**2)   #matriz de HouseHolder

    A = H@D@H.T
    v,l,_ = metpot2k(A, 1e-15, 1e5)
    #max_eigen = abs(D[0][0])
    if abs(l - D[0,0]) < 1e-8:         
        exitos +=1
assert exitos > 95



# Tests diagRH
D = np.diag([1,0.5,0.25])
S = np.vstack([
    np.array([1,-1,1])/np.sqrt(3),
    np.array([1,1,0])/np.sqrt(2),
    np.array([1,-1,-2])/np.sqrt(6)
              ]).T

A = S@D@S.T
SRH,DRH = diagRH(A,tol=1e-15,K=1e5)
assert np.allclose(D,DRH)
assert np.allclose(np.abs(S.T@SRH),np.eye(A.shape[0]),atol=1e-7)
#print(np.round(D, 4))
#print(np.round(DRH, 4))
#print(np.abs(D - DRH).max())


# Pedimos que pase el 95% de los casos
exitos = 0
for i in range(100):
    A = np.random.random((5,5))
    A = 0.5*(A+A.T)
    S,D = diagRH(A,tol=1e-15,K=1e5)
    ARH = S@D@S.T
    e = normaExacta(ARH-A,p='inf')
    if e < 1e-5: 
        exitos += 1
assert exitos >= 95



