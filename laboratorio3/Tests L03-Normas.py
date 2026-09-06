# Tests L03-Normas
import numpy as np
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

def normaliza(X,p):
    normas = []
    for i in range(len(X)):
        vector = X[i].astype(float).copy()
        for j in range(len(vector)):
            vector[j] = vector[j]/norma(X[i],p)
        normas.append(vector)
    return normas

def normaMatMC(A,q,p,Np):
    angulos = np.linspace(0,2*np.pi,Np)
    vectores = np.array([[np.cos(a),np.sin(a)] for a in angulos])
    normalizados = normaliza(vectores,p)
    max = []
    norm = 0
    for vector in normalizados:
        norma_matriz = norma(A@vector,q)
        if norma_matriz > norm:
            norm = norma_matriz
            max = vector
    return norm,max

def traspuesta(A):
    f,c = A.shape
    T = np.zeros((c,f))
    i = 0
    while i < f:
        T[i:,:] = A[:,i]
        i+=1
    return T

def normaExacta(A,p=[1,"inf"]):
    if isinstance(p,int):
        if p == 1:
            norm = 0
            for j in range(len(A[0])):
                norma_columna = 0
                for i in range(len(A)):
                    norma_columna += abs(A[i][j])
                if norma_columna > norm:
                    norm = norma_columna
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

def condMC(A,p,Np):
    A_inversa = np.linalg.inv(A)
    nro_cond = normaMatMC(A,p,p,Np)[0]*normaMatMC(A_inversa,p,p,Np)[0]
    return nro_cond

def condExacta(A,p):
    A_inversa = np.linalg.inv(A)
    nro_cond = normaExacta(A,p)*normaExacta(A_inversa,p)
    return nro_cond

# Tests norma
assert(np.allclose(norma(np.array([1,1]),2),np.sqrt(2)))
assert(np.allclose(norma(np.array([1]*10),2),np.sqrt(10)))
assert(norma(np.random.rand(10),2)<=np.sqrt(10))
assert(norma(np.random.rand(10),2)>=0)


# Tests normaliza
for x in normaliza([np.array([1]*k) for k in range(1,11)],2):
    assert(np.allclose(norma(x,2),1))
for x in normaliza([np.array([1]*k) for k in range(2,11)],1):
    print(not np.allclose(norma(x,2),1) )
for x in normaliza([np.random.rand(k) for k in range(1,11)],'inf'):
    assert( np.allclose(norma(x,'inf'),1) )


# Tests normaExacta

assert(np.allclose(normaExacta(np.array([[1,-1],[-1,-1]]),1),2))
assert(np.allclose(normaExacta(np.array([[1,-2],[-3,-4]]),1),6))
assert(np.allclose(normaExacta(np.array([[1,-2],[-3,-4]]),'inf'),7))
assert(normaExacta(np.array([[1,-2],[-3,-4]]),2) is None)
assert(normaExacta(np.random.random((10,10)),1)<=10)
assert(normaExacta(np.random.random((4,4)),'inf')<=4)

# Test normaMC

nMC = normaMatMC(A=np.eye(2),q=2,p=1,Np=100000)
assert(np.allclose(nMC[0],1,atol=1e-3))
assert(np.allclose(np.abs(nMC[1][0]),1,atol=1e-3) or np.allclose(np.abs(nMC[1][1]),1,atol=1e-3))
assert(np.allclose(np.abs(nMC[1][0]),0,atol=1e-3) or np.allclose(np.abs(nMC[1][1]),0,atol=1e-3))

nMC = normaMatMC(A=np.eye(2),q=2,p='inf',Np=100000)
assert(np.allclose(nMC[0],np.sqrt(2),atol=1e-3))
assert(np.allclose(np.abs(nMC[1][0]),1,atol=1e-3) and np.allclose(np.abs(nMC[1][1]),1,atol=1e-3))

A = np.array([[1,2],[3,4]])
nMC = normaMatMC(A=A,q='inf',p='inf',Np=1000000)
assert(np.allclose(nMC[0],normaExacta(A,'inf'),rtol=2e-1)) 

# Test condMC

A = np.array([[1,1],[0,1]])
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaMatMC(A,2,2,10000)
normaA_ = normaMatMC(A_,2,2,10000)
condA = condMC(A,2,10000)
assert(np.allclose(normaA[0]*normaA_[0],condA,atol=1e-3))

A = np.array([[3,2],[4,1]])
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaMatMC(A,2,2,10000)
normaA_ = normaMatMC(A_,2,2,10000)
condA = condMC(A,2,10000)
assert(np.allclose(normaA[0]*normaA_[0],condA,atol=1e-3))

# Test condExacta

A = np.random.rand(10,10)
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaExacta(A,1)
normaA_ = normaExacta(A_,1)
condA = condExacta(A,1)
assert(np.allclose(normaA*normaA_,condA))

A = np.random.rand(10,10)
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaExacta(A,'inf')
normaA_ = normaExacta(A_,'inf')
condA = condExacta(A,'inf')
assert(np.allclose(normaA*normaA_,condA))
