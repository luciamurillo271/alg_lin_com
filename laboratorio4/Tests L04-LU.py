#from moduloALC import calculaLU, res_tri, inversa, calculaLDV, esSDP
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

def calculaLU(A):
    cant_op = 0

    if A is None:
        return None, None, cant_op
    m, n = A.shape
    Ac = A.copy().astype(float)
    
    if m!=n:
        return None,None,cant_op

    for col in range(n):
        pivo = Ac[col][col]
        if abs(pivo) < 1e-15:
            return None, None, 0

        for fil in range(col+1,m):
            if Ac[fil][col] != 0:
                Ac[fil][col] = Ac[fil][col]/pivo
                cant_op +=1

                Ac[fil][col+1:] -= Ac[col][col+1:]*(Ac[fil][col])

        cant_op += 2*((n - col-1)**2)

    L = np.eye(n)
    U = np.zeros((m,n))
    for f in range(m):
        for c in range(n):
            if f > c:
                L[f][c] = Ac[f][c]
            else:
                U[f][c] = Ac[f][c]
    
    return L, U, cant_op


def es_inferior(A):
    inferior = True
    for i in range(len(A)):
        for j in range(i+1,len(A[0])):
            inferior = inferior and A[i][j] == 0
    return inferior


def res_tri2(L,b,inferior=True):
    x = np.zeros(len(L))
    for i in range(len(L)):
        sum = 0
        for j in range(i):
            sum += L[i][j]*x[j]
        x[i] = (b[i] - sum)/L[i][i]
        
    return x

def res_tri(L,b,inferior=True):
    x = np.zeros(len(L))
    if not inferior:
        for i in range(len(L)-1,-1,-1):
            sum = 0
            for j in range(len(L)-1,i,-1):
                sum += L[i][j]*x[j]
            x[i] = (b[i] - sum)/L[i][i]
        return x
    for i in range(len(L)):
        sum = 0
        for j in range(i):
            sum += L[i][j]*x[j]
        x[i] = (b[i] - sum)/L[i][i]
        
    return x


def inversa(A):
    L, U, nops = calculaLU(A)
    if L is None or U is None:
        return None
    L_inv = np.eye(len(L))
    for i in range(len(L[0])):
        piv = L[i][i]
        for j in range(i+1,len(L)):
            l = L[j][i]/piv
            L[j][:i+1] -= L[i][:i+1]*l
            L_inv[j][:i+1] -= L_inv[i][:i+1]*l

    U_inv = np.eye(len(U))
    for i in range(len(U[0])-1,-1,-1):
        if U[i][i] != 1:
            U_inv[i] = U_inv[i]/U[i][i]
            U[i] = U[i]/U[i][i]
        piv = U[i][i]
        for j in range(i):
            u = U[j][i]/piv
            U[j][i:] -= U[i][i:]*u
            U_inv[j][i:] -= U_inv[i][i:]*u

    
    return multiplicarMatrices(U_inv,L_inv)

def traspuesta(A):
    f,c = A.shape
    T = np.zeros((c,f))
    i = 0
    while i < f:
        T[i:,:] = A[:,i]
        i+=1
    return T

def calculaLDV(A): #mision cumplida
    L,U,nops = calculaLU(A)
    if L is None or U is None:
        return None
    U = traspuesta(U)
    V,D, nops = calculaLU(U)
    V = traspuesta(V)

    return L,D,V

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

def esSDP(A, atol=1e-8):
    if not esSimetrica(A,atol):
        return False
    if calculaLDV(A) is None:
        return False
    L,D,V = calculaLDV(A)
    diag_pos = True
    for i in range(len(D)):
        diag_pos = diag_pos and D[i][i] > 0
    return diag_pos


def calculaCholesky(A,atol=1e-10): #mision cumplida
    if not esSDP(A,atol):
        return None
    L,D,V = calculaLDV(A)

    for i in range(len(D)):
        D[i][i] = np.sqrt(D[i][i])

    R = multiplicarMatrices(L,D)
    return R

# TESTS L04-LU

# TESTS LU
print("TESTS calculaLU")

L0 = np.array([[1,0,0],
               [0,1,0],
               [1,1,1]])

U0 = np.array([[10,1,0],
               [0,2,1],
               [0,0,1]])

A =  L0 @ U0
L,U,nops = calculaLU(A)
assert(np.allclose(L,L0))
assert(np.allclose(U,U0))


L0 = np.array([[1,0,0],
               [1,1.001,0],
               [1,1,1]])

U0 = np.array([[1,1,1],
               [0,1,1],
               [0,0,1]])
A =  L0 @ U0
L,U,nops = calculaLU(A)
assert(not np.allclose(L,L0))
assert(not np.allclose(U,U0))
assert(np.allclose(L,L0,atol=1e-3))
assert(np.allclose(U,U0,atol=1e-3))
assert(nops == 13)

L0 = np.array([[1,0,0],
               [1,1,0],
               [1,1,1]])

U0 = np.array([[1,1,1],
               [0,0,1],
               [0,0,1]])

A =  L0 @ U0
L,U,nops = calculaLU(A)
assert(L is None)
assert(U is None)
assert(nops == 0)

assert(calculaLU(None) == (None, None, 0))

assert(calculaLU(np.array([[1,2,3],[4,5,6]])) == (None, None, 0))

print("-----ÉXITO!!!!\n")


## TESTS res_tri
print("TESTS res_tri")

A = np.array([[1,0,0],
              [1,1,0],
              [1,1,1]])

b = np.array([1,1,1])
assert(np.allclose(res_tri(A,b),np.array([1,0,0])))

b = np.array([0,1,0])
assert(np.allclose(res_tri(A,b),np.array([0,1,-1])))

b = np.array([-1,1,-1])
assert(np.allclose(res_tri(A,b),np.array([-1,2,-2])))

b = np.array([-1,1,-1])
assert(np.allclose(res_tri(A,b,inferior=False),np.array([-1,1,-1])))

A = np.array([[3,2,1],[0,2,1],[0,0,1]])
b = np.array([3,2,1])
assert(np.allclose(res_tri(A,b,inferior=False),np.array([1/3,1/2,1])))

A = np.array([[1,-1,1],[0,1,-1],[0,0,1]])
b = np.array([1,0,1])
assert(np.allclose(res_tri(A,b,inferior=False),np.array([1,1,1])))
print("-----ÉXITO!!!!\n")


# Test inversa
print("TESTS inversa")

def esSingular(A):
    try:
        np.linalg.inv(A)
        return False
    except:
        return True

# Por que no siempre es invertible, hacemos varios tests
ntest = 10
for i in range(ntest):
    A = np.random.random((4,4))
    A_ = inversa(A)
    if not esSingular(A):
        inversaConNumpy = np.linalg.inv(A)
        assert(A_ is not None)
        assert(np.allclose(inversaConNumpy,A_))
    else: 
        assert(A_ is None)

# Matriz singular devería devolver None
A = np.array([[1,2,3],[4,5,6],[7,8,9]])
assert(inversa(A) is None)

print("-----ÉXITO!!!!\n")



# Test LDV:
print("TESTS calculaLDV")

L0 = np.array([[1,0,0],[1,1.,0],[1,1,1]])
D0 = np.diag([1,2,3])
V0 = np.array([[1,1,1],[0,1,1],[0,0,1]])
A =  L0 @ D0 @ V0
L,D,V = calculaLDV(A)
assert(np.allclose(L,L0))
assert(np.allclose(D,D0))
assert(np.allclose(V,V0))


L0 = np.array([[1,0,0],[1,1.001,0],[1,1,1]])
D0 = np.diag([3,2,1])
V0 = np.array([[1,1,1],[0,1,1],[0,0,1.001]])
A =  L0 @ D0  @ V0
L,D,V = calculaLDV(A)
assert(np.allclose(L,L0,1e-3))
assert(np.allclose(D,D0,1e-3))
assert(np.allclose(V,V0,1e-3))

print("-----ÉXITO!!!!\n")

# TESTS SDP
print("TESTS esSDP")

L0 = np.array([[1,0,0],[1,1,0],[1,1,1]])
D0 = np.diag([1,1,1])
A = L0 @ D0 @ L0.T
assert(esSDP(A))

D0 = np.diag([1,-1,1])
A = L0 @ D0 @ L0.T
assert(not esSDP(A))

D0 = np.diag([1,1,1e-16])
A = L0 @ D0 @ L0.T
assert(not esSDP(A))

L0 = np.array([[1,0,0],
               [1,1,0],
               [1,1,1]])
D0 = np.diag([1,1,1])
V0 = np.array([[1,0,0],
               [1,1,0],
               [1,1+1e-3,1]]).T
A = L0 @ D0 @ V0
assert(esSDP(A,1e-3))

print("-----ÉXITO!!!!\n")
print("---FINALIZADO LABO 4!---")

