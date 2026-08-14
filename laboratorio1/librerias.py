import numpy as np

def esCuadrada(A):
    return A.shape[0] == A.shape[1]



def triangSup(A):
    if esCuadrada(A):
        U = np.copy(A)
        f, c = U.shape
        i = 0
        while i < f:
            U[i:,:i+1] = 0
            i+= 1
        return U


def triangInf(A):
    if esCuadrada(A):
            L = np.copy(A)
            f, c = L.shape
            i = 0
            while i < f:
                L[:i+1,i:] = 0
                i+= 1
            return L


def diagonal(A):
    if esCuadrada(A):
        D = np.copy(A)
        f, c = D.shape
        i = 0
        while i < f:
            D[:i,i:] = 0
            D[i:,:i] = 0
            i+=1
        return D


def traza(A):
    if esCuadrada(A):
        tr = 0
        f, c = A.shape
        i = 0
        while i < f:
            tr += A[i,i]
            i += 1
        return tr


def traspuesta(A):
    f,c = A.shape
    T = np.zeros((c,f))
    i = 0
    while i < f:
        T[i:,:] = A[:,i]
        i+=1
    return T


def esSimetrica(A):
    if esCuadrada(A):
        f,c = A.shape
        T = traspuesta(A)
        simetrica = True
        i = 0
        while simetrica and i < f:
            for j in range(c):
                simetrica = A[i,j] == T[i,j]
            i += 1
        
        return simetrica
                   

def calcularAx(A,x):
    f,c = A.shape
    b = np.zeros(f)
    for i in range(f):
        n = 0
        for j in range(c):
            n += A[i][j]*x[j]
        b[i] = n
    return b


def intercambiarFilas(A,i,j):
    filJ = A[j,:].copy()
    A[j,:] = A[i,:]
    A[i,:] = filJ


def sumar_fila_multiplo(A, i, j, s):
    multJ = A[j,:]*s
    A[i,:] += multJ


def esDiagonalmenteDominante(A):
    f,c = A.shape
    diagdom = True
    i=0
    while diagdom and i < f:
        sumaVal = 0
        for j in range(c):
            if j != i:
                sumaVal += abs(A[i][j])
        diagdom = sumaVal < abs(A[i][i])
        i += 1
            
    return diagdom


def matrizCirculante(v):
    A = np.zeros((len(v),len(v)))
    for i in range(len(v)):
        A[i,i:] = v[:len(v)-i]
        A[i,:i] = v[len(v)-i:]
    return A


def matrizVandermonde(v):
    A = np.zeros((len(v),len(v)))
    for i in range(len(v)):
        A[i,:] = v**i
    return A


'''Desarrollar una funcion numeroAureo(n) que estime el numero
aureo ϕ como Fk+1/Fk, siendo Fk el k-esimo numero de la sucesion de Fibonacci.
Para esto, formulen la sucesion de Fibonacci
Fk+1 = Fk +Fk-1
de forma matricial, usando la semilla F0 = 0,F1 = 1. Grafique el valor aproxi
mado de ϕ en funcion del numero de pasos de la sucesion considerado.'''
def numeroAureo(n):
    x0 = 0
    x1 = 1
    return


def matrizFibonacci(n):
    if n == 0:
        A = [[0]]
        return A
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i+j <= 1:
                A[i,j] = i + j
            else:
                A[i,j] = A[i,j-1] + A[i,j-2]
    return A

print(matrizFibonacci(3))

'''
print("es cuadrada")
A = np.array([[1,1],[2,2],[1,1]])
print(esCuadrada(A))

print("triangular superior")
B = np.array([[1,1,1],[2,2,2],[3,3,3]]) 
print(B)
print(triangSup(B))

print("triangular inferior")
B = np.array([[1,1,1],[2,2,2],[3,3,3]]) 
print(B)
print(triangInf(B))

print("diagonal")
B = np.array([[1,1,1],[2,2,2],[3,3,3]]) 
print(B)
print(diagonal(B))

print("traza")
B = np.array([[1,1,1],[2,2,2],[3,3,3]]) 
print(B)
print(traza(B))

print("traspuesta")
B = np.array([[1,1,1],[2,2,2],[3,3,3]]) 
print(B)
print(traspuesta(B))

print("es simetrica")
B = np.array([[1,1,1],[2,2,2],[3,3,3]]) 
print(B)
print(esSimetrica(B))

print("es simetrica")
B = np.array([[1,1,1],[1,1,1],[1,1,1]]) 
print(B)
print(esSimetrica(B))

print("calcular Ax")
B = np.array([[1,1,1],[1,1,1],[1,1,1]]) 
r = np.array([3,3,4])
print(B)
print(r)
print(calcularAx(B,r))

print("intercambiar filas")
B = np.array([[1,1,1],[0,0,0],[1,2,3]])
print(B)
intercambiarFilas(B,1,2)
print(B) 

print("sumar multiplo fila")
B = np.array([[1,1,1],[1,1,0],[1,2,3]])
print(B)
sumar_fila_multiplo(B,1,2,2)
print(B)

print("diagonalmente dom")
B = np.array([[-3,1,1],[1,2,0],[1,2,4]])
print(B)
print(esDiagonalmenteDominante(B))
print(B) 

print("matriz circulante")
b = np.array([1,2,3])
print(matrizCirculante(b))

print("matriz vandermonde")
b = np.array([1,2,3])
print(matrizVandermonde(b))'''