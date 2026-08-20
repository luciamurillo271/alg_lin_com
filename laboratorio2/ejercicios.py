print(format(0.1, '.20f')) #'0.10000000000000000555'

'''Ejercicio 1
print(x == y) da False. Utilizar format() para verificar si la respuesta es correcta '''
x = 0.1 + 0.1 + 0.1
y = 0.3

print(x == y) 
print(format(x))
print(format(y))

'''Ejercicio 2
Correr el siguiente programa en Python. Identificando la fuente
del error, proponer una forma de solucionar su mal funcionamiento'''
a = 1.0
while a != 0.1:
    print(a)
    a = a - 0.1
    a = float(format(a,'.1f'))
#print('fin')

'''Ejercicio 3
1. Comparen el resultado de hacer 0.3 + 0.25 con el de hacer 0.3 - 0.25. ¿En
ambos casos obtienen el resultado esperado? ¿Por qué?
2. Escriban el número 0.25 en base 2. ¿Cómo queda expresado en términos
de su mantisa y exponente?
3. Escriban el número 0.3 en base 2. ¿Qué dificultades aparecen al escribir
0.3 en binario? ¿Se puede escribir exactamente con una mantisa finita'''

x = 0.3 + 0.25
print(x)
y = 0.3 - 0.25
print(y)

x = 0.25
ceros = 0
while x != 1.0:
    x = x*2
    ceros += 1
#print(ceros)

print('0.25 = 0,01')
# (1.0)2 * 2 a la -2

x = 0.3
#0.3 no puede representarse exactamente como una fracción con denominador potencia de 2.
#no se puede escribir exactamente con una mantisa finita

'''Ejercicio 4
(No tan distintos). En este punto exploraremos expresiones que son aparentemente iguales.
1. ¿Cuánto da (√2)^2 - 2? Simbólicamente sabemos que el resultado es 0, pero ¿qué ocurre en python?
Importen la librería numpy (import numpy as np) para emplear la función
np.sqrt y calculen np.sqrt(2)**2-2
2. Para 100 valores de equiespaciados x en el intervalo de [0, 5 × 10^(-8)], eval-
uar las siguientes 2 expresiones que son matematicamente equivalentes
(pruebenlo) y graficarlas usando matplotlib.pyplot. En base al gráfico
obtenido identificar la opción que mejor resiste la pérdida de valores sig-
nificativos.'''

import numpy as np

#from laboratorio1.librerias import traspuesta, esCuadrada,esSimetrica
a = np.sqrt(2)**2 -2
print(a)


fin = 5*(10**(-8))
vector = np.linspace(0,fin,100)

form_a = []
form_b = []
for elemento in vector:
    elem_a = np.sqrt(2*elemento**2 + 1) -1 
    form_a.append(elem_a)

    elem_b = (2*elemento**2)/(np.sqrt(2*elemento**2 + 1) +1 )
    form_b.append(elem_b)

#print(form_a)
#print(form_b)

from matplotlib import pyplot as plt

#plt.plot(form_a)
#plt.plot(form_b)
#plt.show()

'''Ejercicio 5
(Acumulación del error). Calculen algebráicamente el límite cuando
n → ∞ de esta sucesión
x1 = √2
xn+1 = (xn · xn)/ √2
Implementen una rutina que calcule el valor de xi para i = 1,...,100 y grafiquen
sus valores. ¿En qué punto se desestabiliza la sucesión? Tip: pueden guardar
los elementos de la sucesión en una lista l=[] usando l.append(xi) dentro
de un loop for. Luego importar matplotlib.pyplot as plt y graficarlos con
plt.plot(l).'''
#el limite de la sucesion es raiz de 2 pues todos los elementos de la sucesion son iguales a raiz de 2
'''
l = []
l.append(np.sqrt(2))
for i in range(1, 100):
    l.append((l[i-1] * l[i-1]) / np.sqrt(2))
plt.plot(l)

print(l)
plt.show()'''

'''Ejercicio 6
'''
import numpy as np

n = 7
s = np.float32(0)
for i in range(1,10**n+1):
    s = s + np.float32(1/i)
#print('suma = ', s)

s = np.float32(0)
for i in range(1,5*10**n+1):
    s = s + np.float32(1/i)
#print('suma = ', s)

'''Para pensar:
 ¿Cuánto vale 1/i en precisión simple cuando i = 2 · 107
?
 Si escribimos 1/107 usando el mismo exponente que el necesario para representar a P5·106
i=1 1/i, ¿a cuánto equivale 1/i?
 ¿Por qué la siguiente modificación cambia el resultado?
'''
s = np.float32(0)
for i in range(2*10**n,0,-1):
    s = s + np.float32(1/i)
#print('suma = ', s)

'''Ejercicio 7
(Arrastre de error: Descomposición LU). Desarrollar una
función matricesIguales(A, B) que devuelve True si ambas matrices son iguales.
Verificar que la función desarrollada matricesIguales(A, L@U) devuelva
True.'''
def matricesIguales(A, B):
    iguales = True
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        iguales = False
    else:
        iguales =A.all() == B.all()
        #for i in range(len(A)):
    return iguales

A = np.array([[4,2,1],[2,7,9],[0,5,22/3]])
L = np.array([[1,0,0],[0.5,1,0],[0,5/6,1]])
U = np.array([[4,2,1],[0,6,8.5],[0,0,0.25]])
#print(matricesIguales(A,L@U))

'''Ejercicio 8
(Arrastre de error: función esSimetrica). Verificar que devuelve la función esSimetrica() 
desarrollada para el laboratorio pasado en los
siguientes casos. Para una A = np.array(np.random.rand(4,4))
esSimetrica(A.T@A)
esSimetrica(A.T@((A*0.25)/0.25))
esSimetrica(A.T@((A*0.2)/0.2))'''

def esCuadrada(A):
    return A.shape[0] == A.shape[1]

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

A = np.array(np.random.rand(4,4))
print(esSimetrica(A.T@A))
print(esSimetrica(A.T@((A*0.25)/0.25)))
print(esSimetrica(A.T@((A*0.2)/0.2)))
