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
2. Escriban el número 0.25 en base 2. ¿Cómo queda expresado en términos
de su mantisa y exponente?
3. Escriban el número 0.3 en base 2. ¿Qué dificultades aparecen al escribir
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
1. ¿Cuánto da (√2)^2 - 2? Simbólicamente sabemos que el resultado es 0, pero ¿qué ocurre en python?
Importen la librería numpy (import numpy as np) para emplear la función
np.sqrt y calculen np.sqrt(2)**2-2
2. Para 100 valores de equiespaciados x en el intervalo de [0, 5 × 10^(-8)], eval-
uar las siguientes 2 expresiones que son matematicamente equivalentes
(pruebenlo) y graficarlas usando matplotlib.pyplot. En base al gráfico
obtenido identificar la opción que mejor resiste la pérdida de valores sig-
nificativos.'''

import numpy as np
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

print(form_a)
print(form_b)

from matplotlib import pyplot as plt

plt.plot(form_a)
plt.plot(form_b)
plt.show()
