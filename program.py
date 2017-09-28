import math
import sys
import numpy as np
import random


def polynom(num, p, n):
    myArray = [0 for i in range(n)]
    for i in range(n):
        myArray[n - 1 - i] = (num % p)
        num = (num // p)
    return myArray


def value(pol, x, p):
    n = len(pol)
    result = (x**n) % p 
    for i in range(n):
        result += (pol[n - 1 - i] * x**i) % p
    result %= p
    return result


def irreduciblePol(p, n):
    for num in range(p**n):
        pol = polynom(num, p, n)
        count = 0
        for x in range(p):
            if value(pol, x, p) > 0:
                count += 1
                if count == p:
                    result = [1 for j in range(n+1)]
                    for j in range(n):
                        result[j+1] = pol[j]
                    return result


def printPol(pol, n):
    print('����������� ���������:', end='   f(x) = ')
    for i in range(n+1):
        if (pol[i] > 0):
            if (i > 0 ):
                print(' + ', end='')
            if (pol[i] > 1 or i == n):
                print(pol[i], end='')
            if (i < n):
                print('x', end='')
            if (i < n - 1):
                print('^', n - i, sep='', end='')
    print()


def printElement(pol, n):
    first = True
    for i in range(n+1):
        if (pol[i] > 0):
            if not first:
                print(' + ', end='')
            if (pol[i] > 1 or i == n):
                print(int(pol[i]), end='')
                first = False
            if (i < n):
                print('a', end='')
                first = False
            if (i < n-1):
                print('^', n - i, sep='', end='')
                first = False


def generateElementsOfField(p, n):
    print('�������� ����:')
    print('   ', 0, end='')
    for num in range(p**n):
        pol = polynom(num, p, n)
        print('    ', end='')
        printElement(pol, len(pol)-1)
        print()



def polyToInt(pol, p):
    result = 0
    deg = (len(pol) - 1)
    for i in range(deg + 1):
        result += pol[i] * p**(deg - i)
    return result


def multiply(pol1, pol2, minim, p):
    deg1 = len(pol1)
    deg2 = len(pol2)
    result = [0 for j in range(deg1 + deg2 - 1)]
    for i in range(deg1):
        for j in range(deg2):
            result[deg1 + deg2 - 2 - i - j] += pol1[deg1 - 1 - i] * pol2[deg2 - 1 - j]
    for j in range(deg1 + deg2 - 1):
        result[j] %= p
    quotient, remainder = np.polydiv(result, minim)
    n = len(remainder)
    for i in range(n):
        remainder[i] %= p
    return remainder
    

def printMatrix(pol, p, n, A, E):
    first = True
    for i in range(n+1):
        if (pol[i] > 0):
            if not first:
                print(' + ', end='')
            if (pol[i] > 1):
                print(int(pol[i]), end='')
                first = False
            if (i == n):
                print('E', end='')
            if (i < n):
                print('A', end='')
                first = False
            if (i < n-1):
                print('^', n - i, sep='', end='')
                first = False
    mat = pol[n] * E
    for i in range(1, n+1):
        mat += pol[n - i] * A**(i)
    mat %= p
    print(' =\n', mat)   


def generateMatricesOfField(p, n, A, E, O):
    print('�������� ���� � ��������� �������������:')
    print('O =\n', O)
    for num in range(1, p**n):
        pol = polynom(num, p, n)
        printMatrix(pol, p, len(pol)-1, A, E)






print('������� ����� ������� �����')
p = int(input())
flag = True
if p > 1:
    d = math.sqrt(p)
    for i in range(2,math.floor(d) + 1):
        if (p % i) == 0:
            print('����� p �� �������� �������')
            flag = False
else:
    print('����� p �� �������� �������')
    flag = False


if not flag:
    sys.exit(0)


print('������� ����� ����������� �����')
n = int(input())
if n < 1:
    print('����� n �� �������� �����������')
    flag = False

if not flag:
    sys.exit(0)


print('p =', p)
print('n =', n)
print('q = p^n =', p**n)


minPoly = irreduciblePol(p, n)


printPol(minPoly, n)
print('\na - ������ f(x): f(a) = 0\n')



minimal = polyToInt(minPoly, p)


generateElementsOfField(p, n)



print()
print()


primitiveElements = {i for i in range(p, p**n)}


primElem = minPoly
primElemFound = False
while not primElemFound:
    num = (random.sample(primitiveElements, 1))[0]
    poly = polynom(num, p, n)
    prevPoly = poly
    for i in range(2, p**n):
        newPoly = multiply(poly, prevPoly, minPoly, p)
        newInt = polyToInt(newPoly, p)
        primitiveElements.discard(newInt)
        prevPoly = newPoly
        if (len(primitiveElements) == 1):
            primElem = poly
            primElemFound = True

print('\n������� ��������:')
print('\ti\tzeta^i\n')
print('\t', 1, ':\t', sep='', end='')
printElement(primElem, n-1)
print()
prevDeg = primElem
for i in range(2, p**n):
    nextDeg = multiply(primElem, prevDeg, minPoly, p)
    print('\t', i, ':\t', sep='', end='')
    printElement(nextDeg, len(nextDeg) - 1)
    print()
    prevDeg = nextDeg
print('\n')



print('��������� �������������:\n')

A = np.zeros((n, n), int)
for i in range(n):
    for j in range(n):
        if (i == j + 1):
            A[i][j] = 1
        if (j == n - 1):
            A[i][j] =  -minPoly[n - i] % p

print('�������������� ������� �:')
print('A =\n', A)

print('������� � ��������� �������:')
O = np.zeros((n, n), int)
print('O =\n', O)
E = np.zeros((n, n), int)
np.fill_diagonal(E, 1)
print('E =\n', E)


generateMatricesOfField(p, n, A, E, O)
    
