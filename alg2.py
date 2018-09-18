#This is the code for example at 4.1.1 to show the efficiency for algorithm 2.


import math
import random as rd

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import time

a = 1.5
c = 1.

N_MC = 20000


ksi = -a + math.sqrt(a**2 - c**2)
ksibar = -a - math.sqrt(a**2 - c**2)

Ksi = ksi - ksibar
cc = 2 * c
A = 1.
B = - 1.

def phi(x):
    return 0.0
def Psi(x):
    return Ksi * math.sin(x)
def sol(t, x):
    return  (Psi(x +  t) - Psi(x -  t)) / (cc)


# calculating the multi-points
#NT = 11
#h = 5. / NT
#tab = (np.arange(NT) +1) * h

#tab = [5. , 10.]
#tab = [2.5, 5.]
#tab = [1.,2.,3.,4.,5.]
#tab = [1.25, 2.5, 3.75, 5.]
#tab = [0.85, 1.7, 2.55,3.4,4.25, 5.]
#tab = [0.625, 1.25,1.875,2.5,3.125, 3.75, 4.375, 5.]

#tab = [0.5, 1., 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
tab=[1.]
#Calculting the points
#points = (np.arange(2 * nx + 1 ) - nx) * h
#points = np.array([0])
#points= np.array([-1.,0.,1.])
#points= np.array([-2.,-1.,0.,1.,2.])

points= np.array([-3., -2.,-1.,0.,1.,2.,3.])
#points= np.array([-4., -3., -2.,-1.,0.,1.,2.,3.,4.])
#points= np.array([-5., -4., -3., -2.,-1.,0.,1.,2.,3.,4.,5.])

n = len(points)
print points
nt = len(tab)

multiresult = np.zeros((n, nt)) #store the result
starttime = time.time()

for i in range(N_MC):
        T = 0
        u = -1
        j = 0;
        r = 0
        while j < nt:
            s  = np.random.exponential(1./a)
            T += s
            u = -u
            r += u * s
            while (j <  nt) and T > tab[j]:
                newt = r - u * (T - tab[j])
                for k in range (n):
                    multiresult[k, j] += (math.sin(points[k] +
                            newt) - math.sin(points[k] -  newt)) * Ksi / cc
                j += 1


multiresult = multiresult / N_MC

mtime = time.time() - starttime

print "Algorithm2 cost: "
print mtime

_multiresult = np.zeros((n, nt))

print tab
_starttime = time.time()
for k, t  in enumerate(tab):
    for i in range(N_MC):
        for j in range(n):
            T = 0
            u = -1
            r = 0
            jj = 0
            leng = 0
            while (jj < 1):
                s  = np.random.exponential(1./a)
                T += s
                u = -u
                r += u * s
                leng += 1
                if T > t:
                    newtimes = r - u * (T - t)
                    _multiresult[j, k] +=  (
                            math.sin(points[j] + newtimes) - math.sin(points[j]
                                -  newtimes)) * Ksi/cc
                    jj += 1

_multiresult = _multiresult / N_MC

_mtime= time.time() - _starttime

print "Algorithm1 Cost"
print _mtime

print "Speedup: ", _mtime / mtime
