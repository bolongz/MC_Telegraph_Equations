#this is the code for the example at 4.2.1 Figure 7.


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

# Takes a list of times as input
def randtimes(l):
    T = 0
    n = len(l)
    result = -1. * np.ones(n)
    u = -1
    i = 0
    r = 0
    while (i < n ) and (result[i] < 0):
        s  = np.random.exponential(1. /a)
        T += s
        u = -u
        r += u * s
        while (i < n) and T > l[i]:
            result[i] = r - u * (T - l[i])
            i += 1
    return result


ksi = -a + math.sqrt(a**2 - c**2)
ksibar = -a - math.sqrt(a**2 - c**2)
A = 0.3
B = -1.

def psi(X):
    return math.cos((X[0] + math.sqrt(4) * X[1])/ math.sqrt(5))

def wave_sol(t, X):
    alpha = A + B
    beta = (A * ksi  + B * ksibar) / c
    return (alpha * math.cos(c *t) + beta * math.sin(c*t)) * psi(X)

def real_sol(t, X):
   return (A * math.exp(ksi * t) + B * math.exp(ksibar * t) )* psi(X)

def one_eval(times, points):
    newtimes = randtimes(times)
    n = len(times)
    result = np.empty(n)
    for i in range(n):
        result[i] = wave_sol(newtimes[i], points[i])
    return result


def mc_method(times, points, err):
    n = len(times)
    #result = np.empty((n, N_MC))
    stderr = np.zeros(n)
    sum1 = np.zeros(n)
    sum2 = np.zeros(n)
    t = 1000.
    i = 1;
    while( i < N_MC):
    #while( t > err or i < 1000):
    #while( (t > err or i < 1000) and i < 10000):
        res = one_eval(times, points)
        sum1 = sum1 + res
        sum2 = sum2 + np.square(res)
        if( i > 2):
            stderr = np.abs(sum2-  (np.square(sum1) / (i))) / ( i - 1.)
        t = max(np.vectorize(math.sqrt)(stderr/(i)) / (np.abs(sum1/i)))
        if(i % 2000 == 0):
            print i, t
        i = i + 1

    return sum1 / (i - 1), i - 1

#t_tab = 0.1 * np.arange(1)
#t_tab = [1.]
t_tab = [0.1, 0.2, 1., 2.]
deviations = np.empty(len(t_tab))
err = 0.01
for k, t in enumerate(t_tab):

    #plt.clf()
    plt.close()

    nx = 20 # number of positive space points
    ny = 20
    #nx = 0
    #ny = 0

    hx = 0.2 #step between two space points
    hy = 0.2

    X = (np.arange(2 * nx + 1 ) - nx) * hx
    Y = (np.arange(2 * ny + 1 ) - ny) * hy
    #xx = -0.5
    #X = np.ones(1) * xx
    #Y = np.ones(1) * xx
    points = np.empty(((2 * nx + 1) , (2 * ny + 1), 2))

    for i, x in enumerate(X):
        for j, y in enumerate(Y):
            points[i, j, :] = np.array([x, y])

    points = points.reshape(((2 * nx + 1) * (2 * ny + 1), 2))

    #t = t + 1.0
    times = (t) * np.ones((2 * nx + 1) * (2 * ny + 1))

    start_time = time.time()


    mean, steps = mc_method(times, points,err)

    print t, mean

    Z = mean.reshape(2 * nx + 1, 2 * ny + 1)
    time_t = time.time()-start_time
    real_z = np.array([real_sol(t, point) for point in points])
    real_z = real_z.reshape((2*nx + 1, 2*ny +1))

    print t, real_z

    # For plotting purposes only

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')

    ax.set_ylabel('y')
    ax.set_xlabel('x')
    ax.set_zlabel('u(t, x)')

    ax.set_title('t = '+str(t), loc='left', fontweight = 'bold')

    X, Y = np.meshgrid(X, Y)
    surf = ax.plot_surface(X, Y, real_z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, alpha=0.5, antialiased=False)

    ax.set_zlim(-0.6, 0.6)

    fig.colorbar(surf, shrink=0.5, aspect=5)
    surf2 = ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1, linewidth=1, antialiased=False)
    plt.savefig("img/t = "+str(t)+".pdf", bbox_inches='tight')
    real_z = real_z.reshape((2 * nx + 1) * (2 * ny +1), 1)
    Z = Z.reshape((2 * nx + 1) * (2 * ny +1), 1)
    dev = np.max(np.abs(real_z - Z)/real_z)
