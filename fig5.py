# This is an code for the example at 4.1.1 to generate Figure 5


import math
import random as rd

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.legend_handler import HandlerLine2D


import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)
import time

a = 1.5
c = 1.
N_MC = 200001


ksi = -a + math.sqrt(a**2 - c**2)
ksibar = -a - math.sqrt(a**2 - c**2)

A = 1.
B = - 1.

def phi(x):
    return 0.0
# Solution to the wave equation with same conditions

def opsi(x):
    return (ksi - ksibar) * math.cos(x)

def Psi(x):
    return (ksi - ksibar) * math.sin(x)

def sol(t, x):
    return (phi(x + t) + phi(x - t)) / 2. +  (Psi(x + t) - Psi(x - t)) / (2.*c)

# Solution of the telegrapher's equation / For comparing purposes.
def real_sol(t, x):
    return (math.exp(ksi * t) - math.exp(ksibar * t )) * math.cos(x)
# argument is an array of ordered times
def randtimes(l):
    T = 0
    n = 1
    result = -1. * np.ones(n)
    u = -1
    i = 0
    r = 0
    while (i < n ) and (result[i] < 0):
        s  = np.random.exponential(1./a)
        T += s
        u = -u
        r += u * s
        while (i < n) and T > l: # Need to use while loop because we could jump multiple epochs in one step
            result[i] = c*(r - u * (T - l))
            i += 1
    return result

def route(l):
    n = 1
    T = 0
    u = -1
    r = 0
    i = 0
    result = -1. * np.ones(n)
    length = 0.0 * np.ones(n)
    k = 0
    while (i < n)  and (result[i] < 0):
        s  = np.random.exponential(1./a)
        T += s
        u = -u
        r += u * s
        k = k +1
        while( i < n) and T > l:
            result[i] = c * (r - u * (T - l))
            length[i] = k
            i+= 1
    return result, length
def jcp_op(t, x, l):
    f = 0.0
    b = 0.0

    '''
    if l % 2 == 1:

        return (phi(x + c * t) + phi(x - c * t)) / 2. +  (Psi(x + c * t) - Psi(x - c * t)) / (2.*c)
    else:

        return (phi(x + c * t) + phi(x - c * t)) / 2. -  (Psi(x + c * t) - Psi(x - c * t)) / (2.*c)
    '''
    if l % 2 == 1:
        f = phi(x - t) - 1. / c * Psi(x - t)
        b = phi(x + t) + 1. / c * Psi(x + t)
    else:

        f = phi(x + t) - 1. / c * Psi(x + t)
        b = phi(x - t) + 1. / c * Psi(x - t)

    return (f + b)/2

# Evaluates the solution for one random generation of the time !
def one_eval(times, points):
    newtimes = randtimes(times)
    n = len(points)
    result = np.empty(n)
    for i in range(n):
        result[i] = sol(newtimes, points[i])
    return result

def mc_method(times, points, err):
    n = len(points)
    stderr = np.zeros(n)
    sum1 = np.zeros(n)
    sum2 = np.zeros(n)
    t = 1000.
    i = 1;
    while( t > err or i < 1001):
        res = one_eval(times, points)
        sum1 = sum1 + res
        sum2 = sum2 + np.square(res)
        if(i % 2000 == 0):
            stderr = np.abs(sum2-  (np.square(sum1) / (i))) / ( i - 1.)
            t = max(np.vectorize(math.sqrt)(stderr/(i)) / (np.abs(sum1/i)))
            print i, t
        i = i + 1

    return sum1 / (i - 1), i - 1 #, int_conf

def jcp_opsi(times, points):
    n = len(points)
    #n = 1
    result = np.empty(n)
    newtimes, length = route(times)
    for i in range(n):
        result[i]= jcp_op(newtimes, points[i], length)
    return result

def jcp_sol(times, points, err):
    n = len(points)
    stderr = np.zeros(n)
    sum1 = np.zeros(n)
    sum2 = np.zeros(n)
    t = 1000.
    i = 1;
    while( t > err or i < 1001):
        res = jcp_opsi(times, points)
        sum1 = sum1 + res
        sum2 = sum2 + np.square(res)

        if(i % 2000 == 0):
            stderr = np.abs(sum2-  (np.square(sum1) / (i))) / ( i - 1.)
            t = max(np.vectorize(math.sqrt)(stderr/(i)) / (np.abs(sum1/i)))
            print i, t
        i = i + 1


    return sum1/ (i - 1), i - 1 #, int_conf


tab = [0.05, 0.5, 1., 2.]
max_deviations =  0.0
max_deviations2  = 0.0
err = 0.01
for i, t in enumerate(tab):
    print "****************"
    plt.clf()
    nx = 12 # number of positive space points

    h = 0.5 #step between two space points

    times = t  #* np.ones(2 * nx + 1)
    points = (np.arange(2 * nx + 1 ) - nx) * h
    start_time = time.time()
    mean, steps = mc_method(times, points,err)

    kacs_time = time.time() - start_time
    line1 =  plt.scatter(points, mean,s = 40, marker = 'o', color = 'r', label =
            'Presented algorithm' )
    start_time2 = time.time()
    mean1, steps2 = jcp_sol(times, points,err)

    jcp_time = time.time() - start_time2

    line2 =  plt.scatter(points, mean1,s = 40, marker = '^', color = 'b', label = 'NMC' )

    plt.ylim((-2.5, 2.5))
    plt.xlim((-6, 6))
    plt.ylabel('u(t,x)', size = 20)
    plt.xlabel('x', size = 20)
    plt.title('t = '+str(t),fontweight='bold',  size = 20)

    plt.legend(scatterpoints=1,loc='upper right', numpoints = 1,
            prop={'size':15})


    nx = 100
    h = 0.06
    points = (np.arange(2 * nx + 1 ) - nx) * h
    real_points = np.vectorize(lambda x: real_sol(t, x))(points)
    line3, = plt.plot(points, real_points,  color = 'g', label='Accurate solution',
                      linewidth = 2)
    plt.legend(handler_map={line3: HandlerLine2D(numpoints=2)}, prop =
            {'size':15})
    plt.savefig("img/t="+str(t)+".pdf", bbox_inches='tight')

    plt.clf()


    print t
    print max_deviations
    print max_deviations2
    print "kac based time"
    print kacs_time
    print "kac based steps"
    print steps
    print "jcp based time"
    print jcp_time
    print "jcp based steps"
    print steps2
