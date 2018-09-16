#this the code for example at 4.1.2 Figure 6. To generate the data in table 6, we need to set the stoch_points to 0


import math
import random as rd

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from matplotlib.legend_handler import HandlerLine2D

import time

import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)


a = 1.
c = 1.

t_max = 0.2

alpha= -1.
beta = 1.

x_max = beta + c * t_max + 0.1
x_min = alpha - c * t_max - 0.1

#change the value for dx
#dx = 0.002
dx = 0.005
dt = 0.001

measuring_time = time.time()

K = int(t_max / dt)
N = int((x_max - x_min) / dx)

print "time division"
print K
print "x division"
print N

space_points = np.linspace(x_min, x_max, N)

# Initial conditions

const = 5.

def phi(x):
    if (alpha < const * x <  beta):
        return const *  math.exp(-1. / (1. -(const*x)**2))
    else:
        return 0.

def Psi(x):
    if (0. < const*x <  beta):
        return const * math.exp(-1. / (const*x*(1. - const*x)))
    else:
        return 0.

def phi_stable(x):
    if (alpha <= x <= beta):
        return phi(x)
    else:
        return 0.

def psi(x):
    if (0. < const*x <  beta):
        return const * (1.- 2.*const*x)  / (const*x*(1-const*x))**2 * math.exp(-1. / (const*x*(1. - const*x)))
    else:
        return 0.

X0 = np.zeros((2*N, 1))
for i in range(N):
    xi = space_points[i]
    X0[i] = phi_stable(xi)
    X0[N + i] = psi(xi)

# Definition of matrices used for the scheme

A = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if (j == i+1)or(j == i-1):
            A[i, j] = -1
        if (j == i):
            A[i, j] = 2
A = 1. / dx**2 * A

B = np.zeros((2*N, 2*N))
B[N:, :N] = A
B[:N, N:] = - np.eye(N)
B[N:, N:] = 2 * a * np.eye(N)
B = np.asmatrix(B)
B2 = np.dot(B, B)
I = np.asmatrix(np.eye(2*N))

mat = np.dot(12 * I  - 6*dt*B + dt**2 * B2, np.linalg.linalg.inv(12*I  + 6*dt*B + dt**2 * B2))

# Iterations

results = np.empty(K, dtype=object)
results[0] = X0


for k in range(1,K):
    X = results[k-1]
    results[k] = np.dot(mat, X)

measuring_time = time.time() - measuring_time

print "Time taken first method : "+str(measuring_time)

line1, = plt.plot(space_points, results[K-1][:N], color = 'g', label = "FDM[19]", lw=2)


########################### Stochastic solver : ##############################

N_MC =10000

# Solution to the wave equation with same conditions
def sol(t, x):
    return (phi(x + t) + phi(x -  t)) / 2. +  (Psi(x +  t) - Psi(x - t)) / (2.*c)

def randtimes(l):
    T = 0
    #n = len(l)
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
        while (i < n) and T > l:
            result[i] = c * (r - u * (T - l))
            i += 1
    return result

# Evaluates the solution for one random generation of the time !
def one_eval(times, points):
    newtimes = randtimes(times)
    n = len(points)
    #n = 1
    result = np.empty(n)
    for i in range(n):
        #result[i] = sol(newtimes[i], points[i])
        result[i] = sol(newtimes, points[i])
    return result

def F_(times, points):
    #n = 1
    n = len(points)
    result = np.empty((n, N_MC))
    for i in range(N_MC):
        result[:, i] = one_eval(times, points)
    return result

def treat_result(times, points, r):
    mean = np.sum(r, axis=1) / N_MC
    int_conf = 1.96 * np.vectorize(math.sqrt)(np.var(r, axis=1) / N_MC)
    return mean, int_conf

nx = N / 2 + 1   # number of positive space points
h = dx #step between two space points

#stoch_times = t_max * np.ones(N)
stoch_times = t_max
stoch_points = space_points
#stoch_points = np.array([0.0])
measuring_time = time.time()

r = F_(stoch_times, stoch_points)
measuring_time = time.time() - measuring_time
mean, int_conf = treat_result(stoch_times, stoch_points, r)


plt.ylim((-0.1, 1.4))
line2, = plt.plot(stoch_points, mean, color='b', label = "Presented algorithm", lw=2)
line3, = plt.plot(stoch_points, mean + int_conf, '--', color='r', label = "95% confidence interval", lw=2)
plt.plot(stoch_points, mean - int_conf, '--', color='r',lw=2)
plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)}, prop={'size':15})

plt.ylabel('u(t,x)', size = 20)
plt.xlabel('x', size =20)
plt.title('t = '+str(t_max), fontweight='bold', size = 20)

plt.savefig("labimg/t="+str(t_max)+".pdf", bbox_inches='tight', size = 20)
plt.clf()

print "2nd method time taken : "+str(measuring_time)

deterministic_r = np.asarray(results[K-1][:N].transpose())[0]

diff=abs(mean - deterministic_r[N/2]) / deterministic_r[N/2]

#print "error ", diff

u = np.abs((mean - deterministic_r ) / deterministic_r)
u = u * (u < 1)

## Comparison
max_diff = np.max(u)
argmax_diff = np.argmax(u)

print "maximum_difference : "+ str(max_diff)

plt.plot(stoch_points, u)

plt.ylabel('Relative difference')
plt.xlabel('x')
plt.title('t = '+str(t_max))

plt.savefig("labimg/diff_for_t="+str(t_max)+".pdf", bbox_inches='tight')

print "Max half-width of the confidence interval : " + str(np.max(int_conf))

u = u * (u > 0.1)
print "Proportion of values of x with D > 10 % : " + str(float(sum((u>0))) / N)
