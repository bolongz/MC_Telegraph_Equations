#This is the code for the example at 4.2.2. Figure 8

import math
import random as rd

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import cm
import time
from matplotlib.legend_handler import HandlerLine2D
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)


from mpmath import *

a = 1.5
c = 1.

NUM = 50

#chagne t_max = 0.3
t_max = 0.5
r_max = 1.

dt = 0.001
dx = 0.002

K = int(t_max / dt)
N = int(r_max / dx)

S = int(0.1/dx)

values = np.empty((K+1, N+1))

for i in range(N+1):
    rr = i * dx
    values[0, i] = r_max -  (i * dx ) * (i * dx) * (i * dx ) * (i * dx)
    values[1, i] = values[0,i] #1. -  (i * dx ) * (i * dx)
for i in range(K+1):
    values[i, N] = 0.


_j = 1

dtt = dt * dt
dxx = dx * dx
while(_j != K):

    for i in range(N):

        if i == 0:
            values[_j+1, i] = ((4 * dtt/dxx) * values[_j, i+1] + (2 - 4*
                dtt/dxx) * values[_j, i] + (3./2. *dt -1.) * values[_j-1, i]) / (1 + 3./2. *dt)
        else:
            ri = i * dx
            ii = i
            values[_j +1, ii] = (( 2 * dtt / dxx * ri + dtt / dx ) * values[_j,
                ii+1] + (2 *dtt/dxx * ri - dtt/dx) * values[_j, ii-1] - (
                        dtt/dxx * ri * 4) * values[_j, ii] - (2 *ri - 3 *dt *
                            ri) * values[_j - 1, ii] + 4 * ri * values[_j, ii]) / (2 *ri + 3 *dt * ri)

    _j = _j + 1

real_res = []
r_table = []
print 0, values[K][0]
print S * dx, values[K][S]
print 2 * S * dx, values[K][2*S]
print 3 * S * dx, values[K][3*S]
print 4 * S * dx, values[K][4*S]
print 5 * S * dx, values[K][5*S]
print 6 * S * dx, values[K][6*S]
print 7 * S * dx, values[K][7*S]
print 8 * S * dx, values[K][8*S]
print 9 * S * dx, values[K][9*S]

r_table.append(0)
real_res.append(values[K][0])
real_res.append(values[K][S])
r_table.append(S * dx)
real_res.append(values[K][2*S])
r_table.append(2 * S * dx)
real_res.append(values[K][3*S])
r_table.append(3 * S * dx)
real_res.append(values[K][4*S])
r_table.append(4 * S * dx)
real_res.append(values[K][5*S])
r_table.append(5 * S * dx)
real_res.append(values[K][6*S])
r_table.append(6 * S * dx)
real_res.append(values[K][7*S])
r_table.append(7 * S * dx)
real_res.append(values[K][8*S])
r_table.append(8 * S * dx)
real_res.append(values[K][9*S])
real_res.append(0)
r_table.append(9 * S * dx)
r_table.append(1)
r_res = np.asarray(real_res)


r_tab = np.asarray(r_table)

# Takes a list of times as input
def randtimes(l):
    T = 0
    #n = len(l)
    n = 1
    #result = -1. * np.ones(n)
    result = -1
    u = -1
    i = 0
    r = 0
    while (i < n ) and (result  < 0):
        s  = np.random.exponential(1. /a)
        T += s
        u = -u
        r += u * s
        while (i < 1) and T > l: # Need to use while loop because we could jump multiple epochs in one step
            result = c*(r - u * (T - l))
            i += 1
    return result



#def psi_2D(X):
#    return math.cos((X[0] + math.sqrt(4) * X[1])/ math.sqrt(5))


j0 = lambda x: besselj(0,x)
j1 = lambda x: besselj(1,x)
j2 = lambda x: besselj(2,x)
j3 = lambda x: besselj(3,x)



def one_item(bz, r, t):
    _j0 = j0(bz * r)
    _j1 = j1(bz)
    _j2 = j2(bz)
    _j3 = j3(bz)
    return (8.*(bz * _j2 - 2. * _j3 )/ (math.pow(bz, 3) * math.pow(_j1, 2))) * _j0 *math.cos(bz *t)



def wave_sol(r,t):
    ans = 0.0
    for i in range(NUM):
        _i = i + 1
        bz = besseljzero(0,_i)
       # _x = one_item(bz, r,t)
       # print i, _x
        ans = ans + one_item(bz, r,t)
    return ans



def one_eval(r, t):
    newtimes = randtimes(t)
    #print newtimes
    #n = len(times)

    n = 1
    result = np.empty(n)
    #for i in range(1):
    result[0] = wave_sol(r, newtimes)
    return result

# retursn a gigantic array with all the calculations inside. For memory efficieny, not necessary
def mc_method(r, t):
    '''
    n = 1
    result = np.empty((1, N_MC))
    for i in range(N_MC):
        result[0, i] = one_eval(r, t)
    return result
    '''
    n = 1
    #result = np.empty((n, N_MC))
    stderr = np.zeros(n)
    sum1 = np.zeros(n)
    sum2 = np.zeros(n)
    _t = 1000.
    i = 1;
    err = 0.005
    while( _t > err or i < 100):
        res = one_eval(r, t)
        sum1 = sum1 + res
        sum2 = sum2 + np.square(res)
        if( i > 2):
            stderr = np.abs(sum2-  (np.square(sum1) / (i))) / ( i - 1.)
        _t = math.sqrt((stderr/(i)) / (np.abs(sum1/i)))
        #if(i % 200 == 0):
        #    print i, _t
        i = i + 1

    return sum1 / (i - 1)



#deviations = np.empty(len(t_tab))
mp.dps = 25
mp.pretty = True
r_tabb = [0.1, 0.2, 0.3, 0.4, .5,.6,.7,.8,.9]
#r_tabb = [0.1]
e_res = []
for k, r in enumerate(r_tabb):


    plt.clf()
    plt.close()

    res = mc_method(r, t_max)
    #mean = np.sum(res, axis=1) / N_MC
    print r, res
    e_res.append(res)
    #nx = 20 # number of positive space points
    #ny = 20
    #hx = 0.5 #step between two space points
    #hy = 0.5
    #hx = _a / (2. *nx)
    #hy = _b / ( 2. * ny)
    #X = (np.arange(2 * nx + 1 ) - nx) * hx

    #Y = (np.arange(2 * ny + 1 ) - ny) * hy


print ""
#print t_tab, deviations
res_1 = mc_method(1,t_max)
e_res.append(res_1)

m_res = np.asarray(e_res)
#r_tabb = np.asarray(r_tabb)

plt.clf()
line1 =  plt.plot(r_tab, r_res,marker = 'o', color = 'r',markersize = 10,  linestyle = 'dashed', label =
            'FDM' )
line2 =  plt.plot(r_tab, m_res,marker = 'd', color = 'b',markersize = 10,  linestyle = 'dashed', label = 'Presented Algorithm')

plt.ylim((0, 1))
plt.xlim((0, 1))
plt.ylabel('u(t,r, $\\theta$)', size = 20)
plt.xlabel('x', size = 20)
plt.title('t = '+str(t_max),fontweight='bold',  size = 20)

plt.legend(scatterpoints=1,loc='lower left', numpoints = 1,
           prop={'size':15})


plt.show()
plt.savefig("img/t="+str(t_max)+".pdf", bbox_inches='tight')

plt.clf()
