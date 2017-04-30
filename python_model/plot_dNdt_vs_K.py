#!/usr/bin/python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import demography as d

R = [0.1, 0.25, 0.5,1]

N = [0,1,2,3,5,10,15,20,50,100]

c = ['black', 'grey', 'red', 'magenta', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple']

K = np.arange(0,100,0.01)


fig = plt.figure()
fig.suptitle('dN/dt as a function of K for varying values of N (within plot) and r (across plots)')
subplot_num = 100*len(R)
for i, r in enumerate(R):
    s_n = subplot_num+11+i
    plt.subplot(s_n)
    plt.title('r = %0.4f' % r)
    if i == 2:
        plt.ylabel('log(abs(dN/dt))')
    if i == len(R) - 1:
        plt.xlabel('K')
    for ii, n in enumerate(N):
        dndt = [log(abs(d.logistic_eqxn(r,n,k))) for k in K]
        #dndt = [log(-val) if val < 0 else log(val) for val in dndt]
        plt.plot(K, dndt, label = str(n), color = c[ii])
        plt.ylim(-10,15)
    if i == len(R) - 1:
        plt.legend(loc = 4, title = 'N', prop = {'size':6})

