#!/usr/bin/python

import demography

def plot_logistic_eqxn(R, N):

    K = np.arange(0,10*N, 0.01)

    plt.plot(K, [demography.logistic_eqxn(R, N, k) for k in K], color = 'black') 

    plt.plot([0,10*N],[0,0], linestyle = ':', color = 'black')
    plt.plot([0,10*N], [R]*2, linestyle = ':', color = 'red')
    plt.plot([N]*2, [-R-3,0], linestyle = ':', color = 'black')
    
    plt.ylim([-R-3,R+3])
    plt.xlabel('K')
    plt.ylabel('dN/dt')
