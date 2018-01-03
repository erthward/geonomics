#!/usr/bin/python

'''
Plotting both the fitness function used by Yeaman & Whitlock 2011 and Yeaman 2015 (as modified Spichtig &
Kawecki 2004), and the fitness function I'm using.

My fitness function is effectively the same, but will behave correctly in an environmental gradient because 
no division by zeros (i.e. my function's implicit normalizing denominator is 1, which is the max differerence between an
individual's phenotype and it's cell's optimal phenotype (and my environmental optima range from 0 to 1), whereas Yeaman & Whitlock's denominator is 2*theta,
which behaves fine in their divergent patch-optima of -1 and +1, but which would cause division by zero in my
spatial gradients, where a cell could have an optimum of 0)
'''

import numpy as np
import matplotlib.pyplot as plt

def yw_w(s, theta, z, gamma):
    return(1-s*((theta-z)/(2*theta))**gamma)

def my_w(s, e, z, gamma):
    return(1-s*(np.abs(e-z))**gamma)

yw_e = np.arange(-2,2,0.01)

my_e = np.arange(-0.5,1.5,0.01)


s = 0.1

gamma = 1


fig = plt.figure()
fig.suptitle('Fitness [w_z(e)] as a function of environmental optimal phenotype [e] and phenotype [z] for min, mean, and max phenotypes [red, black, and blue]')

ax = fig.add_subplot(2,2,1)
yw_ws = {-1 : [yw_w(s, theta, -1, gamma) for theta in yw_e],
         +1 : [yw_w(s, theta, +1, gamma) for theta in yw_e],
         0  : [yw_w(s, theta, 0, gamma) for theta in yw_e]
         }
ax.set_title('Yeaman, gamma = 1 (i.e. linear)')
plt.plot(yw_e, yw_ws[-1], color = 'red')
plt.plot(yw_e, yw_ws[+1], color = 'blue')
plt.plot(yw_e, yw_ws[0], color = 'black')
plt.plot([-1,-1], [.8,1], color = 'red', linestyle = ':')
plt.plot([1,1], [.8,1], color = 'blue', linestyle = ':')
plt.ylim([0.8,1])
plt.ylabel('w_z(e)')


ax = fig.add_subplot(2,2,2)
ax.set_title('Mine, gamma = 1 (i.e. linear)')
my_ws = {0 : [my_w(s, e, 0, gamma) for e in my_e],
         1 : [my_w(s, e, 1, gamma) for e in my_e],
         0.5  : [my_w(s, e, 0.5, gamma) for e in my_e]
         }
plt.plot(my_e, my_ws[0], color = 'red')
plt.plot(my_e, my_ws[1], color = 'blue')
plt.plot(my_e, my_ws[0.5], color = 'black')
plt.plot([0,0], [.8,1], color = 'red', linestyle = ':')
plt.plot([1,1], [.8,1], color = 'blue', linestyle = ':')
plt.ylim([0.8,1])
plt.ylabel('w_z(e)')


gamma = 2

ax = fig.add_subplot(2,2,3)
yw_ws = {-1 : [yw_w(s, theta, -1, gamma) for theta in yw_e],
         +1 : [yw_w(s, theta, +1, gamma) for theta in yw_e],
         0  : [yw_w(s, theta, 0, gamma) for theta in yw_e]
         }
ax.set_title('Yeaman, gamma = 2 (i.e. Gaussian)')
plt.plot(yw_e, yw_ws[-1], color = 'red')
plt.plot(yw_e, yw_ws[+1], color = 'blue')
plt.plot(yw_e, yw_ws[0], color = 'black')
plt.plot([-1,-1], [.8,1], color = 'red', linestyle = ':')
plt.plot([1,1], [.8,1], color = 'blue', linestyle = ':')
plt.ylim([0.8,1])
plt.ylabel('w_z(e)')
plt.xlabel('e')


ax = fig.add_subplot(2,2,4)
ax.set_title('Mine, gamma = 2 (i.e. Gaussian)')
my_ws = {0 : [my_w(s, e, 0, gamma) for e in my_e],
         1 : [my_w(s, e, 1, gamma) for e in my_e],
         0.5  : [my_w(s, e, 0.5, gamma) for e in my_e]
         }
plt.plot(my_e, my_ws[0], color = 'red')
plt.plot(my_e, my_ws[1], color = 'blue')
plt.plot(my_e, my_ws[0.5], color = 'black')
plt.plot([0,0], [.8,1], color = 'red', linestyle = ':')
plt.plot([1,1], [.8,1], color = 'blue', linestyle = ':')
plt.ylim([0.8,1])
plt.ylabel('w_z(e)')
plt.xlabel('e')

plt.show()
