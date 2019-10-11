#!/usr/bin/python
import matplotlib.pyplot as plt

def panel_plot(i,j):
    fig = plt.figure()
    n = 1
    for I in range(i):
        for J in range(j):
            fig.add_subplot(i,j,n)
            n+=1
    fig.add_subplot(i,j,1)
    return(fig)


