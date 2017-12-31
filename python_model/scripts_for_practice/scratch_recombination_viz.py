#!/usr/bin/python


def vis_recomb(individ, chrom_num, genomic_arch):

    import numpy as np
    import gametogenesis
    import matplotlib.pyplot as plt

    g = list(individ.genome.genome[chrom_num][:,0])
    n_loc = len(g)
    chrom_plot_width = int(0.025*n_loc)
    cmap = 'cool'

    recomb_g = list(gametogenesis.gametogenerate(individ, genomic_arch)[chrom_num])

    fig = plt.figure()

    ax = fig.add_subplot(2,1,1)
    
    plt.imshow(np.vstack((np.vstack([np.array(g) for i in range(chrom_plot_width)]), np.vstack([np.array([float(loc == 0) for loc in g]) for i in range(chrom_plot_width)]))), cmap = cmap)

    plt.title('Non-recombinant')
    plt.yticks([],[])
    plt.xticks([],[])
    
    ax = fig.add_subplot(2,1,2)

    plt.imshow(np.vstack((np.vstack([np.array(recomb_g) for i in range(chrom_plot_width)]), np.vstack([np.array([float(loc == 0) for loc in recomb_g]) for i in range(chrom_plot_width)]))), cmap = cmap)

    plt.title('Recombinant')
    plt.yticks([],[])
    plt.xticks([],[])

    plt.plot(range(len(pop.genomic_arch.r[0][1:])), max(plt.ylim())-pop.genomic_arch.r[0][1:]*max(plt.ylim())*3, linewidth = 3) #exclude first value, which is always 0.5 to effectuate separation of homologous chromosomes

