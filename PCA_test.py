#!/usr/bin/python
# PCA_test.py

# import geonomics
import geonomics as gnx

# other imports
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os

# set some plotting params
img_dir = ('/home/drew/Desktop/stuff/berk/research/sim/methods_paper/'
                               'img/final/')
ax_fontdict = {'fontsize': 20,
                                                 'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 15,
                                                'name': 'Bitstream Vera Sans'}


# function for running and plotting genetic PCA
def plot_genetic_PCA(species):
    # get array of resulting genomic data (i.e. 'speciome'), genotypes meaned
    # by individual
    speciome = np.mean(np.stack([i.g for i in mod.comm[0].values()]), axis=2)
    # run PCA on speciome
    pca = PCA(n_components=3)
    PCs = pca.fit_transform(speciome)
    # normalize the PC results
    norm_PCs = (PCs - np.min(PCs,
                             axis=0)) / (np.max(PCs,
                                                axis=0) - np.min(PCs, axis=0))
    # use first 3 PCs to get normalized values for R, G, & B colors
    PC_colors = norm_PCs * 255
    # scatter all individuals on top of landscape, colored by the
    # RBG colors developed from the first 3 geonmic PCs
    xs = mod.comm[0]._get_x_coords()
    ys = mod.comm[0]._get_y_coords()
    mod.plot()
    plt.scatter([*xs], [*ys], c=PC_colors/255.0)


# make empty figure
fig = plt.figure()
# make model
mod = gnx.make_model('./geonomics/tests/validation/PCA/PCA_params.py')

# define number of timesteps
T = 1000

# burn model in
mod.walk(20000, 'burn')

# plot genetic PCA before genomic evolution begins
ax1 = fig.add_subplot(121)
#ax1.set_title('Before evolution')
plot_genetic_PCA(mod.comm[0])

# run model for T timesteps
mod.walk(T)

# plot genetic PCA after 1/4T timesteps
ax2 = fig.add_subplot(122)
#ax2.set_title('After %i timesteps' % T)
plot_genetic_PCA(mod.comm[0])

# add title
#plt.suptitle(('Neutral genomic evolution across complex landscape with '
#              '_MovementSurface,\n(for a ~%i-individual species with 100 '
#              'loci') % int(np.mean(mod.comm[0].Nt)))
plt.show()
plt.savefig(os.path.join(img_dir, 'PCA_before_after_plot.pdf'))
