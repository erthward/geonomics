#!/usr/bin/python
# IBD_IBE_test.py

# flake8: noqa

# import geonomics
import geonomics as gnx

# other imports
from copy import deepcopy
# from itertools import chain
import numpy as np
from sklearn.decomposition import PCA
# import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
import os
import time
import statsmodels.api as sm


# set flag to indicate whether or not to use the barrier
# (this will typically be used, but can be turned off to run the model without
# a barrier and thus determine, comparatively, how much the barrier increases
# landscape resistance)
use_barr = True

# set flag to indicate whether or not I'm timing this demo
# (because if I am then I don't want to be calculating the mean z-e diff
# at each timestep)
not_timing = True


# set some plotting params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 12,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 15,
                'name': 'Bitstream Vera Sans'}
mark_size = 15


def map_genetic_PCA(species, land):
    """run genetic PCA and mapping individuals colored in RGB
    space by the first 3 PCs"""
    # get array of resulting genomic data (i.e. 'speciome'),
    # genotypes meaned by individual
    speciome = np.mean(np.stack([i.g for i in species.values()]), axis=2)
    # run PCA on speciome
    pca = PCA(n_components=3)
    PCs = pca.fit_transform(speciome)
    # normalize the PC results
    norm_PCs = (PCs - np.min(PCs,
                             axis=0)) / (np.max(PCs,
                                                axis=0) - np.min(PCs,
                                                                 axis=0))
    # use first 3 PCs to get normalized values for R, G, & B colors
    PC_colors = norm_PCs * 255
    # scatter all individuals on top of landscape, colored by the
    # RBG colors developed from the first 3 geonmic PCs
    xs = mod.comm[0]._get_x()
    ys = mod.comm[0]._get_y()
    # get environmental raster, with barrier masked out
    masked_env = deepcopy(mod.land[0].rast)
    masked_env[mod.land[1].rast == 0] = np.nan
    # create light colormap for plotting landscape
    # bot = plt.cm.get_cmap('Blues', 256)(np.linspace(0.4, 0.45, 2))[0]
    # top = plt.cm.get_cmap('Reds', 256)(np.linspace(0.4, 0.45, 2))[0]
    # cols = np.vstack((top, bot))
    # cmap = mpl.colors.ListedColormap(cols, name='OrangeBlue')
    cmap = plt.cm.coolwarm
    cmap.set_bad(color='#8C8C8C')
    # plot landscape
    # plt.imshow(masked_env, cmap=cmap, alpha=0.8)
    mod.plot(0, 0, color=PC_colors/255.0, mask_rast=masked_env,
             size=mark_size, edge_color='black')
    #plt.pcolormesh(land._x_cell_bds, land._y_cell_bds, masked_env, cmap=cmap)
    # scatter plot of individuals, colored by composite PC score
    #plt.scatter(xs, ys, c=PC_colors/255.0, s=mark_size, edgecolors='black')
    # fix x and y limits
    #[f([dim - 0.5 for dim in (0, mod.land.dim[0])]) for f in (plt.xlim,
    #                                                          plt.ylim)]
    # get rid of x and y ticks
    #[f([]) for f in (plt.xticks, plt.yticks)]

def plot_genetic_PCA(mod):
    """run a genetic PCA and plot individuals on first 2 PCs, colored by which
    side of the barrier they are on"""
    from copy import deepcopy
    from sklearn.decomposition import PCA
    figsize = 6
    species = mod.comm[0]
    # get array of resulting genomic data (i.e.
    # 'speciome'), genotypes meaned by individual
    speciome = np.mean(np.stack([i.g for i in species.values()]), axis=2)
    # run PCA on speciome
    pca = PCA(n_components=2)
    PCs = pca.fit_transform(speciome)
    # normalize the PC results
    norm_PCs = (PCs - np.min(PCs,
                    axis=0)) / (np.max(PCs, axis=0) - np.min(PCs, axis=0))
    # assign a value to each species, 0 or 1, indicating whether
    # they're located on the left or right vertical half of the landscape
    ind_colors = [0 if (
                ind.x < mod.land.dim[0]/2) else 1 for ind in species.values()]
    # plot individuals on PCs 1 and 2, colored by their landscape half
    fig = plt.figure(figsize=(figsize, figsize), dpi= 80,
                     facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    PC_dict = dict(zip(range(len(ind_colors)), ind_colors))
    left = [k for k,v in PC_dict.items() if v == 0]
    right = [k for k,v in PC_dict.items() if v == 1]
    patches = []
    patches.append(ax1.scatter(norm_PCs[left, 0], norm_PCs[left, 1],
                               color = '#00ffff'))
    patches.append(ax1.scatter(norm_PCs[right, 0], norm_PCs[right, 1],
                               color = '#ff00ff'))
    ax1.set_xlabel('genetic PC 1')
    ax1.set_ylabel('genetic PC 2')
    ax1.legend(patches, ['left of barrier', 'right of barrier'])



# calculate euclidean distance from two n-length vectors
def calc_euc(x, y):
    euc = np.sqrt(sum([(n-m)**2 for n, m in zip(x, y)]))
    return euc


# calculate lower-triangular of PCA-bsaed Euclidean genetic distances between
# all individuals, using a 'speciome' (2d array of all individs' genomes)
def calc_dists(species, dist_type='gen', env_lyrs=None, return_flat=True):
    # calculate genetic distance as the euclidean distance between individuals
    # in genetic PC space
    if dist_type == 'gen':
        speciome = np.mean(np.stack([i.g for i in species.values()]), axis=2)
        pca = PCA()
        vals = pca.fit_transform(speciome)
    # calculate geographic distance as the linear euclidean distance between
    # individuals
    elif dist_type == 'geo':
        vals = np.stack([np.array((i.x, i.y)) for i in species.values()])
    # calculate environmental distance as the euclidean distance between
    # individuals' environments, for all environmental layers specified by
    # the 'env_lyrs' argument
    elif dist_type == 'env':
        vals = np.stack([np.array(i.e)[env_lyrs] for i in species.values()])
    # print(vals)
    # print(vals.shape)
    n_ind = vals.shape[0]
    dist_mat = np.ones([n_ind] * 2) * -999
    # dist_vals = [[calc_euc(i, j) for j in vals[n:,
    #                                        :]] for n, i in enumerate(vals)]
    for i in range(n_ind):
        for j in range(0, i+1):
            dist_mat[i, j] = calc_euc(vals[i, :], vals[j, :])
    # check that all diagonal values are 0
    assert np.all(np.diag(dist_mat) == 0), "Not all diagonal values are 0!"

    if return_flat:
        # flatten the lower triangle, if return 1-d of values
        dists = dist_mat[np.tril_indices(dist_mat.shape[0], -1)]
        assert dists.size == (n_ind**2 - n_ind)/2, ("Length not equal "
                                                    "to n(n-1)/2!")
    else:
        # make it a symmetric dist matrix, if returning the matrix
        dist_mat[dist_mat == -999] = 0
        dists = dist_mat + dist_mat.T
        assert dists.size == n_ind**2, "Size not equal to n*n!"
    # dist_vals = [item[1:] for item in dist_vals]
    # dists = [*chain.from_iterable(dist_vals)]
    # assert that the length is correct
    return dists


# define a function to calculate the mean difference between phenotype
# and environment for a species
def calc_mean_z_e_diff(spp, trait_num=0):
    zs = spp._get_z().ravel()
    es = spp._get_e(lyr_num=spp.gen_arch.traits[trait_num].lyr_num)
    mean_diff = np.mean(np.abs(zs - es))
    return mean_diff


def track_horiz_crossing(mod, zone_edges, tracker, count):
    """track all horizontal crossing of some vertical zone on the landscape,
    delineated by x-axis bounds in zone_edges argument, so that the model can
    be run with and without the barrier mask, to quantify the extent to which
    the barrier increases landscape resistance
    """
    for idx, ind in mod.comm[0].items():
        # individ is assigned value 0 if on left side of the crossing zone...
        if ind.x <= zone_edges[0]:
            curr_side = 0
        #...or 1 if on right side...
        elif ind.x >= zone_edges[1]:
            curr_side = 1
        # or else nan if they are inside the zone
        else:
            curr_side = np.nan
        # if individ is in the tracker then they were alive last time step,
        # so grab and compare their previous value to this time step's value
        if idx in tracker:
            diff = curr_side - tracker[idx]
            # ignore nans, becuase they are currently inside the zone, so
            # whether or not they are crossing cannot be decided
            if np.isnan(diff):
                pass
            # increment the crossing counter, if a crossing has occurred
            elif diff != 0:
                count +=1
        # now update the individ's value in the tracker (unless they have a
        # nan, because they're currently inside the zone, so maintain the
        # value from their last time outside the zone until they exit
        # the zone and it can be determined whether or not they crossed
        if not np.isnan(curr_side):
            tracker[idx] = curr_side
    return count


# make empty figure
fig = plt.figure(figsize=(6.75, 9.25))
gs = gridspec.GridSpec(4, 2,
                       width_ratios=[1, 1],
                       height_ratios=[1, 1, 0.05, 0.5])

# create objects for crossing of barrier, and of equal-width vertical areas
# within each of the two sides for comparison
cross_count = 0
cross_tracker = {}

# start timer
start = time.time()

# make model
params = gnx.read_parameters_file(('./geonomics/demos/IBD_IBE/'
                                   'IBD_IBE_params.py'))

# get barrier-zone edges, for tracking crossings
barr_rast = params['landscape']['layers']['barrier']['init']['defined']['rast']
zone_edges = np.where(barr_rast[0,:] == 0)[0]
zone_edges = (zone_edges.min(), zone_edges.max() + 1)

# then get rid of barrier, if not to be used for this run
if not use_barr:
    params['landscape']['layers']['barrier']['init']['defined']['rast'] = np.ones(
                                        params['landscape']['main']['dim'])

mod = gnx.make_model(params)

# define number of timesteps
T = 1000

# burn model in
mod.walk(20000, 'burn')

# plot genetic PCA before evolution begins
# ax1 = fig.add_subplot(321)
ax1 = plt.subplot(gs[0])
ax1.set_aspect('equal')
# ax1.set_title('t = 0')
map_genetic_PCA(mod.comm[0], mod.land)
ax1.set_title('colored by genetic PCA', fontsize=16)
ax1.set_ylabel('before simulation', fontsize=14)

# plot phenotypes before evolution begins
# ax3 = fig.add_subplot(323)
ax2 = plt.subplot(gs[1])
ax2.set_title('colored by phenotype', fontsize=16)
mask = np.ma.masked_where(mod.land[1].rast == 0, mod.land[1].rast)
mod.plot_phenotype(0, 0, mask_rast=mask, size=mark_size, cbar=False)
[f((-0.5, 39.5)) for f in [plt.xlim, plt.ylim]]

# create data structure to save z-e diff values
if not_timing:
    z_e_diffs = []

# run model for T timesteps
for t in range(T):
    if not_timing:
        z_e_diffs.append(calc_mean_z_e_diff(mod.comm[0]))
        cross_count = track_horiz_crossing(mod, zone_edges, cross_tracker,
                                           cross_count)
    mod.walk(1)
if not_timing:
    z_e_diffs.append(calc_mean_z_e_diff(mod.comm[0]))
    cross_count = track_horiz_crossing(mod, zone_edges, cross_tracker,
                                       cross_count)

# finish timer
stop = time.time()
tot_time = stop - start


# plot genetic PCA after 1/4T timesteps
# ax2 = fig.add_subplot(322)
ax3 = plt.subplot(gs[2])
ax3.set_aspect('equal')
# ax3.set_title('t = %i' % T)
map_genetic_PCA(mod.comm[0], mod.land)
ax3.set_ylabel('after simulation', fontsize=14)

# plot the individuals' phenotypes
# ax4 = fig.add_subplot(324)
ax4 = plt.subplot(gs[3])
mod.plot_phenotype(0, 0, mask_rast=mask, size=mark_size, cbar=False)
[f((-0.5, 39.5)) for f in [plt.xlim, plt.ylim]]

# plot IBD and IBE
spp_subset = {ind: mod.comm[0][ind] for ind in np.random.choice([*mod.comm[0]],
                                                                100)}
gen_dists = calc_dists(spp_subset)
scaled_gen_dists = gen_dists/gen_dists.max()
assert (np.all(scaled_gen_dists >= 0)
        and np.all(scaled_gen_dists <= 1)), ('Scaled genetic dist is outside '
                                             '0 and 1!')
geo_dists = calc_dists(spp_subset, 'geo')
env_dists = calc_dists(spp_subset, 'env', [0])


# run multiple linear regression of gen on geo and env dists
mlr_est = sm.Logit(endog=np.array(scaled_gen_dists.T),
                   exog=np.vstack((geo_dists, env_dists)).T).fit()
#mlr_est = sm.GLM(endog=np.array(scaled_gen_dists.T),
#                 exog=np.vstack((geo_dists, env_dists)).T,
#                 family=sm.families.Binomial(sm.families.links.cloglog)).fit()


# ax5 = fig.add_subplot(325)
ax5 = plt.subplot(gs[6])
ax5.scatter(geo_dists, scaled_gen_dists, alpha=0.05, c='black')
# plt.scatter(geo_dists, gen_dists, alpha=0.05, c='black')
# add a regression line (NOTE: doesn't include an intercept by default,
# so I need to include one manually in the design matrix)
x_preds = np.arange(0, 50.1, 0.1)
y_preds = np.linspace(0, 1, len(x_preds))
z_preds = mlr_est.predict(np.vstack((x_preds, y_preds)).T)
# y_preds = est.predict(np.vstack(([1] * len(x_preds), x_preds)).T)
ax5.plot(x_preds, z_preds, color='#C33B3B')
ax5.text(min(geo_dists) + 0.4 * (max(geo_dists) - min(geo_dists)),
         0.15, 'slope:    %0.3f' % mlr_est.params[0],
         color='#C33B3B', size=9)
p_val_lt = mlr_est.pvalues[0] < 0.001
assert p_val_lt, 'p-value not less than 0.001!'
ax5.text(min(geo_dists) + 0.4 * (max(geo_dists) - min(geo_dists)),
         0.05, 'p-value < 0.001',
         color='#C33B3B', size=9)
#ax5.text(min(geo_dists) + 0.4 * (max(geo_dists) - min(geo_dists)),
#         0.45, 'Pseudo-$R^{2}$:   %0.2f' % mlr_est.prsquared,
#         color='#C33B3B', size=9)
ax5.set_title('IBD', fontsize=16)
ax5.set_xlabel('geographic distance', fontsize=14)
ax5.set_ylabel('genetic distance', fontsize=14)
# ax5.set_ylim(int(np.floor(min(gen_dists))), int(np.ceil(max(gen_dists))))
ax5.set_ylim((0, 1))
# ax5.set_aspect('equal')
# ax6 = fig.add_subplot(326, sharey=ax5)


ax6 = plt.subplot(gs[7], sharey=ax5)
ax6.scatter(env_dists, scaled_gen_dists, alpha=0.05, c='black')
# plt.scatter(env_dists, gen_dists, alpha=0.05, c='black')
# add a regression line (a quadratic regression, to allow for a saturating
# pattern within the range of x values)
# x_preds = np.arange(0, 1, 0.01)
# z_preds = est.predict(x_preds.T)
# y_preds = est.predict(np.vstack(([1] * len(x_preds), x_preds)).T)
ax6.plot(y_preds, z_preds, color='#C33B3B')
ax6.text(min(env_dists) + 0.4 * (max(env_dists) - min(env_dists)),
         0.15, 'slope:    %0.3f' % mlr_est.params[1],
         color='#C33B3B', size=9)
p_val_lt = mlr_est.pvalues[0] < 0.001
assert p_val_lt, 'p-value not less than 0.001!'
ax6.text(min(env_dists) + 0.4 * (max(env_dists) - min(env_dists)),
         0.05, 'p-value < 0.001',
         color='#C33B3B', size=9)
#ax6.text(min(env_dists) + 0.4 * (max(env_dists) - min(env_dists)),
#         0.45, 'Pseudo-$R^{2}$:   %0.2f' % mlr_est.prsquared,
#         color='#C33B3B', size=9)
ax6.set_title('IBE', fontsize=16)
ax6.set_xlabel('environmental distance', fontsize=14)
# ax6.set_yticks([])
# ax6.set_aspect('equal')

# TODO: add regression line to each plot
fig.tight_layout()
plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.96, wspace=0.07,
                    hspace=0.16)
plt.savefig(os.path.join(img_dir, 'IBD_IBE.pdf'), format='pdf', dpi=1000)


# create 3d axes, and label for geo, env, and gen
fig3d = plt.figure()
ax3d = fig3d.add_subplot(111, projection='3d')
# get 3d-scatter colors
col3d = ((geo_dists / geo_dists.max()) + (env_dists / env_dists.max())) / 2
# scatter all 3 vars on those axes
ax3d.scatter(geo_dists, env_dists, scaled_gen_dists, alpha=0.5, c=col3d,
             cmap='plasma')
ax3d.set_xlabel('geo', size=15)
ax3d.set_ylabel('env', size=15)
ax3d.set_zlabel('gen', size=15)
# create predicted surface, and add to 3d plot
y_vals = np.arange(0, 1.01, 0.01)
ys = np.hstack([list(y_vals) for _ in range(len(y_vals))])
xs = np.hstack([[n] * len(y_vals) for n in np.linspace(0, 50, len(y_vals))])
zs = mlr_est.predict(np.vstack((xs, ys)).T)
xs = xs.reshape([len(y_vals)] * 2)
ys = ys.reshape([len(y_vals)] * 2)
zs = zs.reshape([len(y_vals)] * 2)
# surf_cols = ((xs / xs.max()) + (ys / ys.max())) / 2
# surf_cols = np.int64(surf_cols * 255) + 1
# ax3d.plot_surface(xs, ys, zs, facecolors=surf_cols, alpha=0.5, cmap='plasma')
ax3d.plot_surface(xs, ys, zs, color='black', alpha=0.4)


# create plot of z-e diffs
z_e_fig = plt.figure()
ax = z_e_fig.add_subplot(111)
plt.plot(range(len(z_e_diffs)), z_e_diffs, color='#096075')
ax.set_xlabel('time')
ax.set_ylabel(('mean difference between individuals\' phenotypes and '
               'environmental values'))
plt.show()


# print out the zone-crossing results
cross_rate = cross_count / T
print(("\nThe barrier-zone crossing rate was %0.6f individuals "
      "per time step.\n") % cross_rate)

# print out time
print("\n\nModel ran in %0.2f seconds." % tot_time)
