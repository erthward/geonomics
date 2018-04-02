#!/usr/bin/python


'''NOTE: We need to rethink how the density raster is calculated. Right now, moving-window counts are taken,
then those are interpolated to a raster, but this performs poorly at the small window sizes needed to get
reasonable density estimates. I tried (below) using shapely to count points within windows, rather than the <
and > x and y checks I used before, but it performed about 10x more slowly. Might need to rethink completely.
Ideas include:

    - porting a simplistic scatterplot to a buffer, then overlaying that on the lanscape and scaling the
      values from light to dark to turn them into density estimates (but still doesn't give us number of
      individuals, unless we can back out number of individuals from pixel values)
      example code:

      frame = plt.scatter([v[0] for v in c], [v[1] for v in c], color = '#000000', alpha = 0.2, s = 280,
      edgecolor = None, marker = 's'); plt.xlim((0,50)); plt.ylim((0,50))
      frame.axes.xaxis.set_ticklabels([]); frame.axes.yaxis.set_ticklabels([])
      frame.axes.xaxis.set_ticks([]); frame.axes.yaxis.set_ticks([])
      ax.axis('off')
      buf = io.BytesIO()
      plt.savefig(buf, format = 'PNG')
      im = Image.open(buf)
      im.show()

    - building a matrix for each individual, where the minimum total difference between window centers and the
      individual's coordinates is a 1, the rest are 0s, then summing all individual's matrices. Steps:

      pts = np.linspace(0.5, 49.5, 50)
      c = pop.get_coords().values()
      window_width = 5
      y_yes = sqrt((p[1]-pts)**2)<=(window_width/2)
      x_yes = sqrt((p[0]-pts)**2)<=(window_width/2)
      - somehow turn those two into an array with 1s where both x and y are True
      - then generate a stack of such arrays, one for each individual
      - then sum that stack across indviduals to get a density raster

      - or instead, use the x_yes and y_yes to generate points to then interpolate to a grid


'''




def alt_calc_density(pop, land, grid_mag = 1, count_type = 0, window_width = None, normalize_by = 'none', min_0 = True, max_1 = False, max_val = None, set_N = False): 

    from scipy import interpolate
    from shapely import geometry as g

    #window width defaults to 1/10 the maximum landscape dimension
    if window_width == None:
        window_width = max(land.dims)*0.1


    #shorthand
    dims = land.dims
    mag_dims = tuple([grid_mag * n for n in dims])

    #get a list of pop's coord-tuples
    c = pop.get_coords().values() 
    pt_c = g.MultiPoint(c)
    tot = len(c)

    #multiply by grid-magnification, and make window_width a float, to avoid Py2 integer-division issues
    #window_width = grid_mag * window_width
    window_width = float(window_width)
    
    #create meshgrid using window_width/2 as step size
    grid_j, grid_i = np.mgrid[0:mag_dims[0]:complex("%ij" % (mag_dims[0]/(window_width/2))), 0:mag_dims[1]:complex("%ij" % (mag_dims[1]/(window_width/2)))]

    #grid_j, grid_i = np.mgrid[0+(window_width/2):dims[0]-(window_width/2):complex("%ij" % (dims[0]/(window_width/2))), 0+(window_width/2):dims[1]-(window_width/2):complex("%ij" % (dims[1]/(window_width/2)))] 

    #flatten the arrays, so that I can run over them in a single for loop
    gj = grid_j.ravel()
    gi = grid_i.ravel()

    #make lists of tuples, of same length as gj, containing the window ll and ur coords
    window_ll = [(max(gj[n]-(window_width/2), 0), max(gi[n]-(window_width/2), 0)) for n in range(len(gj))]   #constrain min window vals to 0
    window_ur = [(min(gj[n]+(window_width/2), mag_dims[0]), min(gi[n]+(window_width/2), mag_dims[1])) for n in range(len(gj))] #constrain max window vals to each respective land dimension
    assert len(window_ll) == len(gj)
    assert len(window_ur) == len(gj)

    poly_window = [g.asPolygon((window_ll[i], (window_ll[i][0], window_ur[i][1]), window_ur[i], (window_ur[i][0], window_ll[i][1]))) for i in range(len(window_ll))]

    #make a list of the sizes of each window
    window_size = [(window_ur[n][0] - window_ll[n][0]) * (window_ur[n][1] - window_ll[n][1]) for n in range(len(gj))]#calculate size of this window (not always the same because of edge effects
    assert len(window_size) == len(gj)
   
    if count_type == 0:
        #make a list of the counts of all individs within each window
        #NOTE: DEH 03-17-18: THIS IS THE LINE THAT NEEDS TO BE SPED UP!
        window_ct = [len([ind for ind in range(len(c)) if (c[ind][0]*grid_mag>window_ll[n][0] and c[ind][0]*grid_mag<=window_ur[n][0]) and (c[ind][1]*grid_mag>window_ll[n][1] and c[ind][1]*grid_mag<=window_ur[n][1])]) for n in range(len(gj))] 
        assert len(window_ct) == len(gj)

    elif count_type == 1:
        window_ct = [tot - len(pt_c.difference(wind)) for wind in poly_window]
        assert len(window_ct) == len(gj)


    #divide window counts by window sizes
    window_dens = [window_ct[n]/window_size[n] for n in range(len(window_ct))] #divide by window size
    assert len(window_dens) == len(gj)

    #if normalize_by == census, then divide each density by total pop census size
    if normalize_by == 'census':
        N = pop.census()
        window_dens = [dens/N for dens in window_dens]
    elif normalize_by == 'none':
        pass

    else:  #POTENTIALLY ADD OTHER OPTIONS HERE, i.e. to normalize by starting size?
        pass 

    #interpolate resulting density vals to a grid equal in size to the landscape
    new_gj, new_gi = np.mgrid[0:dims[0]-1:complex("%ij" % (dims[0])), 0:dims[1]-1:complex("%ij" % (dims[1]))]
    dens = interpolate.griddata(np.array(zip(list(gi/grid_mag), list(gj/grid_mag))), window_dens, (new_gj, new_gi), method = 'cubic')
    dens = np.flipud(dens)



    if normalize_by <> 'none':

        #if max_1 == True, set max_val to dens.max(), such that the density raster output will be normalized to
        #its own max, and thus vary between 0 and 1; else set to 1, and the output raster will vary between 0 and the current max value
        if max_1 == True:
            max_val = dens.max()
        elif max_1 == False:
            max_val = 1

        #Use max_val to normalize the density raster to either 0 to its current max val or
        #0 to 1, to make sure the interpolation didn't generate any values slightly outside this range
        norm_factor = max_val - dens.min()
        dens = (dens - dens.min())/norm_factor

    if min_0 == True:
        dens[dens<0] = 0

    if max_val <> None:
        dens[dens>max_val] = max_val

    if set_N == True:
        pop.set_N(landscape.Landscape(dims, dens))
    
    else:
        return(landscape.Landscape(dims, dens))



def create_circle(r, center):
    rads = np.arange(0,2*pi,0.01)
    xs = [center[0] + np.cos(rad)*r for rad in rads]
    ys = [center[1] + np.sin(rad)*r for rad in rads]
    circle = Polygon([xs, ys])
    return(circle)


def create_rectangle(dims):
    from shapely import geometry as g
    ll = (0,0)
    lr = (dims[1],0)
    ur = dims
    ul = (0, dims[0])
    rect = g.Polygon([ll, lr, ur, ul])
    return(rect)



def create_search_area_array(r, dims):
    from shapely import geometry as g
    landscape_rect = create_rectangle(dims)
    i_pts = np.arange(0 + r/2., dims[0]-(r/2.), r)
    j_pts = np.arange(0 + r/2., dims[1]-(r/2.), r)
    search_areas = []
    max_search_area = pi*(r**2)
    for i in i_pts:
        for j in j_pts:
            #create search-area circle around the point, take its difference with the landscape rectangle, and
            #subtract from the full search radius area to get the effective search radius at each cell-center
            #point, appending that to search_areas
            search_areas.append(max_search_area - g.Point((j, i)).buffer(r,128).difference(landscape_rect).area)

    search_area_array = np.reshape(np.asarray(search_areas), (len(i_pts), len(j_pts)))
    return(search_area_array)







def alt2_calc_density(pop, land, search_area_array, window_width = None, grid_mag = 1, count_type = 0, normalize_by = 'none', min_0 = True, max_1 = False, max_val = None, set_N = False): 

    #NOTE: REALLY, I COULD CREATE THE SEARCH_AREA_ARRAY AT THE BEGINNING OF THE MODEL AND ATTACH IT AS AN
    #ATTRIBUTE OF LANDSCAPE_STACK, THEN JUST GET COUNTS AND DIVIDE THEM AT EACH TIMESTEP

    from scipy import interpolate
    from shapely import geometry as g

    #window width defaults to 1/10 the maximum landscape dimension
    if window_width == None:
        window_width = max(land.dims)*0.1
    hww = window_width/2.

    #shorthand
    dims = land.dims

    #create a list of points at cell centers
    i_pts = np.arange(0 + hww/2., dims[0]-(hww/2.), hww)
    j_pts = np.arange(0 + hww/2., dims[1]-(hww/2.), hww)

    
    
    grid_j, grid_i = np.mgrid[0+hww/2.:dims[0]-(hww/2.):hww, 0+hww/2.:dims[1]-(hww/2.):hww]
    gj = grid_j.ravel()
    gi = grid_i.ravel()


    #get a list of pop's coord-tuples
    c = pop.get_coords().values() 

    #get cellwise local population counts by generating a Boolean array indicating whether an individual is within a window surrounding each cell's centerpoint, for each individual, then summing across individuals
    counts = np.sum([np.vstack([abs(p[0]-j_pts)<=hww * row for row in abs(p[1]-i_pts)<=hww]) for p in c], axis = 0)

    #then divide each cell by the area of the window surrounding that cell 
    dens = counts/search_area_array
    #assert dens.shape == dims


    new_gj, new_gi = np.mgrid[0:dims[0]-1:complex("%ij" % (dims[0])), 0:dims[1]-1:complex("%ij" % (dims[1]))]
    dens = interpolate.griddata(np.array(zip(list(gj), list(gi))), dens.flatten(), (new_gj, new_gi), method = 'cubic')
    dens = np.flipud(dens)



    #if normalize_by == census, then divide each density by total pop census size
    if normalize_by == 'census':
        N = pop.census()
        dens = dens/N
    elif normalize_by == 'none':
        pass

    else:  #POTENTIALLY ADD OTHER OPTIONS HERE, i.e. to normalize by starting size?
        pass 

    if normalize_by <> 'none':

        #if max_1 == True, set max_val to dens.max(), such that the density raster output will be normalized to
        #its own max, and thus vary between 0 and 1; else set to 1, and the output raster will vary between 0 and the current max value
        if max_1 == True:
            max_val = dens.max()
        elif max_1 == False:
            max_val = 1

        #Use max_val to normalize the density raster to either 0 to its current max val or
        #0 to 1, to make sure the interpolation didn't generate any values slightly outside this range
        norm_factor = max_val - dens.min()
        dens = (dens - dens.min())/norm_factor

    if min_0 == True:
        dens[dens<0] = 0

    if max_val <> None:
        dens[dens>max_val] = max_val

    if set_N == True:
        pop.set_N(landscape.Landscape(dims, dens))
    
    else:
        return(landscape.Landscape(dims, dens))






def calc_kde_pop_density(pop, land, bw = 0.1, buff = 5):
    from scipy.stats import gaussian_kde as kde
    dims = land.dims
    xmin = 0.5
    ymin = 0.5
    xmax = dims[1]-0.5 + (buff*2)
    ymax = dims[0]-0.5 + (buff*2)
    xx, yy = np.mgrid[xmin:xmax:complex('%ij' % (dims[0]+(2*buff))), ymin:ymax:complex('%ij' % (dims[0]+(2*buff)))]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    c = np.array(pop.get_coords().values()).T+buff
    #NOTE: NEED TO BETTER UNDERSTAND, EXPLORE, AND CHOOSE THE BANDWIDTH OPTIONS
    res = kde(c, bw)
    dens = np.reshape(res(positions), xx.shape).T[buff:(dims[0]+buff), buff:(dims[1]+buff)]

    #NOTE: NEEDS TO BE SCALED TO GIVE A GOOD APPROXIMATION OF THE INDIVIDS-PER-CELL VALUES!
        #DON'T KNOW IF THIS APPROACH IS ACTUALLY JUSTIFIED...
    dens = dens*pop.Nt[::-1][0]

    #NOTE: NEED TO FIGURE OUT HOW TO AVOID EDGE EFFECTS!
   
    assert dens.shape == dims, 'dens.shape is %s' % str(dens.shape)
    return(landscape.Landscape(dims, dens))






class MidpointNormalize(mpl.colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed
    midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge
        # cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def compare_orig_and_kde(pop, land, window_width = None, bw = 0.15):
    c = np.array(pop.get_coords().values())
    orig = pop.calc_density(land, window_width = window_width)
    kde = calc_kde_pop_density(pop, land, bw = bw)
    fig = figure()
    fig.suptitle('window_width = %i, bw = %s' % (window_width, str(bw)))
    ax = fig.add_subplot(131)
    plt.imshow(orig.raster, interpolation='nearest', cmap = 'terrain')
    plt.colorbar() 
    plt.scatter(c[:,0]-0.5, c[:,1]-0.5)
    ax.set_title('Original')
    ax2 = fig.add_subplot(132)
    plt.imshow(kde.raster, interpolation='nearest', cmap = 'terrain')
    plt.colorbar()
    plt.scatter(c[:,0]-0.5, c[:,1]-0.5)
    ax2.set_title('KDE')
    ax3 = fig.add_subplot(133)
    max_val = (orig.raster - kde.raster).max()
    plt.imshow(orig.raster - kde.raster, interpolation='nearest', cmap = 'PiYG', norm=MidpointNormalize(midpoint=0,vmin=-1*max_val, vmax=max_val))
    plt.colorbar()
    plt.scatter(c[:,0]-0.5, c[:,1]-0.5, color = 'gray', alpha = 0.25)
    plt.show()
    return


