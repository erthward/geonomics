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

from collections import Counter as C
from operator import itemgetter as ig
from scipy import interpolate
from shapely import geometry as g




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
    dens = interpolate.griddata(np.array(list(zip(list(gi/grid_mag), list(gj/grid_mag)))), window_dens, (new_gj, new_gi), method = 'cubic')
    dens = np.flipud(dens)



    if normalize_by != 'none':

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

    if max_val != None:
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
    dens = interpolate.griddata(np.array(list(zip(list(gj), list(gi)))), dens.flatten(), (new_gj, new_gi), method = 'cubic')
    dens = np.flipud(dens)



    #if normalize_by == census, then divide each density by total pop census size
    if normalize_by == 'census':
        N = pop.census()
        dens = dens/N
    elif normalize_by == 'none':
        pass

    else:  #POTENTIALLY ADD OTHER OPTIONS HERE, i.e. to normalize by starting size?
        pass 

    if normalize_by != 'none':

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

    if max_val != None:
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




def cell_string_method(pop, land, window_width = None):
    from collections import Counter as C
    from operator import itemgetter as ig
    from scipy import interpolate

    #use ww as shorthand
    ww = window_width

    #set ww to 1 if None
    if ww == None:
        ww = 1
    #otherwise make sure that both landscape dimensions are evenly divisible by ww
    else:
        assert (land.dims[0]%ww == 0 and land.dims[1]%ww == 0), 'land.dims not evenly divisible by window_width'

    #get the order of magnitude of the landscape dimension (to decide length argument for zfill, below)
    ord_mag_dims = max([len(str(land.dims[0])), len(str(land.dims[1]))])

    #tuple of dimensions for the counting grid
        #adding 2 to each dimension creates a grid that represents the actual landscape grid, with a 1-cell
        #buffer around it, from which the actual landscape grid can be extracted at the end, to avoid
        #interpolation edge effects
    ct_grid_dims = (int(land.dims[1]//ww)+2, int(land.dims[0]//ww)+2)
    
    #create meshgrid for the counting grid
    xx,yy = np.meshgrid(range(ct_grid_dims[1]), range(ct_grid_dims[0]))
    #create a list of the string versions of the zfill'd x and y coordinates for each cell
    xx_str = [str(x).zfill(ord_mag_dims) for x in xx.flatten()]
    yy_str = [str(y).zfill(ord_mag_dims) for y in yy.flatten()]
    #then concatenate the strings of the x and y coordinates for each cell, to create unique cell IDs
    cells = [''.join(c) for c in list(zip(yy_str,xx_str))]

    #create a function to extract the dictionary value for each key in an input list (to be called on the
    #Counter object below)
    f = ig(*cells) 
        # the post here: https://stackoverflow.com/questions/18453566/python-dictionary-get-list-of-values-for-list-of-keys
        #suggests this is the fastest approach

    #get population's coordinates
    c = np.array(pop.get_coords().values())
    #add one window's width to the coordinates, to shift them into their correct landscape locations within
    #the nested grid at the center of the buffered count grid
    true_x = c[:,0] + ww
    true_y = c[:,1] + ww

    offsets = [(0, 0), (0.5, 0), (0, 0.5), (0.5, 0.5)]

    #array to hold the results for all four coordinate offsets
    res_array = np.zeros((len(offsets), land.dims[0], land.dims[1]))

    all_x_pts = []
    all_y_pts = [] 
    all_res = []

    #loop over offsets
    for n, off in enumerate(offsets):

        x = true_x + (off[0]*ww)
        y = true_y + (off[1]*ww)

        #and get zfill'd strings of these coordinates as well  
        X = [str(int(int(j)//ww)).zfill(ord_mag_dims) for j in x]
        Y = [str(int(int(i)//ww)).zfill(ord_mag_dims) for i in y]
        locs = [''.join(i) for i in list(zip(Y,X))]


        #count the number of individuals within each grid cell
        cts = C(locs)

        #then use the function defined above to extract the counts each cell's unique string ID
        res = f(cts)

        #convert the res vector to an array, of shape determined by ct_grid_dims
        res = np.reshape(res, ct_grid_dims)
        #divide each cell's count by the cell area
        #res = res/(ww**2)

        #if the count grid has a resolution greater than 1 landscape cell, interpolate to the landscape resolution
        if ww > 1:
            #create a meshgrid of points at the centers of the landscape cells (plus 1 ww, to account for the
            #1-cell-buffer around the count grid)
            out_x, out_y = np.meshgrid(np.arange(land.dims[0])+0.5+ww, np.arange(land.dims[1])+0.5+ww)
           
            #get an array of the x,y coordinates associated with each count in the count grid
            x_pts = [(j + 0.5 + off[0]) * ww for j in xx.flatten()]
            y_pts = [(i + 0.5 + off[1]) * ww for i in yy.flatten()]
            pts = np.array(list(zip(y_pts, x_pts)))

            all_x_pts.extend(list(x_pts)) 
            all_y_pts.extend(list(y_pts)) 
            all_res.extend(list(res.flatten()))


            #then interpolate to the landscape cells
            dens = interpolate.griddata(pts, res.flatten(), (out_y, out_x), method = 'cubic')

            #floor the array at 0
            dens[dens<0] = 0

            #add to res_array
            res_array[n,:,:] = dens

        #otherwise, if the count grid has a 1-landscape-cell resolution, just return the center-nested count grid
        else:
            res_array[n,:,:] = res[1:ct_grid_dims[0]-1, 1:ct_grid_dims[1]-1]

    all_dens = interpolate.griddata(np.array(list(zip(all_y_pts, all_x_pts))), np.array(all_res), (out_y, out_x), method = 'cubic')
    
    dens = np.mean(res_array, axis = 0)
    return(res_array, dens, all_dens)





################################################################################################
#DEH: 04/07/18: ATTEMPTING TO USE THE CELL-STRING COUNTING APPROACH TO CREATE A Grid_Stack class
################################################################################################

class Grid:
    def __init__(self, dims, ww, gi, gj, cells, areas, corner_cells):

        self.dims = dims
        self.ww = ww
        
        self.dim_om = max([len(str(d)) for d in self.dims])

        self.corner_cells = corner_cells  #True if this is the outer grid (i.e. contains cells centered
                                          #on the landscape corners); else False for the inner grid (i.e.
                                          #cells centered on points half-window-width in from landscape corners 

        self.gi = gi
        self.gj = gj

        self.grid_coords = np.array(list(zip(self.gi.flatten(), self.gj.flatten())))
        #self.grid_coords, self.grid_areas = get_grid_coords_and_areas(self.gi, self.gj, self.dims)

        self.cells = cells
        self.areas = areas

        self.grid_counter = ig(*self.cells) 


    def get_dens(self, x, y):
        x = np.array(x)
        y = np.array(y)
        x_cells = (x - self.corner_cells*self.ww/2.)//self.ww + self.corner_cells
        y_cells = (y - self.corner_cells*self.ww/2.)//self.ww + self.corner_cells
        cells = get_cell_strings(y_cells, x_cells, self.dim_om)
        counts = C(cells)
        grid_counts = self.grid_counter(counts)
        grid_counts = np.reshape(grid_counts, self.gi.shape)
        grid_dens = grid_counts/self.areas
        return(grid_dens)



class Grid_Stack:
    def __init__(self, land, ww):

        self.dims = land.dims
        self.ww = ww
        self.dim_om = max([len(str(d)) for d in self.dims])

        self.land_gi, self.land_gj = np.meshgrid(np.arange(0, self.dims[0])+0.5, np.arange(0, self.dims[1])+0.5)

        self.grids = dict([(n,g) for n,g in enumerate(make_grids(land, self.ww))])


    def calc_density(self, x, y):
        pts = np.vstack([grid.grid_coords for grid in self.grids.values()])
        vals = np.hstack([grid.get_dens(x, y).flatten() for grid in self.grids.values()])
        dens = interpolate.griddata(pts, vals, (self.land_gi, self.land_gj), method = 'cubic')
        return(dens)



def get_cell_strings(gi, gj, dim_om):
    i_strs = [str(int(i)).zfill(dim_om) for i in gi.flatten()]
    j_strs = [str(int(j)).zfill(dim_om) for j in gj.flatten()]
    cells = [''.join(c) for c in list(zip(i_strs,j_strs))]
    return(cells)



def make_grids(land, ww):

    hww = ww/2.
    dims = land.dims
    
    dim_om = max([len(str(d)) for d in dims])

    gj1, gi1 = np.meshgrid(np.arange(0,dims[0]+ww,ww), np.arange(0,dims[1]+ww,ww))
    gj2, gi2 = np.meshgrid(np.arange(0+hww,dims[0]+hww,ww), np.arange(0+hww,dims[1]+hww,ww))
   
    land_poly_coords = ((0,0), (dims[0], 0), (dims[0], dims[1]), (0, dims[1]))
    land_poly = g.Polygon(land_poly_coords)
    
    j1 = gj1.flatten()
    i1 = gi1.flatten()
    j2 = gj2.flatten()
    i2 = gi2.flatten()

    polys1 = [g.Polygon(((j1[n]-hww, i1[n]-hww), (j1[n]-hww, i1[n]+hww), (j1[n]+hww, i1[n]+hww), (j1[n]+hww, i1[n]-hww))) for n in range(len(j1))]
    polys2 = [g.Polygon(((j2[n]-hww, i2[n]-hww), (j2[n]-hww, i2[n]+hww), (j2[n]+hww, i2[n]+hww), (j2[n]+hww, i2[n]-hww))) for n in range(len(j2))]
    #polys1 = [g.asPolygon(((j1[n]-hww, i1[n]-hww), (j1[n]-hww, i1[n]+hww), (j1[n]+hww, i1[n]+hww), (j1[n]+hww, i1[n]-hww))) for n in range(len(j1))]
    #polys2 = [g.asPolygon(((j2[n]-hww, i2[n]-hww), (j2[n]-hww, i2[n]+hww), (j2[n]+hww, i2[n]+hww), (j2[n]+hww, i2[n]-hww))) for n in range(len(j2))]

    areas1 = np.reshape([p.intersection(land_poly).area for p in polys1], gj1.shape)
    areas2 = np.reshape([p.intersection(land_poly).area for p in polys2], gj2.shape)

    i_cells_1 = (i1 - ww/2.)//ww + 1
    j_cells_1 = (j1 - ww/2.)//ww + 1
    i_cells_2 = (i2)//ww 
    j_cells_2 = (j2)//ww 

    cells1 = get_cell_strings(i_cells_1, j_cells_1, dim_om)
    cells2 = get_cell_strings(i_cells_2, j_cells_2, dim_om)

    g1 = Grid(dims, ww, gi1, gj1, cells1, areas1, corner_cells = True)
    g2 = Grid(dims, ww, gi2, gj2, cells2, areas2, corner_cells = False)

    return(g1, g2)




def compare_original_new(land, pop, grid_stack):
    x = list(pop.get_x_coords().values())
    y = list(pop.get_y_coords().values())
    d1 = pop.calc_density(land, window_width = 5).raster
    d2 = grid_stack.calc_density(x,y)

    fig = plt.figure()
    fig.suptitle('Comparison of pop-density calculation approaches')

    ax = fig.add_subplot(131)
    ax.set_title('ORIG:\n 50x50 x 5-cell wind\nx 1675inds: 256ms')
    plt.imshow(d1, interpolation = 'nearest', cmap = 'PiYG')
    plt.colorbar()
    plt.xlim((0,50))
    plt.ylim((0,50))
    plt.scatter(np.array(x)-0.5, np.array(y)-0.5, color = 'black', s = 5)

    ax = fig.add_subplot(132)
    ax.set_title('ORIG - NEW')
    plt.imshow(d1 - d2.T, interpolation = 'nearest', cmap = 'PiYG')
    plt.colorbar()
    plt.xlim((0,50))
    plt.ylim((0,50))
    plt.scatter(np.array(x)-0.5, np.array(y)-0.5, color = 'black', s = 5)

    ax = fig.add_subplot(133)
    ax.set_title('NEW:\n 50x50 x 5-cell wind\nx 1675inds: 12.4ms')
    plt.imshow(d2.T, interpolation = 'nearest', cmap = 'PiYG')
    plt.colorbar()
    plt.xlim((0,50))
    plt.ylim((0,50))
    plt.scatter(np.array(x)-0.5, np.array(y)-0.5, color = 'black', s = 5)

    plt.show()
