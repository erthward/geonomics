#!/usr/bin/python
#landscape.py

'''
##########################################

Module name:          landscape


Module contains:
                      - definition of the Landscape type
                      - function for creating a random landscape, based on input parameters
                      - associated functions


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           12-28-15
Documentation:        URL


##########################################
'''

from scipy import interpolate
import numpy as np
import numpy.random as r
import matplotlib as mpl
import matplotlib.pyplot as plt



#------------------------------------
# CLASSES ---------------------------
#------------------------------------


class Landscape:
    def __init__(self, dims, raster):
        self.dims = dims
        self.raster = raster
        self.island_mask = False
        self.mask_island_vals = False
        assert type(self.dims) in [tuple, list], "dims must be expressed on a tuple or a list"
        assert type(self.raster) == np.ndarray, "raster should be a numpy.ndarray"


    #####################
    ### OTHER METHODS ###
    #####################


    def show(self, colorbar = True, im_interp_method = 'nearest', pop = False):
        if plt.get_fignums():
            colorbar = False
        cmap = plt.cm.terrain
        vmin = 0
        if self.mask_island_vals == True:
            cmap.set_under(color = 'black')
            vmin = 1e-7
        #TODO: REALLY GET TO THE BOTTOM OF THE RASTER-PLOTTING ISSUE!! 04/16/17: While working on the 
            #pop-density raster function, I was plotting the results and realized that the function seemed to
            #be working, but that the results seemed reflected about the diagonal. I am doing this, for the
            #moment, as a TEMPORARY fix, but will have to make sure that I haven't already 'fixed' this (and
            #forgotten) elsewhere, such that this is actually introducing more problems in the long-run
        if pop == True:
            plt.imshow(self.raster, interpolation = im_interp_method, cmap = cmap, vmin = vmin)
        else:
            plt.imshow(np.flipud(self.raster), interpolation = im_interp_method, cmap = cmap, vmin = vmin)
        if colorbar:
            if self.raster.max() > 1:
                plt.colorbar(boundaries=np.linspace(0,self.raster.max(),51))
            else:
                plt.colorbar(boundaries=np.linspace(0,1,51))



    def zoom(self, min_i, max_i, min_j, max_j, colorbar = True, im_interp_method = 'nearest', pop = False):
        if plt.get_fignums():
            colorbar = False
        cmap = 'terrain'
        zoom_rast = np.array([row[min_j:max_j] for row in self.raster[min_i:max_i]]) 
        if pop == True:
            plt.imshow(zoom_rast, interpolation = im_interp_method, cmap = cmap)
        else:
            plt.imshow(np.flipud(zoom_rast), interpolation = im_interp_method, cmap = cmap)
        if colorbar:
            if zoom_rast.max() > 1:
                plt.colorbar(boundaries=np.linspace(0,zoom_rast.max(),51))
            else:
                plt.colorbar(boundaries=np.linspace(0,1,51))







class Landscape_Stack:
    def __init__(self, raster_list):
        self.scapes = dict(zip(range(len(raster_list)), raster_list))  
        assert False not in [raster.__class__.__name__ == 'Landscape' for raster in raster_list], 'All landscapes supplied in raster_list must be of type landscape.Landscape.'
        self.n_scapes = len(raster_list)
        assert len(set([land.dims for land in self.scapes.values()])) == 1, 'Dimensions of all landscapes must be even.'
        self.dims = self.scapes.values()[0].dims
        self.movement_surf = None
        self.n_movement_surf_surf = None
        self.n_island_mask_scape = None



    def show(self, scape_num = None, colorbar = True, im_interp_method = 'nearest', pop = False):
        if plt.get_fignums():
            colorbar = False
        cmaps = ['terrain'] + ['bone'] * 10
        alphas = [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

        if scape_num <> None:

            cmap = plt.cm.terrain
            vmin = 0
            if self.scapes[scape_num].mask_island_vals == True:
                cmap.set_under(color = 'black')
                vmin = 1e-7

            if pop == True:
                plt.imshow(self.scapes[scape_num].raster, interpolation = im_interp_method, cmap = cmap, vmin = vmin)
            else:
                plt.imshow(np.flipud(self.scapes[scape_num].raster), interpolation = im_interp_method, cmap = cmap, vmin = vmin)


            if colorbar:
                if self.scapes[scape_num].raster.max() > 1:
                    plt.colorbar(boundaries=np.linspace(0,self.scapes[scape_num].raster.max(),51))
                else:
                    plt.colorbar(boundaries=np.linspace(0,1,51))


        else:
            for n, scape in self.scapes.items():
                if pop == True:
                    plt.imshow(scape.raster, interpolation = im_interp_method, alpha = alphas[n], cmap = cmaps[n] )
                else:
                    plt.imshow(np.flipud(self.scapes[scape_num].raster), interpolation = im_interp_method, alpha = alphas[n], cmap = cmaps[n])

                if colorbar:
                    if self.scapes[n].raster.max() > 1:
                        plt.colorbar(boundaries=np.linspace(0,self.scapes[n].raster.max(),51))
                    else:
                        plt.colorbar(boundaries=np.linspace(0,1,51))



    def zoom(self, min_i, max_i, min_j, max_j, scape_num = None, colorbar = True, im_interp_method = 'nearest', pop = False):
        if plt.get_fignums():
            colorbar = False
        cmaps = ['terrain', 'bone']
        alphas = [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        cmap = 'terrain'
        if scape_num <> None:
            zoom_rast = np.array([row[min_j:max_j] for row in self.scapes[scape_num].raster[min_i:max_i]]) 
            if pop == True:
                plt.imshow(zoom_rast, interpolation = im_interp_method, cmap = cmap)
            else:
                plt.imshow(np.flipud(zoom_rast), interpolation = im_interp_method, cmap = cmap)
            if colorbar:
                if zoom_rast.max() > 1:
                    plt.colorbar(boundaries=np.linspace(0,zoom_rast.max(),51))
                else:
                    plt.colorbar(boundaries=np.linspace(0,1,51))
        else:
            for n, scape in self.scapes.items():
                zoom_rast = np.array([row[min_j:max_j] for row in self.scapes[n].raster[min_i:max_i]]) 
                if pop == True:
                    plt.imshow(zoom_rast, interpolation = im_interp_method, alpha = alphas[n], cmap = cmaps[n] )
                else:
                    #NOTE: FIGURE OUT WHY FLIPUD SEEMED GOOD WHEN PLOTTING POP, BUT BAD FOR CIRC_HIST FN BELOW
                    plt.imshow(np.flipud(zoom_rast), interpolation = im_interp_method, alpha = alphas[n], cmap = cmaps[n])

                if colorbar:
                    if zoom_rast.max() > 1:
                        plt.colorbar(boundaries=np.linspace(0,zoom_rast.max(),51))
                    else:
                        plt.colorbar(boundaries=np.linspace(0,1,51))






    def plot_movement_surf_vectors(self, params, circle = False):
        if params['movement_surf'] == True and self.movement_surf <> None:
            import movement
            movement.plot_movement_surf_vectors(self, params, circle = circle)
        else:
            print('\n\nThis Landscape_Stack appears to have no movement surface.\n\n')
            pass



    def plot_vonmises_mix_circ_draws(self, i, j, params):
        if self.movement_surf is None:
            print('This landscape stack appears to have no movement surface layer. Function not valid.')
            return
        else:
            scape_num = params['n_movement_surf_scape']
            pts = [(np.cos(a), np.sin(a)) for a in [self.movement_surf[i][j]()[0] for n in range(1000)]]
            plt.scatter([pt[0]*0.5+i for pt in pts], [pt[1]*0.5+j for pt in pts], color = 'red', alpha = 0.2)
            self.scapes[scape_num].zoom(max(i-10, 0), min(i+10, self.dims[0]), max(j-10, 0), min(j+10, self.dims[1]))



    def plot_vonmises_mix_hist(self, i, j):
        if self.movement_surf is None:
            print('This landscape stack appears to have no movement surface layer. Function not valid.')
            return
        else:
            plt.hist([self.movement_surf[i][j]()[0] for n in range(10000)], bins = 100, normed = True, alpha = 0.5)


    def plot_vonmises_mix_circ_hist(self, i, j, zoom, params, scale_fac = 4.5, color = 'black'):
        if self.movement_surf is None:
            print('This landscape stack appears to have no movement surface layer. Function not valid.')
            return
        else:
            scape_num = params['n_movement_surf_scape']
            v,a = np.histogram([self.movement_surf[i][j]()[0] for n in range(7500)], bins = 15)
            v = v/float(v.sum())
            a = [(a[n]+a[n+1])/2 for n in range(len(a)-1)]
            x = [np.cos(a[n])*0.5 for n in range(len(a))]
            y = [np.sin(a[n])*0.5 for n in range(len(a))]
            x = np.array(x)*v*scale_fac
            y = np.array(y)*v*scale_fac
            [plt.plot((j,(j+x[n])), (i,(i+y[n])), linewidth = 2, color = color) for n in range(len(x))] 
            self.scapes[scape_num].zoom(max(i-zoom, 0), min(i+zoom, self.dims[0]), max(j-zoom, 0), min(j+zoom, self.dims[1]), pop =True)




    #method for pickling a landscape stack
    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)













#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------


def random_surface(dims, n_rand_pts, interp_method = "cubic", island_val = 0, num_hab_types = 2, dist = 'beta', alpha = 0.05, beta = 0.05):  #requires landscape to be square, such that dim = domain = range
    #NOTE: can use "nearest" interpolation to create random patches of habitat (by integer values); can change num_hab_types to > 2 to create a randomly multi-classed landscape
    #NOTE: I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions
    #NOTE: This seems to generate decent, believable random rasters! 
        # n_rand_pts/dim ratio:
            # ~0.01-0.05 --> broad, simple gradients
            # ~0.05-0.10 --> slightly more complex, occasionally landscape of low regions divided by a major high region (or vice versa)
            # ~0.10-0.50 --> landscape can be broken up into numerous fragmented patches (though sometimes still relatively homogeneous, with one or two small, extremely different patches
            # ~0.50-1.00 --> more highly fragmented version of above
    max_dim = max(dims)
    if interp_method == 'nearest':
        if dist == 'unif':
            vals = r.rand(n_rand_pts) * (num_hab_types-1)
        elif dist == 'beta':
            vals = r.beta(alpha, beta, n_rand_pts) * (num_hab_types-1)
    else:
        if dist == 'unif':
            vals = r.rand(n_rand_pts)
        elif dist == 'beta':
            vals = r.beta(alpha, beta, n_rand_pts)
    pts = r.normal(max_dim/2, max_dim*2,[n_rand_pts,2]) #selects seed points from well outside the eventaul landscape, to ensure interpolation across area of interest
    grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim), 1:max_dim:complex("%ij" % max_dim)] #by this function's definition, the use of complex numbers in here specifies the number of steps desired
    I = interpolate.griddata(pts, vals, (grid_x, grid_y), method = interp_method)
    if interp_method == 'nearest':  #i.e., if being used to generate random habitat patches...
        I = I.round().astype(float)
    if interp_method == 'cubic':  #transform to constrain all values to 0 <= val <= 1
        I = I + np.abs(I.min())+(0.01*r.rand()) #NOTE: adding a small jitter to keep values from reaching == 0 or == 1, as would likely be the case with linear interpolation
        I = I/(I.max()+(0.01*r.rand()))
    if dims[0] <> dims[1]:
        pass #NOTE: figure out how to get it to use the dims tuple to subset an approriate size if dims not equal

    return Landscape(dims,I)



def defined_surface(dims, pts, vals, interp_method = "cubic", num_hab_types = 2):  #pts should be provided as n-by-2 Numpy array, vals as a 1-by-n Numpy array

    #NOTE: There seem to be some serious issues with all of this code, because the resulting landscapes are not quite symmetrical; and making small tweaks (e.g. adding 5 to all input points' coordinates) doesn't just tweak the output landscape but instead completely gets ride of it; I have an intuition that it comes from the code that coerces all raster values to 0 <= val <= 1, becuase that doesn't look it does quite what I originally intended for it to do, but I'm not really sure... anyhow, for now it works for my initial testing purposes
    

    #NOTE: like the random_surface function, this also requires landscape to be square, but I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions


    #NOTE: if discrete habitat patches desired, values should still be fed in as proportions, and will then be multipled by num_hab_types to develop habitat class values
    if interp_method == 'nearest':
        vals = vals * (num_hab_types-1)


    #add 0.5 to all pts, to center them with respect to the raster display, to make intended symmetry actually symmetrical
    #pts = pts + 0.5


    max_dim = max(dims)
    grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim), 1:max_dim:complex("%ij" % max_dim)]
    I = interpolate.griddata(pts, vals, (grid_x, grid_y), method = interp_method)
    if interp_method == 'nearest':  #i.e., if being used to generate random habitat patches...
        I = I.round().astype(float)
    if interp_method == 'cubic':  #transform to constrain all values to 0 <= val <= 1
        I = I + np.abs(I.min())+(0.01*r.rand()) #NOTE: adding a small jitter to keep values from reaching == 0 or == 1, as would likely be the case with linear interpolation
        I = I/(I.max()+(0.01*r.rand()))
    if dims[0] <> dims[1]:
        pass #NOTE: figure out how to get it to use the dims tuple to subset an approriate size if dims not equal
    return Landscape(dims,I)




def build_scape_stack(params, num_hab_types = 2):

    #NOTE: If a multi-class (rather than binary) block-habitat raster would be of interest, would need to make num_hab_types customizable)


    #grab necessary parameters from the params dict

    if params['num_scapes'] == None:
        num_scapes = 1
    else:
        num_scapes = params['num_scapes']


    if params['interp_method'] == None:
        interp_method = ['cubic'] * num_scapes
    else:
        interp_method = params['interp_method']


    dims = params['dims']






    #create rasters for random landscape, if params['rand_land'] == True
    if params['rand_land'] == True:

        n_rand_pts = params['n_rand_pts']

        #if only a single integer provided for n_rand_pts, then use that for all landscape layers (i.e. scale of spatial heterogeneity will be roughly equal for all layers); otherwise, a list or tuple of n_rand_pts could create different scales of heterogeneity for different landscape layers
        if type(n_rand_pts) == int:
            n_rand_pts = [n_rand_pts] * num_scapes




        land = Landscape_Stack([random_surface(dims, n_rand_pts[n], interp_method = interp_method[n], num_hab_types = num_hab_types) for n in range(num_scapes)])





    #or create rasters for defined landscape, if params['rand_land'] == False
    elif params['rand_land'] == False:

        #get pts
        pts = params['landscape_pt_coords']
        #if only a single array of pts provided, multiply them into a list to use for each separate landscape layer
        if type(pts) == np.ndarray:
            pts = [pts]*num_scapes


        #get vals
        vals = params['landscape_pt_vals']


        land = Landscape_Stack([defined_surface(dims, pts[n], vals[n], interp_method = interp_method[n], num_hab_types = num_hab_types) for n in range(num_scapes)])

   


    #if params['island_val'] > 0, then use the movement_surf raster to create T/F mask-raster, added as an
    #additional landscape, to be used to kill all individuals straying off 'islands'
    
    #create a movement surface, if params call for it  
    #NOTE: THIS WILL TAKE A WHILE TO COMPUTER UP-FRONT!
    if params['movement_surf'] == True:

        land.n_movement_surf_scape = params['n_movement_surf_scape']

        if params['islands'] == True:

            if params['island_val'] > 0:
                iv = params['island_val']
         
                #zero out the appropriate parts of the movement raster
                movement_rast = land.scapes[land.n_movement_surf_scape].raster
                #zero out the appropriate parts of the movement raster (but not exactly 0, because 
                #creates divison-by-zero problems in the pop_dynamics calculations)
                movement_rast[movement_rast < iv] = 1e-8
               

            elif ('island_mask' in params.keys() and params['island_mask'] <> None):
                im_file = params['island_mask']
                ma = np.loads(im_file)
                assert type(ma) == np.ndarray, "The pickled file located at the path provided in params['island_mask'] does not appear to be a numpy ndarray."

                #get the movement raster scape number
                movement_rast = land.scapes[land.n_movement_surf_scape].raster
                #zero out the appropriate parts of the movement raster (but not exactly 0, because 
                #creates divison-by-zero problems in the pop_dynamics calculations)
                movement_rast[ma == 0] = 1e-8
                
               
            #replace the movement_surf_scape raster with the updated raster with all outside-habitat values set to 1e-8
            land.scapes[land.n_movement_surf_scape].raster = movement_rast
            #set the Landscape.mask_island_vals flag on the movement surf to True, for plotting purposes
            land.scapes[land.n_movement_surf_scape].mask_island_vals  = True


            #then create an island mask and add as the last landscape (to use to quickly cull individuals outside the 'islands')
            island_mask = create_island_mask(dims, land.scapes[land.n_movement_surf_scape].raster)
            island_mask.island_mask = True  #set Lanscape.island_mask attribute to True
            land.scapes[land.n_scapes] = island_mask #set as the last scape
            land.n_island_mask_scape = land.n_scapes #the Landscape_Stack n_island_mask_scape attribute to last scape num
            land.n_scapes += 1 #then increment total scape nums by 1



 
        #create the movement surface, and set it as the land.movement_surf attribute
        import movement
        land.movement_surf = movement.create_movement_surface(land, params)   

   
    return land



def create_island_mask(dims, scape):
    mask = np.ma.masked_less_equal(scape,1e-8).mask
    return Landscape(dims,mask)



#function for reading in a pickled landscape stack
def load_pickled_scape_stack(filename):
    import cPickle
    with open(filename, 'rb') as f:
        land = cPickle.load(f)

    return land
