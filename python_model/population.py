#!/usr/bin/python 
#population.py

'''
##########################################

Module name:              population


Module contains:
                          - definition of the Population class
                          - function for creating a population of Individuals (incl. their genomes and associated data)
                          - associated functions


Author:                    Drew Ellison Hart
Email:                     drew.hart@berkeley.edu
Github:                    URL
Start date:                12-28-15
Documentation:             URL


##########################################
'''

import genome
import individual
import mating
import dispersal
import selection
import mutation
import landscape

import numpy as np
from numpy import random as r
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import cKDTree

import sys

#------------------------------------
# CLASSES ---------------------------
#------------------------------------


class Population:



    '''Population(self, N, individs, genomic_arch, size, birth_rate = None, death_rate = None, T = None)'''


     
    def __init__(self, N, individs, genomic_arch, size, T): 

        self.Nt = []                                        #list to record population size at each step (with starting size as first entry)

        self.N = None                                       #slot to hold an landscape.Landscape object of the current population density 

        self.K = None                                       #slot to hold an landscape.Landscape object of the local carrying capacity (i.e. 'target' dynamic equilibrium population density)

        self.individs = individs                            #dict of all individuals, as instances of individual.Individual class


        self.initial_size = len(self.individs)


        
        if size.__class__.__name__ in ['float', 'int'] and T <> None:
            self.size = [size] * T              
        elif size.__class__.__name__ in ['list']:
            assert T<> None and len(size) == T, "if expressing population size as a list, total model runtime must be provided and must equal the list length"
            self.size = size
        else:
            self.size = size



        self.genomic_arch = genomic_arch
        #Add other demographic parameters?? Other stuff in general??



        assert type(N) == int, "N must be an integer"
        assert type(self.individs) == dict, "self.individs must be a dictionary"
        assert list(set([i.__class__.__name__ for i in self.individs.values()])) == ['Individual'], "self.individs must be a dictionary containing only instances of the individual.Individual class"
        assert self.genomic_arch.__class__.__name__ == 'Genomic_Architecture', "self.genomic_arch must be an instance of the genome.Genomic_Architecture class"





    #####################
    ### OTHER METHODS ###
    #####################

    def census(self):
        return len(self.individs)


    #method to set self.K
    def set_K(self, K): #NOTE: Requires a landscape.Landscape instance
        self.K = K


    #method to set self.N
    def set_N(self, N): #NOTE: Requires a landscape.Landscape instance
        self.N = N


    #method to increment all population's age by one (also adds current pop size to tracking array)
    def birthday(self):
        #add current pop size to pop.N (for later demographic analysis)
        self.Nt.append(self.census())
        #increment age of all individuals
        [ind.birthday() for ind in self.individs.values()];



    #method to extract habitat values at each individual's coordinates for each land.scapes raster
    def query_habitat(self, land):
        [ind.query_habitat(land) for ind in self.individs.values()];
 




    #method to move all individuals simultaneously
    def move(self, land, params):
        [ind.move(land, params) for ind in self.individs.values()];

        self.query_habitat(land)



    #function for finding all the mating pairs in a population
    def find_mating_pairs(self, land, params):
        
        mating_pairs = mating.find_mates(self, params) 
        return(mating_pairs)




    #function for executing mating for a population
    def mate(self, land, params, mating_pairs):

        #pull necessary parameters from params dict
        mu_dispersal = params['mu_dispersal']
        sigma_dispersal = params['sigma_dispersal']
        sex = params['sex']
        repro_age = params['repro_age']

        
        num_offspring = 0


        for pair in mating_pairs:

            parent_centroid_x = np.mean((self.individs[pair[0]].x, self.individs[pair[1]].x))
            parent_centroid_y = np.mean((self.individs[pair[0]].y, self.individs[pair[1]].y))

            zygotes = mating.mate(self, pair, params)

            
            for zygote in zygotes:

                num_offspring += 1
                
                offspring_x, offspring_y = dispersal.disperse(land, parent_centroid_x, parent_centroid_y, mu_dispersal, sigma_dispersal)

                if sex == True:
                    offspring_sex = r.binomial(1, 0.5)

                age = 0

                offspring_key = max(self.individs.keys()) + 1

                if sex == True:
                    self.individs[offspring_key] = individual.Individual(zygote, offspring_x, offspring_y, offspring_sex, age)
                else:
                    self.individs[offspring_key] = individual.Individual(zygote, offspring_x, offspring_y, age)


        #sample all individuals' habitat values, to initiate for offspring
        self.query_habitat(land)

        print '\n\t%i individuals born' % num_offspring






    #method to carry out selection
    def select(self, t, params):
        selection.select(self, t, params, sigma_deaths = params['sigma_deaths'], density_dependent_fitness = params['density_dependent_fitness'])






    #method to carry out mutation
    def mutate(self, params, t):
        for ind in [ind for ind, individ in self.individs.items() if individ.age == 0]:
            mutation.mutate(self, self.individs[ind], self.genomic_arch, t, alpha_mut_s = params['alpha_mut_s'], beta_mut_s = params['beta_mut_s'])





    def check_extinct(self):
        if len(self.individs.keys()) == 0:
            print '\n\nYOUR POPULATION WENT EXTINCT!\n\n\t(Press <Enter> to exit.)'
            raw_input()
            sys.exit()
    



    def calc_density(self, land, window_width = None, normalize_by = 'none', min_0 = True, max_1 = False, max_val = None, set_N = False):

        '''
        Calculate an interpolated raster of local population density, using a window size of window_width.
        Valid values for normalize_by currently include 'census' and 'none'. If normalize_by = 'census', max_1 =
        True will cause the output density raster to vary between 0 and 1, rather than between 0 and the current
        max normalized density value. Window width will default to 1/10 of the larger raster dimension.
        '''

        #window width defaults to 1/10 the maximum landscape dimension
        if window_width == None:
            window_width = max(land.dims)*0.1

        #shorthand
        dims = land.dims

        #get a list of pop's coord-tuples
        c = self.get_coords().values() 

        #make window_width a float, to avoid Py2 integer-division issues
        window_width = float(window_width)
        
        #create meshgrid using window_width/2 as step size
        grid_j, grid_i = np.mgrid[0:dims[0]:complex("%ij" % (dims[0]/(window_width/2))), 0:dims[1]:complex("%ij" % (dims[1]/(window_width/2)))]

        #grid_j, grid_i = np.mgrid[0+(window_width/2):dims[0]-(window_width/2):complex("%ij" % (dims[0]/(window_width/2))), 0+(window_width/2):dims[1]-(window_width/2):complex("%ij" % (dims[1]/(window_width/2)))] 

        #flatten the arrays, so that I can run over them in a single for loop
        gj = grid_j.ravel()
        gi = grid_i.ravel()

        #make lists of tuples, of same length as gj, containing the window ll and ur coords
        window_ll = [(max(gj[n]-(window_width/2), 0), max(gi[n]-(window_width/2), 0)) for n in range(len(gj))]   #constrain min window vals to 0
        window_ur = [(min(gj[n]+(window_width/2), land.dims[0]), min(gi[n]+(window_width/2), land.dims[1])) for n in range(len(gj))] #constrain max window vals to each respective land dimension
        assert len(window_ll) == len(gj)
        assert len(window_ur) == len(gj)

        #make a list of the sizes of each window
        window_size = [(window_ur[n][0] - window_ll[n][0]) * (window_ur[n][1] - window_ll[n][1]) for n in range(len(gj))]#calculate size of this window (not always the same because of edge effects
        assert len(window_size) == len(gj)
        
        #make a list of the counts of all individs within each window
        window_ct = [len([ind for ind in range(len(c)) if (c[ind][0]>window_ll[n][0] and c[ind][0]<=window_ur[n][0]) and (c[ind][1]>window_ll[n][1] and c[ind][1]<=window_ur[n][1])]) for n in range(len(gj))] 
        assert len(window_ct) == len(gj)

        #divide window counts by window sizes
        window_dens = [window_ct[n]/window_size[n] for n in range(len(window_ct))] #divide by window size
        assert len(window_dens) == len(gj)

        #if normalize_by == census, then divide each density by total pop census size
        if normalize_by == 'census':
            N = self.census()
            window_dens = [dens/N for dens in window_dens]
        elif normalize_by == 'none':
            pass

        else:  #POTENTIALLY ADD OTHER OPTIONS HERE, i.e. to normalize by starting size?
            pass 

        #interpolate resulting density vals to a grid equal in size to the landscape
        new_gj, new_gi = np.mgrid[0:dims[0]-1:complex("%ij" % (dims[0])), 0:dims[1]-1:complex("%ij" % (dims[1]))]
        dens = interpolate.griddata(np.array(zip(list(gi), list(gj))), window_dens, (new_gj, new_gi), method = 'cubic')



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
            self.set_N(landscape.Landscape(dims, dens))
        
        else:
            return(landscape.Landscape(dims, dens))



    

    #function to discover loci above, below, or between threshold selection coefficients
    def find_loci(self, min_s=None, max_s = None):
        assert min_s<> None or max_s <> None, "No parameters provided. Must provide at least one value (either a max or min)."
        if max_s and not min_s:
            loci = dict([(i, np.array(range(len(self.genomic_arch.s[i])))[self.genomic_arch.s[i] <= max_s]) for i in range(len(self.genomic_arch.s))])
            #print "\n%i loci found\n" % sum([len(chrom_set) for chrom_set in loci.values()])
            return loci


        elif min_s and not max_s:
            loci = dict([(i, np.array(range(len(self.genomic_arch.s[i])))[self.genomic_arch.s[i] >= min_s]) for i in range(len(self.genomic_arch.s))])
            #print "\n%i loci found\n" % sum([len(chrom_set) for chrom_set in loci.values()])
            return loci
 

        else:
            loci = dict([(i, np.array(range(len(self.genomic_arch.s[i])))[np.array(self.genomic_arch.s[i] >= min_s) & np.array(self.genomic_arch.s[i] <= max_s)]) for i in range(len(self.genomic_arch.s))])
            #print "\n%i loci found\n" % sum([len(chrom_set) for chrom_set in loci.values()])
            return loci
 

   


   
 

    def get_habitat(self, individs = None):
        if individs <> None:
            return dict([(k, ind.habitat) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, ind.habitat) for k, ind in self.individs.items()])





    def get_age(self, individs = None):
        if individs <> None:
            return dict([(k, ind.age) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, ind.age) for k, ind in self.individs.items()])


    

    def get_genotype(self, chromosome, locus, format = 'mean', individs = None):

        dom_funcs = { 0 : np.mean,                          #codominance
                      1 : lambda x: np.ceil(np.mean(x)),    #dominant (with respect to relationship between allele 1 and habitat values -> 1
                      2 : lambda x: np.floor(np.mean(x))    #recessive (with respect to relationship between allele 1 and habitat values -> 1
                    }


        if individs == None:
            individs = self.individs.keys()
            #individs = range(len(self.genomic_arch.s[chromosome]))

        if format == 'biallelic':
            return dict([(i, self.individs[i].genome.genome[chromosome][locus, :]) for i in self.individs.keys() if i in individs]) 

        elif format == 'mean':
            return dict([(i, [dom_funcs[self.genomic_arch.dom[chromosome][locus]](self.individs[i].genome.genome[chromosome][locus,:]) ]) for i in self.individs.keys() if i in individs])


    def get_fitness(self):
        return selection.get_fitness(self)


    def hist_fitness(self):
        plt.hist(selection.get_fitness(self).values())



    def get_dom(self, chromosome, locus):
        return {locus: self.genomic_arch.dom[chromosome][locus]} 



    def get_env_var(self, chromosome, locus):
        return {locus: self.genomic_arch.env_var[chromosome][locus]} 





    def get_coords(self, individs = None):
        if individs <> None:
            return dict([(k, (float(ind.x), float(ind.y))) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, (float(ind.x), float(ind.y))) for k, ind in self.individs.items()])


    def get_x_coords(self, individs = None):
        if individs <> None:
            return dict([(k, (float(ind.x))) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, (float(ind.x))) for k, ind in self.individs.items()])


    def get_y_coords(self, individs = None):
        if individs <> None:
            return dict([(k, (float(ind.y))) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, (float(ind.y))) for k, ind in self.individs.items()])














    def show(self, land, scape_num = None, color = 'black', colorbar = True, markersize = 25, im_interp_method = 'nearest', alpha = False):
		#if land <> None:
		if scape_num <> None:
			land.scapes[scape_num].show(colorbar = colorbar, im_interp_method = im_interp_method, pop = True)
		else:
			land.show(colorbar = colorbar, im_interp_method = im_interp_method, pop = True)

		c = np.array(self.get_coords().values())
        #NOTE: subtract 0.5 to line up the points with the plt.imshow() grid of the land; imshow plots each pixel centered on its index, but the points then plot against those indices, so wind up shifted +0.5 in each axis
		x = c[:,0]-0.5
		y = c[:,1]-0.5
		if alpha == True:
			alpha = 0.6
		else:
			alpha = 1.0
	
		plt.scatter(x,y, s = markersize, c = color, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)
        #mpl.pyplot.plot([n[0] for n in coords], [n[1] for n in coords], 'ko', scalex = False, scaley = False, color = color, markersize = markersize, alpha = alpha)





    def show_individs(self, individs, land, scape_num = None, color = 'black', im_interp_method = 'nearest', markersize = 40, alpha = 0.5):
		#if land <> None and scape_num <> None:
		land.scapes[scape_num].show(im_interp_method = im_interp_method, pop = True)

		#coords = dict([(k, (ind.x, ind.y)) for k, ind in self.individs.items() if k in individs])
		c = np.array(self.get_coords(individs).values())
		#NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
		x = c[:,0]-0.5
		y = c[:,1]-0.5
		plt.scatter(x,y, s = markersize, c = color, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)
        #NOTE: perhaps worth figuring out how to label with the individual number!!
        #for k, coord_pair in coords.items():
            #ax = mpl.pyplot.plot(coord_pair[0], coord_pair[1], 'ko', scalex = False, scaley = False, color = color, markersize = 8.5)


    
    
    
    def show_density(self, land, window_width = None, normalize_by = 'census', max_1 = False, color = 'black'):
		dens = self.calc_density(land, window_width = window_width, normalize_by = normalize_by, max_1 = max_1)
		dens.show(im_interp_method = 'nearest', pop = True)
		
		c = np.array(self.get_coords().values())
		#NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
		x = c[:,0]-0.5
		y = c[:,1]-0.5
		plt.scatter(x,y, s = markersize, c = color, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)
        #ax = mpl.pyplot.plot([i[0] - 0.5 for i in c.values()], [i[1] - 0.5 for i in c.values()], 'ko', scalex = False, scaley = False, color = color, markersize = 8.5)

        #NOTE: perhaps worth figuring out how to label with the individual number!!







    def show_locus(self, chromosome, locus, land, scape_num = None, im_interp_method = 'nearest'): 
		

		if scape_num <> None:
			land.scapes[scape_num].show(im_interp_method = im_interp_method, pop = True)

		else:
			land.show(im_interp_method = im_interp_method, pop = True)
		
		genotypes = self.get_genotype(chromosome, locus) 

		colors = ['#3C22B4', '#80A6FF', '#FFFFFF'] # COLORS TO MATCH LANDSCAPE PALETTE EXTREMES, BUT WITH HYBRID A MIX OF THE EXTREMES RATHER THAN THE YELLOW AT THE MIDDLE OF THE PALETTE, FOR NICER VIEWING: blue = [0,0], light blue = [0,1], white = [1,1]
        #colors = ['#ff4d4d', '#ac72ac', '#4d4dff'] # red = [0,0], purple = [0,1], blue = [1,1] 

		for n, genotype in enumerate([0.0, 0.5, 1.0]):
			inds = [i for i, g in genotypes.items() if g[0] == genotype]
			c = np.array(self.get_coords(inds).values())
			x = c[:,0]-0.5
			y = c[:,1]-0.5
			plt.scatter(x,y, s = markersize, c = colors[n], alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)
    
            #mpl.pyplot.plot([coord[0] for coord in coords], [coord[1] for coord in coords], 'o', markersize = 11, scalex = False, scaley = False, color = colors[n], alpha = 0.8)



    #method for plotting a population pyramid
    #NOTE: NEED TO FIX THIS SO THAT EACH HALF PLOTS ON OPPOSITE SIDES OF THE Y-AXIS
    def show_pyramid(self):
        plt.hist([ind.age for ind in self.individs.values() if ind.sex == 0], orientation = 'horizontal', color = 'pink', alpha = 0.6)
        plt.hist([ind.age for ind in self.individs.values() if ind.sex == 1], orientation = 'horizontal', color = 'skyblue', alpha = 0.6)








    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)











#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def create_population(genomic_arch, land, params):



    #grab necessary params from params dict

    N = params['N']

    dims = params['dims']

    size = params['size']

    T = params['T']






    assert dims.__class__.__name__ in ['tuple', 'list'], "dims should be expressed as a tuple or a list"
    individs = dict()
    for i in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        individs[i] = individual.create_individual(genomic_arch, dims)

    pop = Population(N = N, individs = individs, genomic_arch = genomic_arch, size = size, T = T)
    

    #get initial habitat values
    pop.query_habitat(land)

    return pop






#function for reading in a pickled pop
def load_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)
    
    return pop


