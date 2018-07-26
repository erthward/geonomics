#!/usr/bin/python
#data.py


'''
##########################################

Module name:              data

Module contents:          - definition of data class (i.e. a structured container for returning data)
                          - definition of data formatting functions
                          - definition of functions for sampling data, 
                            at specified frequencies and with specified arguments, 
                            according to the contents of the params['model']['data'] section


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               01-01-18
Documentation:            URL


##########################################
'''


import numpy as np    
from random import sample as rand_sample
from numpy import random as r
import os, sys
import datetime





#------------------------------------
# CLASSES ---------------------------
#------------------------------------

class Data: 
    def __init__(self, params): 
        
        #some lookup dicts for writing data 
        self.extension_dict =   {'VCF': 'vcf',
                            'FASTA': 'fasta',
                            'ms': 'ms',
                            'CSV': 'csv',
                            'Shapefile': 'shp',
                            'Geotiff': 'tif'
                            }
       
        self.geo_data_write_fn_dict = {'CSV': write_csv,
                             'Shapefile': write_shapefile,
                             'Geotiff': write_geotiff
                             }

        #pull necessary stuff from params
        T = params.model.its.main.T
        data_params = params.model.data

        #create a Data.data object, where all of the collected data will be stored
        self.data = {'0':{
                        'genomic_arch':None,   #NOTE: all data will be structured like this, to be easily parseable into various formats
                         'individs':None,
                         'land':None,
                         'params':None
                         }
                    }
     
        
        #grab the params['data'] content into a self.params attribute
        self.params = data_params

        #and set the sampling times in it accordingly
        assert type(self.params.freq in (list, float, int))

        if type(self.params.freq) == list:
            if 0 not in self.params.freq:
                self.params.freq = [0] + self.params.freq
            if T not in self.params.freq:
                self.params.freq = self.params.freq + [T]
            
        elif type(self.params.freq) in (float, int): 
            self.params.freq = range(0, T, int(self.params.freq)) + [T]





    #a do_collection method, to be called each timestep, which will collect needed data and then calls the functions to calculate them all and adds the results to self.stats
    def do_collection(self, pop, land, params, drop_after_write = True):

        #if it's the appropriate time
        if pop.t in self.params.freq:

            #sample data according to the scheme defined 
            self.data[pop.t] = sample_data(pop, land, **self.params.sampling_args)

            if self.params.write_intermittent:
                self.write(pop, drop_after_write = drop_after_write)




        





    #a method to write data to the data_directory, then drop it from the Data object
    def write(self, pop, drop_after_write = True): 
        
        #get data_directory_name
        data_dir = os.path.join(os.getcwd(), 'GEONOMICS_%s' % 'test')

        #create the data directory if it doesn't exist yet
        if not os.path.exists(data_dir):
                os.makedirs(data_dir)


        #get filename components
        rn = self.params.run_name        #run_name (could be used to identify separate paramterizations
                                            #and/or separate populations/species modeled simultaneously)
        it = pop.it                         #iteration number
        t = pop.t                           #timestep
        #t = max(self.data.keys())           #timestep
        gen_df = self.params.gen_data_format    #gen data format
        geo_df = self.params.geo_data_format    #geo data formats

        #get filenames
        gen_data_file = '%s_it%i_t%i_.%s' % (rn, it, t, self.extension_dict[gen_df])
        geo_vect_data_file = '%s_it%i_t%i_.%s' % (rn, it, t, self.extension_dict[geo_df[0]])
        if self.params.include_land:
            geo_rast_data_file = '%s_it%i_t%i_.%s' % (rn, it, t, self.extension_dict[gen_df[1]])


        #get data
        
        
        #format the data as stipulated in params.model.data
        gen_data = format_gen_data(self.data[t], self.params.gen_data_format)

        geo_data = format_geo_data(self.data, self.params.geo_data_format)



        #write gen_data file:
        with open(os.path.join(data_dir, gen_data_file), 'w') as f:
            f.write(gen_data)


        #TODO: GET RID OF THIS if STATEMENT AND MOVE BACK OUT ONE INDENT ONCE GEO-DATA WRITING WORKING
        if True == False:
            #write geo_data
            self.geo_data_write_fn_dict[geo_data[0]](geo_vect_data_file)
        
            if self.params.include_land == True:
                self.geo_data_write_fn_dict[geo_data[1]](geo_rast_data_file)
            
            if  drop_after_write == True:
                self.data[t] = None


        return()










#----------------------------------
# FUNCTIONS -----------------------
#----------------------------------


def sample_data(pop, land, scheme, n = None, points = None, radius = None, transect_endpoints = None, n_transect_points = None):

    '''<scheme> can be: 
                    'all'       --> takes whole population
                    'random'    --> takes random sample of size <n>
                                    from anywhere on landscape
                    'points'    --> takes random sample of size <n> 
                                    within <radius> (in units of 
                                    cells) of each point
                    'transect'  --> takes random sample of size <n> 
                                    within <radius> at 
                                    <n_transect_points> evenly spaced
                                    points along a transect 
                                    between each tuple of endpoints 
                                    listed in <transect_endpoints>

    '''

    #get a dict of the individuals to be sampled
    sample = {}
    
    if scheme == 'all':
        sample.update(pop.items())
    
    elif scheme == 'random':
        ids = r.choice([*pop], size = n, replace = False)
        sample.update({i:v for i,v in pop.items() if i in ids})

    elif scheme == 'points': 
        sample.update(point_sample(pop, n, points, radius))

    elif scheme == 'transect':
        transect_pts = [get_transect_points(eps, n_transect_points) for eps in transect_endpoints]
        [sample.update(point_sample(pop, n, points, radius)) for points in transect_pts]



    #now sample data from those individuals

    sampled_data = {'genomic_arch': pop.genomic_arch,
                    'individs': sample,
                    'land': land
                    }


    
    return(sampled_data)





def get_transect_points(eps, n_pts):
    x_pts = np.linspace(eps[0][0] , eps[1][0], n_pts)
    y_pts = np.linspace(eps[0][1] , eps[1][1], n_pts)
    return(list(zip(x_pts, y_pts)))
     




def point_sample(pop, n, points, radius):
    from shapely import geometry
   
    sample = {}

    #create lists of points and then buffers
    pts = [geometry.Point(p[0], p[1]) for p in points]
    buffs = [p.buffer(radius) for p in pts]

    #NOTE: I'm sure there's a faster way than the below (and one that would also assure that 
    #the same individual isn't sampled at 2 points that overlap), but for now not worrying 
    #about better code, not mission-critical
    
    for b in buffs:
        inds = {i:v for i,v in pop.items() if b.contains(geometry.Point(v.x, v.y))}
        if len(inds) > n:
            point_sample = rand_sample(list(inds.keys()), n)
            sample.update({i:v for i,v in inds.items() if i in point_sample}) 
        else:
            sample.update(inds)


    return(sample)

        









def format_gen_data(data, data_format):

    '''<data_format> can be:
                    'FASTA'
                    'VCF'
                    'ms'
'''


    format_fns = {'FASTA': format_fasta,
                  'VCF': format_vcf,
                  'ms': format_ms
                  }

    formatted_data = format_fns[data_format](data)
    

    return(formatted_data)





def format_geo_data(data, format):

    
    '''
        <geo_data_format> can be:
                1.) 'CSV'
                    'Shapefile'
                2.) 'Geotiff'
    '''
pass



#TODO: WRITE THIS
def format_fasta(gen_data):
    pass




def format_vcf(data):


    #template header
        #NOTE: HAS 1 STRING SLOT FOR:
            #date

        #TODO: TOOK THIS HEADER FROM A RANDOM SLiM FILE I HAD ON HAND; ADJUST AND/OR REMOVE/ADD STUFF FROM IT AS NEEDED
    header = '''##fileformat=VCFv4.2
##fileDate=%s
##source=Geonomics
##INFO=<ID=MID,Number=1,Type=Integer,Description="Mutation ID in SLiM">
##INFO=<ID=S,Number=1,Type=Float,Description="Selection Coefficient">
##INFO=<ID=DOM,Number=1,Type=Float,Description="Dominance">
##INFO=<ID=PO,Number=1,Type=Integer,Description="Population of Origin">
##INFO=<ID=GO,Number=1,Type=Integer,Description="Generation of Origin">
##INFO=<ID=MT,Number=1,Type=Integer,Description="Mutation Type">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description="Multiallelic">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'''

    #template column-header row
    col_header_row = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n'
        #NOTE: HAS 1 STRING SLOT FOR:
                #tab-separated list of all individ ids

    #and template data row
    data_row = '%i\t%i\t.\tA\tT\t1000\tPASS\tMID=216754706;S=0;DOM=0.5;PO=1;GO=118200;MT=1;DP=1000\tGT\t%s\n'
        #NOTE: HAS 2 INTEGER SLOTS, THEN 1 STRING SLOT FOR:
            #chrom number (NOTE: unpythonically, starts from 1)
            #locus number (NOTE: reported cumulative from locus 0, not from start of each chrom)
            #tab-separated list of individs' genotypes at this locus

            #TODO: UPDATE/CHANGE INFO PORTION OF THIS TEMPLATE



    #create col_header_row for this data
    inds = sorted(data['individs'].keys())
    ind_cols = '\t'.join([str(i) for i in inds])
    cols = col_header_row % (ind_cols)


    #get chroms and loci of all segregating sites
    samplome = np.array([data['individs'][i].genome for i in inds])
    #get segregating sites
    max_val = 2 * len(data['individs'])
    segs = np.where(samplome.sum(axis = 2).sum(axis = 0) > 0)[0]
    segs2 = np.where(samplome.sum(axis = 2).sum(axis = 0) < max_val)[0]
    segs = sorted(list(set(segs).intersection(set(segs2))))
    #and get their chrom nums
    chroms = np.cumsum(data['genomic_arch'].l_c)
    chroms = [list((site - chroms) < 0).index(True) for site in segs]

    rows = ''
    for n, site in enumerate(segs):
        genotypes = samplome[:,site,:]
        genotypes = '\t'.join(['%i|%i' % (genotypes[i,0], genotypes[i,1]) for i in range(np.shape(genotypes)[0])])

        rows = rows + data_row % (chroms[n], site, genotypes) 
       


    #and get the date
    now = datetime.datetime.now()
    month = str(now.month).zfill(2) 
    day = str(now.day).zfill(2) 
    date = '%d%s%s' % (now.year, month, day)


    #paste all the VCF content together
    out_vcf = ''.join([header % date, cols, rows])

    #and return it
    return(out_vcf)





#TODO: WRITE THIS
def format_ms(gen_data):
    pass






#TODO: use existing basic packages/functionalities to write these functions (geopandas; what for the raster?)
def write_shapefile(filename, data):
    pass

def write_csv(filename, data):
    pass

def write_geotiff(filename, data):
    pass
