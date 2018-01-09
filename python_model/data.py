#!/usr/bin/python
#data.py


'''
##########################################

Module name:              data

Module contents:          - definition of data class (i.e. a structured container for returning data)
                          - definition of data formatting functions
                          - definition of functions for sampling data, 
                            at specified frequencies and with specified arguments, 
                            according to the contents of the params['data'] section


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               01-01-18
Documentation:            URL


##########################################
'''


import numpy as np    
from numpy import random as r
import os, sys





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



        #create a Data.data object, where all of the collected data will be stored
        self.data = {'0':{
                         'genomic_arch':None,
                         'individs':None,
                         'land':None,
                         'params':None
                         }
                    }
     
        
        #grab the params['data'] content into a self.params attribute
        self.params = params['data']

        #and set the sampling times in it accordingly
        assert type(self.params['freq'] in (list, float, int))
        if type(self.params['freq']) == list:
            if 0 not in self.params['freq']:
                self.params['freq'] = [0] + self.params['freq']
            if params['T'] not in self.params['freq']:
                self.params['freq'] = self.params['freq'] + [params['T']]
            
            elif type(self.params['freq'] in (float, int)): 
                self.params['freq'] = range(0,params['T'], int(self.params['freq']))+params['T']



    #a collection method, to be called each timestep, which will collect needed data and then calls the functions to calculate them all and adds the results to self.stats
    def collect(self, pop):

        #if it's the appropriate time
        if pop.t in self.params['freq']:

            #sample data according to the scheme defined 
            data = sample_data(self.params['sampling_scheme'], **self.params['sampling_args'])

        self.data[t] = data

        if 'write_intermittent':
            self.write(self, drop_after_write = self.params['drop_after_write'])


        return()


        





    #a method to write data to the data_directory, then drop it from the Data object
    def write(self, drop_after_write = False): 
        
        #get data_directory_name
        data_dir = os.path.join(os.getcwd(), 'GEONOMICS_%s' % 'test')

        #get filename components
        rn = self.params['run_name']        #run_name
        i = 'ITERATION_NUMBER_HERE'         #NOTE: figure out some way for easily keeping track of iteration numbers
        t = max(self.data.keys())           #timestep
        gen_df = self['gen_data_format']    #gen data format
        geo_df = self['geo_data_format']    #geo data formats

        #get filenames
        gen_data_file = '%s_n%i_t%i_.%s' % (rn, i, t, self.extension_dict[self.params[gen_data_format]])
        geo_vect_data_file = '%s_n%i_t%i_.%s' % (rn, i, t, self.extension_dict[self.params[geo_data_format[0]]])
        if self.params['include_land']:
            geo_rast_data_file = '%s_n%i_t%i_.%s' % (rn, i, t, self.extension_dict[self.params[gen_data_format[1]]])


        #get gen_data, formatted as stipulated in params['data']
        gen_data, geo_data = format_data(self.params['gen_data_format'], self.params['geo_data_format'])

        #write gen_data file:
        with open(os.path.join(data_dir, gen_data_file), 'w') as f:
            f.write(gen_data)


        #write geo_data
        self.geo_data_write_fn_dict[geo_data[0]](geo_vect_data_file)
        
        if self.params['include_land'] == True:
            self.geo_data_write_fn_dict[geo_data[1]](geo_rast_data_file)
        
        if  drop_after_write == True:
            self.data[t] = None


        return()










#----------------------------------
# FUNCTIONS -----------------------
#----------------------------------


def sample_data(pop, sampling_scheme, n = None, points = None, radius = None, transect_endpoints = None, n_transect_points = None):

    '''<sampling_scheme> can be: 
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
    sample = []
    
    if sampling_scheme == 'all':
        sample = pop.individs.values()
    
    elif sampling_scheme == 'random':
        sample = list(r.choice(pop.individs.keys(), size = n, replace = False))

    elif sampling_scheme == 'points': 
        [sample.extend(point_sample(pop, n, p, radius)) for p in points]

    elif sampling_scheme == 'transect':

        transects = [np.linspace(ep[0], ep[1], n_transect_points) for ep in transect_endpoints] 
        [[sample.extend(point_sample(pop, n, p, radius)) for p in ts] for ts in transects]


    #make sure didn't somehow get duplicates
    sample = list(set(sample))

  

    #now sample data from those individuals

    #TODO WRITE SECOND HALF, WHICH POP.<methods> TO COLLECT DATA FOR ALL INDIVIDUALS


    
    return(sampled_data)






def point_sample(pop, n, point, radius):

#TODO: WRITE FUNCTION TO FIND <= n INDIVIDS WITHIN radius OF EACH POINT IN points, WHICH RETURNS A LIST OF INDIVID IDS

    pass




def format_data(gen_data_format, geo_data_format):

    '''<gen_data_format> can be:
                    'FASTA'
                    'VCF'
                    'ms'

       <geo_data_format> can be:
                1.) 'CSV'
                    'Shapefile'
                2.) 'Geotiff'
    '''

    #TODO: WRITE FUNCTION TO FORMAT DATA, USING THE SUBSIDIARY FUNCTIONS THAT I NEED TO WRITE BELOW





def format_fasta(gen_data):
    pass

def format_vcf(gen_data):
    pass

def format_ms(gen_data):
    pass






#TODO: use existing basic packages/functionalities to write these functions (geopandas; what for the raster?)
def write_shapefile(filename, data):
    pass

def write_csv(filename, data):
    pass

def write_geotiff(filename, data):
    pass
