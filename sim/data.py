#!/usr/bin/python
#data.py


'''
##########################################

Module name:              sim/data

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

#geonmics imports
from utils import io

#other imports
import numpy as np
from random import sample as rand_sample
from numpy import random as r
import os, sys
import datetime
import re
from shapely.geometry import Point, MultiPolygon
import pandas as pd
import geopandas as gpd
from osgeo import gdal
from itertools import chain


#------------------------------------
# CLASSES ---------------------------
#------------------------------------

class DataSampler:
    def __init__(self, model_name, T, params):

    #some lookup dicts for writing data 
        self.file_extension_dict =   {'VCF': 'vcf',
                            'FASTA': 'fasta',
                            'CSV': 'csv',
                            'Shapefile': 'shp',
                            'GeoJSON': 'json',
                            'Geotiff': 'tif'
                            }

        self.write_geodata_fn_dict = {'CSV': io.write_csv,
                             'Shapefile': io.write_shapefile,
                             'GeoJSON': io.write_geojson,
                             'Geotiff': io.write_geotiff,
                             'txt': io.write_array
                            }

        #set other attributes
        self.model_name = model_name
        self.T = T

        #grab the params['data'] contents into objects
        sampling_params = params.model.data.sampling
        format_params = params.model.data.format

        #TODO: decide what to do about the stats stuff!
        #grab the params['data'] content into a self.data_params attribute
        self.stats_params = deepcopy(params.model.data)

        #get the sampling scheme
        self.scheme = sampling_params.scheme

        #and run some asserts
        assert self.scheme in ['all', 'random',
            'points', 'transect'], ("The sampling scheme provided in the "
            "parameters must be one of the following values: 'all', 'random', "
            "'points', or 'transect'.")

        if sampling_params.scheme != 'all':
            assert 'n' in sampling_params.keys(), ("If the "
            "sampling scheme is not 'all' then the 'n' parameter must be "
            "defined, indicating the number of individuals to be sampled "
            "each time data is collected.")
            assert type(sampling_params.n) is int, ("The "
            "'n' data-sampling parameter must be an integer.")
        #TODO: Add more assert statements here to check that only the right
        #combinations of parameters can be provided

        #get the number of individuals to sample, if applicable
        self.n = None
        if sampling_params.scheme != 'all':
            self.n = sampling_params.n

        #calculate the transect points, if 'transect' is the chosen sampling
        #scheme, or just get the points if 'points' is the scheme
        self.pts = None
        if sampling_params.scheme == 'points':
            self.pts = sampling_params.points
        if sampling_params.scheme == 'transect':
            endpts = sampling_params.transect_endpoints
            n_transect_pts = sampling_params.n_transect_points
            pts = get_transect_points(endpoints = endpts, n = n_transect_pts)
            self.pts = chain.from_iterable(pts)

        #create the point buffers, if either 'transect' or 'points' is the 
        #chosen sampling scheme
        self.pts_buff = None
        if sampling_params.scheme in ['transect', 'points']:
            self.pts_buff = make_point_buffers(self.pts,
                                        sampling_params.radius)

        #get the 'include_land' param (defaults False)
        self.include_land = False
        if ('include_land' in sampling_params.keys() and
            type(sampling_params.include_land) is bool):
            self.include_land = sampling_params.include_land

        #get the 'include_fixed_sites' param (defaults False)
        self.include_fixed_sites = False
        if ('include_fixed_sites' in sampling_params.keys() and
            type(sampling_params.include_land) is bool):
            self.include_fixed_sites = sampling_params.include_fixed_sites

        #get the 'when' parameter
        self.when = sampling_params.when

        #check type- and value-validity of self.when, and update its value
        #as necessary
        assert type(self.when in (list, float, int, type(None)))
        #if it's a list, make sure no values are greater than final timestep
        if type(self.when) is list:
            assert ([n < self.T for n in self.when]).all(), ('ERROR:'
            ' Values provided for sampling times must be less '
            'than total model run-time.')
            #add the last timestep, if not already included
            if self.T-1 not in self.when:
                self.when = self.when + [self.T-1]
        #if it's a float, int, or None
        elif type(self.when) in (float, int, type(None)):
            #check value is less than or equal to last timestep (or None)
            assert self.when is None or self.when < self.T, ('ERROR: Values '
            'provided for sampling times must be less than total '
            'model run-time.')
            #make it a list containing just last timestep, if 0 or None
            if self.when in (0, None):
                self.when = [T-1]
            #make it a stepwise timestep list, if integer other than 0
            else:
                self.when = [*range(0, T, int(self.when))]
                if self.T-1 not in self.when:
                    self.when.append(self.T -1)

        #now turn the when attribute into an iterator
        self.when = iter(self.when)
        #and grab the next timestep into next_t
        self.next_t = None
        self.set_next_t()

        #grab the genetic data formats as a separate attribute
        self.gen_formats = format_params.gen_format
        #change the gen_formats attribute to a list if it came in as a string
        if type(self.gen_formats) == str:
            self.gen_formats = [self.gen_formats]
        #also grab the geographic data formats as a separate attribute
        self.geo_formats = [format_params.geo_vect_format]
        #and grab the raster format, if a raster is required
        #NOTE: added to a separate attribute because this is written per
        #timestep, not per population within timestep
        self.rast_format = None
        if sampling_params.include_land and 'geo_rast_format' in format_params.keys():
            self.rast_format = format_params.geo_rast_format


    #method to set self.next_t
    def set_next_t(self):
        self.next_t = next(self.when)


    #method to create filenames for genetic and geographic datafiles
    def make_filenames(self, iteration, pop_name):
        filenames = []
        for att_name in ['gen_formats', 'geo_formats']:
            filenames.append(['m%s_i%i_t%i_p%s.%s' % (self.model_name, iteration,
                        self.next_t, pop_name, self.extension_dict[fmt])
                        for fmt in getattr(self, att_name)])
        return(filenames)


    #a method to be called each timestep, which will collect needed
    #data and then write the data (if write_intermittent == True) if it's
    #the right timestep
    def do_collection(self, community, iteration):
        #tracker to determine whether or not to update self.next_t
        update_next_t = False

        #for each population
        for pop in community.values():

            #if this timestep is scheduled for sampling
            if pop.t == self.next_t:

                #get the data directory name for this timestep
                dirname = os.path.join(os.getcwd(),
                   'gnx_m%s_i%i_t%i' % (self.model_name,
                    iteration, self.next_t))

                #get the subdirectory for this population
                subdirname = os.path.join(dirname, 'p', pop.name)

                #and create (and its parent data directory, if needed)
                os.makedirs(subdirname)

                #get filenames
                gen_files, geo_files = make_filenames(iteration = iteration,
                                                        pop_name = pop.name)

                #sample data according to the scheme defined 
                individs = self.get_individuals(pop)
                #for each genetic data format to be written
                for n, data_format in enumerate(self.gen_formats):

                    #format the data accordingly
                    data = format_gen_data(data_format = data_format,
                                        individs = individs, pop = pop)

                    #then write it to disk
                    self.write_gendata(subdirname, gen_files[n])

                #also write the geodata for this pop
                for n, data_format in enumerate(self.geo_formats):
                    #write the geodata to disk
                    self.write_geodata_fn_dict[data_format](dirname =subdirname,
                                    filename = geo_files[n], individuals = individs)

                #set the update_next_t tracker to True
                update_next_t = True

        #write the raster, if necessary
        if self.rast_format is not None:
            #for each Scape
            for scape in community.land.values():
                #get the raster filename
                filename = 'm%s_i%i_t%i_s%s.%s' % (self.model_name, iteration,
                                    self.next_t, scape.name, 
                                    self.file_extension_dict[self.rast_format])
                #and write it to disk
                self.write_geodata_fn_dict[self.rast_format](dirname, filename, scape)

        #update self.next_t, if required
        if update_next_t:
            self.set_next_t()


    def do_random_sample(self, individs):
        inds = r.choice(individs, size = self.n, replace = False)
        return(inds)


    def do_point_sample(self, pop):
        #TODO: see if this needs to be sped up any more
        inds = [i for i,v in pop.items() if self.pts_buff.contains(Point(v.x, v.y))]
        if len(inds) > self.n:
            inds = do_random_sample(individs = inds)
        return(inds)


    def get_individuals(self, pop):

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

        #get a set of indices for the individuals in the sample
        sample = set()
        #take the whole population, if scheme is 'all'
        if self.scheme == 'all':
            sample.update([*pop])
        #or take a random sample, if scheme is 'random'
        elif self.scheme == 'random':
            inds = do_random_sample([*pop])
            sample.update(inds)
        #or take individuals within a given radius of self.buffs (which could
        #have come from a set of input points or from a calculated set of
        #transect points), if scheme is 'points' or 'transect'
        elif scheme in ['points', 'transect']:
            inds = do_point_sample(pop)
            sample.update(inds)
        #convert sample to a dict of individuals
        sample = {i:pop[i] for i in sample}
        return(sample)


    def format_gen_data(self, data_format, individs, pop):
        '''<data_format> can be:
                            'FASTA'
                            'VCF'
        '''
        if data_format == 'fasta':
            formatted_data = format_fasta(individs)
        elif data_format == 'vcf':
            formatted_data = format_vcf(individs, pop.gen_arch,
                            include_fixed_sites = self.include_fixed_sites)
        return(formatted_data)


    def write_gendata(self, dirname, filename, gen_data):
        io.write_file(dirname, filename, gen_data)


    def write_geodata(self, dirname, filename, data):
        ext = os.splitext(filename)[1]
        write_fn = self.write_geodata_fn_dict[ext]
        write_fn(filename, data)

#----------------------------------
# FUNCTIONS -----------------------
#----------------------------------

#a function to get a set of n evenly spaced points between endpoints
def get_transect_points(endpoints, n):
    x_pts = np.linspace(endpoints[0][0] , endpoints[1][0], n)
    y_pts = np.linspace(endpoints[0][1] , endpoints[1][1], n)
    return(list(zip(x_pts, y_pts)))


#a function to make shapely geometry buffers around a set of points
def make_point_buffers(points, radius):
    pts = [Point(p[0], p[1]) for p in points]
    buffs = [p.buffer(radius) for p in pts]
    buff_poly = MultiPolygon(buffs)
    return(buff_poly)


def format_fasta(individs):

    '''
    FASTA FORMAT:

    >idx:haploid_num|x_location|y_location|phenotype0;phenotype1;...;phenotypeN|env_var0;env_var1;...;env_varN|age|sex
    001110101010101010010101011101010110.....01011110

    '''
    row1 = '>%s:HAP;%s;%s;%s;%s;%s;%s\n'
    file_text = ''

    for ind in individs:
        for hap in range(2):
            ind_row1 = re.sub('HAP', str(hap), row1)
            replace = tuple(map(lambda att: re.sub(',', '|', re.sub('[\[\] ]', '', str(getattr(ind, att)))),
                                ['idx', 'x', 'y', 'age', 'sex', 'phenotype', 'habitat']))
            ind_row1 = ind_row1 % replace
            ind_row2 = ''.join([str(base) for base in ind.genome[:,hap]]) + '\n'

            file_text = file_text + ind_row1 + ind_row2

    return(file_text)


def format_vcf(individs, gen_arch, include_fixed_sites=False):

    #create a template header
        #NOTE: has 1 string slot for a date

        #TODO: DECIDE ON NECESSARY INFO AND FORMAT CONTENTS, THEN ADD METADATA ROWS HERE
    header = '''##fileformat=VCFv4.2
##fileDate=%s
##source=Geonomics
'''

    #template column-header row
    col_header_row = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n'
        #NOTE: this has 1 string slot for a tab-separated list of all individ ids

    #template data row
    #TODO: UPDATE/CHANGE THE INFO AND FORMAT PORTIONS OF THIS TEMPLATE, AFTER I DECIDE ON THEIR CONTENTS (above)
    data_row = '%i\t%i\t.\tA\tT\t1000\tPASS\tMID=216754706;S=0;DOM=0.5;PO=1;GO=118200;MT=1;DP=1000\tGT\t%s\n'
        #NOTE: this has 2 integer slots, then 1 string slot for:
            #- chrom number (NOTE: unpythonically, starts from 1)
            #- locus number (NOTE: reported cumulative from locus 0, not from start of each chrom)
            #- a tab-separated list of individs' genotypes at this locus

    #create a col_header_row for this data
    inds = sorted(individs.keys())
    ind_cols = '\t'.join([str(i) for i in inds])
    cols = col_header_row % (ind_cols)

    #get a list of the chromosome numbers
    chroms = np.cumsum(gen_arch.l_c)

    #and get all individuals' genomic data in a 'samplome' object (a 3-d array)
    samplome = np.array([individs[i].genome for i in inds])

    #get loci of all segregating sites, if not_include_fixed_sites
    if not include_fixed_sites:
        #get segregating sites
        max_val = 2 * len(individs)
        segs = np.where(samplome.sum(axis = 2).sum(axis = 0) > 0)[0]
        segs2 = np.where(samplome.sum(axis = 2).sum(axis = 0) < max_val)[0]
        loci = sorted(list(set(segs).intersection(set(segs2))))
    #or else get all loci
    else:
        loci = range(gen_arch.L)

    #and get the sites' chrom nums
    chroms = [list((locus - chroms) < 0).index(True) for locus in loci]

    #build all the VCF data rows
    rows = ''
    for n, locus in enumerate(loci):
        genotypes = samplome[:,locus,:]
        genotypes = '\t'.join(['%i|%i' % (genotypes[i,0], genotypes[i,1]) for i in range(np.shape(genotypes)[0])])

        rows = rows + data_row % (chroms[n], locus, genotypes)

    #get the date
    now = datetime.datetime.now()
    month = str(now.month).zfill(2)
    day = str(now.day).zfill(2)
    date = '%d%s%s' % (now.year, month, day)

    #paste all the VCF content together
    out_vcf = ''.join([header % date, cols, rows])

    #return it
    return(out_vcf)

