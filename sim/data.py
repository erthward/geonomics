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


import numpy as np
from random import sample as rand_sample
from numpy import random as r
import os, sys
import datetime
import re
from shapely.geometry import Point
import pandas as pd
import geopandas as gpd
from osgeo import gdal


#TODO: 
    # - map this all out, and decide if it makes more sense to change Data to
    # Sampler class, and just have it sample the pops and land directly, rather
    # than copying all their objects to Data.data first; seems unnecessarily
    # complicated now

    # - consider if/how the emerging Changer/Sampler/___? standardization would
    # fit with the stats module. Is it reasonable to just keep the stats small
    # and simple, and to include their collection within the Sampler class?
    # Or should it should be a Collector or DataCollector class instead?

    # - then rewrite however makes the most sense (and that may actually be
    # simplest and easiest if I truly just rewrite the whole class, keeping and
    # implementing the data-formatting and -writing functions below it

    # - toy with and figure out the raster-writing functionality (i.e. need to
    # also collect/keep the PROJ and CRS of a raster I've read in, and also set
    # default EPSG codes of some sort???)

    # - decide if (and if so, how) any of this should go into the utils.io
    # module (at the least I think I could write a writer_raster() function
    # there and then just call that from here)

    # - finalize the initial go at the code (and the stats, while I'm at it),
    # then test and debug!


#------------------------------------
# CLASSES ---------------------------
#------------------------------------

class DataCollector:
    def __init__(self, params):

    #some lookup dicts for writing data 
        self.file_extension_dict =   {'VCF': 'vcf',
                            'FASTA': 'fasta',
                            'ms': 'ms',
                            'CSV': 'csv',
                            'Shapefile': 'shp',
                            'Geotiff': 'tif'
                            }

        self.write_geodata_fn_dict = {'CSV': write_csv,
                             'Shapefile': write_shapefile,
                             'GeoJSON': write_geojson,
                             'Geotiff': write_geotiff
                             }
 
        #grab the params['data'] content into a self.data_params attribute
        self.data_params = deepcopy(params.model.data)
        #grab the params['data'] content into a self.data_params attribute
        self.stats_params = deepcopy(params.model.data)




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
                             'GeoJSON': write_geojson,
                             'Geotiff': write_geotiff
                             }

        #pull necessary stuff from params
        T = params.model.time.T
        
        #grab the params['data'] content into a self.params attribute
        self.params = deepcopy(params.model.data)

        #create a Data.data object, where all of the collected data will be stored
        #NOTE: all data will be structured like this, to be easily parseable into various formats 
        self.data = {'0':{
                        'gen_arch':None,
                         'inds':None,
                         'land':None,
                         'params':None
                         }
                    }

        #check that frequency values are less than T
        assert type(self.params.freq in (list, float, int, type(None)))
        if type(self.params.freq) is list:
            assert ([n < T for n in self.params.freq]).all(), ('ERROR: Values '
            'provided for data-sampling frequency must be less than total model run-time.')
        elif type(self.params.freq) in (float, int):
            assert n < T, ('ERROR: Values provided for data-sampling frequency '
                           'must be less than total model run-time.')

        #and set the sampling times in it accordingly
        if type(self.params.freq) == list:
            if T-1 not in self.params.freq:
                self.params.freq = self.params.freq + [T-1]

        elif type(self.params.freq) in (float, int, type(None)):
            if self.params.freq in (0, None):
                self.params.freq = [T-1]
            else:
                self.params.freq = range(0, T, int(self.params.freq)) + [T-1]
                if T-1 not in self.params.freq:
                    self.params.freq = self.params.freq + [T-1]


    #a do_collection method, to be called each timestep, which will collect needed
    #data and then write the data (if write_intermittent == True)
    def do_collection(self, pop, land, params, drop_after_write = True):

        #if it's the appropriate time
        if pop.t in self.params.freq:

            #sample data according to the scheme defined 
            self.data[pop.t] = sample_data(pop, land, **self.params.sampling_args)

            if self.params.write_intermittent:
                self.write(pop, drop_after_write = drop_after_write)


    #a method to write data to the data_directory, then drop it from the Data object
    def write(self, model_name, pop, drop_after_write = True):

        #get data_directory_name
        data_dir = os.path.join(os.getcwd(), 'GEONOMICS_%s' % 'test')

        #create the data directory if it doesn't exist yet
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)

        it = pop.it                         #iteration number
        t = pop.t                           #timestep
        #t = max(self.data.keys())           #timestep
        gen_df = self.params.gen_data_format    #gen data format
        geo_df = self.params.geo_data_format    #geo data formats

        #get filenames
        gen_data_file = '%s_it%i_t%i_.%s' % (model_name, it, t, self.extension_dict[gen_df])
        geo_vect_data_file = '%s_it%i_t%i_.%s' % (model_name, it, t, self.extension_dict[geo_df[0]])
        if self.params.include_land:
            geo_rast_data_file = '%s_it%i_t%i_.%s' % (model_name, it, t, self.extension_dict[gen_df[1]])


        #TODO: get data


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

    sampled_data = {'gen_arch': pop.gen_arch,
                    'inds': sample,
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


def format_fasta(data):

    '''
    FASTA FORMAT:

    >idx:haploid_num|x_location|y_location|phenotype0;phenotype1;...;phenotypeN|env_var0;env_var1;...;env_varN|age|sex
    001110101010101010010101011101010110.....01011110

    '''
    row1 = '>%s:HAP;%s;%s;%s;%s;%s;%s\n'
    file_text = ''

    for ind in data['inds']:
        for hap in range(2):
            ind_row1 = re.sub('HAP', str(hap), row1)
            replace = tuple(map(lambda att: re.sub(',', '|', re.sub('[\[\] ]', '', str(getattr(ind, att)))),
                                ['idx', 'x', 'y', 'age', 'sex', 'phenotype', 'habitat']))
            ind_row1 = ind_row1 % replace
            ind_row2 = ''.join([str(base) for base in ind.genome[:,hap]]) + '\n'

            file_text = file_text + ind_row1 + ind_row2

    #with open(filename, 'w') as f:
    #    f.write(file_text)
    return(file_text)


def format_vcf(data, include_fixed_sites=False):

    #create a template header
        #NOTE: has 1 string slot for a date

        #TODO: DECIDE ON NECESSARY INFO AND FORMAT CONTENTS AND ADD METADATA ROWS HERE
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
    inds = sorted(data['inds'].keys())
    ind_cols = '\t'.join([str(i) for i in inds])
    cols = col_header_row % (ind_cols)

    #get a list of the chromosome numbers
    chroms = np.cumsum(data['gen_arch'].l_c)

    #and get all individuals' genomic data in a 'samplome' object (a 3-d array)
    samplome = np.array([data['inds'][i].genome for i in inds])
    
    #get loci of all segregating sites, if not_include_fixed_sites
    if not include_fixed_sites:
        #get segregating sites
        max_val = 2 * len(data['inds'])
        segs = np.where(samplome.sum(axis = 2).sum(axis = 0) > 0)[0]
        segs2 = np.where(samplome.sum(axis = 2).sum(axis = 0) < max_val)[0]
        loci = sorted(list(set(segs).intersection(set(segs2))))
    #or else get all loci
    else:
        loci = range(data['gen_arch'].L)

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


#TODO: WRITE THIS
def format_ms(gen_data):
    pass


#TODO: NEED TO READ IN ALL GEODATA FROM THE RASTER THAT WAS READ, IF
#APPLICABLE, OR ELSE SET DEFAULT PROJ AND CRS, AND THEN USE THOSE TO WRITE THIS OUT
def write_geopandas(filename, data, driver):
    attributes = ['idx', 'phenotype', 'habitat', 'age', 'sex'] 
    #FIXME: replace the call to str() below with something more sophisticated
    #that will actually separate distinct phenotype and habitat values into
    #separate, labeled columns
    df_dict = {att: [str(getattr(ind, att)) for ind in data['inds']] for att in attributes}
    pts = [Point(ind.x, ind.y) for ind in data['inds']]
    df_dict['pt'] = pts
    df = pd.DataFrame.from_dict(df_dict)
    gdf = gpd.GeoDataFrame(df, geometry = 'pt')
    if driver == 'CSV':
        gdf['x'] = gdf.pt.x
        gdf['y'] = gdf.pt.y
        gdf = gdf.drop(labels = ('pt'), axis = 1)
        gdf.to_csv(filename, index = False)
    else:
        gdf.to_file(filename, driver = driver)

def write_shapefile(filename, data):
    write_geopandas(filename, data, driver='ESRI Shapefile')

def write_geojson(filename, data):
    write_geopandas(filename, data, driver='GeoJSON')

def write_csv(filename, data):
    write_geopandas(filename, data, driver = 'CSV')


#TODO: NEED TO READ IN ALL THE GEODATA FROM THE RASTER THAT WAS READ, IF
#APPLICABLE, OR ELSE SET DEFAULT PROJ AND CRS, AND THEN USE THOSE TO WRITE THIS OUT
def write_geotiff(filename, data):

    #TODO: tweak code taken from https://gis.stackexchange.com/questions/58517/python-gdal-save-array-as-raster-with-projection-from-other-file

    #get the driver
    driver = gdal.GetDriverByName('GTiff')

    #get values
    x_pixels = data.dim[0]  # number of pixels in x
    y_pixels = data.dim[1]  # number of pixels in y
    PIXEL_SIZE = data.res[0]  # size of the pixel...        
    x_min = data.ulc[0]
    y_max = data.ulc[1] + (data.dim[1]*data.res[1]) # x_min & y_max are like the "top left" corner.
    wkt_projection = 'a projection in wkt that you got from other file'

    dataset = driver.Create(
       filename+'.tif',
       x_pixels,
       y_pixels,
       1,
       gdal.GDT_Float32, )

    dataset.SetGeoTransform((
       x_min,      # 0
       PIXEL_SIZE, # 1
       0,          # 2
       y_max,      # 3
       0,          # 4
       -PIXEL_SIZE))

    dataset.SetProjection(wkt_projection)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()  # Write to disk.

