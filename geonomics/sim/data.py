#!/usr/bin/python
#data.py


'''
Classes and functions for implementing data collection and output
'''

#geonmics imports
from geonomics.utils.io import (_write_csv, _write_shapefile, _write_geojson,
                                _write_file)

#other imports
import numpy as np
from random import sample as rand_sample
from numpy import random as r
import os, sys
import datetime
import re
from shapely.geometry import Point
from shapely.ops import cascaded_union
import pandas as pd
import geopandas as gpd
from itertools import chain


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#a DataCollector class, to paramtereize and manage data 
#collection, and write data to disk
class _DataCollector:
    def __init__(self, model_name, params):

    #some lookup dicts for writing data 
        self.file_extension_dict =   {'vcf': 'vcf',
                            'fasta': 'fasta',
                            'csv': 'csv',
                            'shapefile': 'shp',
                            'geojson': 'json',
                            'geotiff': 'tif'
                            }

        self.write_geodata_fn_dict = {'csv': _write_csv,
                            'shapefile': _write_shapefile,
                            'geojson': _write_geojson,
                            }

        #set model name and T
        self.model_name = model_name
        self.T = params.model.T

        #grab the params['data'] contents into objects
        sampling_params = params.model.data.sampling
        format_params = params.model.data.format

        #get the sampling scheme
        self.scheme = sampling_params.scheme

        #and run some asserts
        assert self.scheme in ['all', 'random',
            'point', 'transect'], ("The sampling scheme provided in the "
            "parameters must be one of the following values: 'all', 'random', "
            "'point', or 'transect'.")

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
        #scheme, or just get the points if 'point' is the scheme
        self.pts = None
        if sampling_params.scheme == 'point':
            self.pts = sampling_params.points
        if sampling_params.scheme == 'transect':
            endpts = sampling_params.transect_endpoints
            n_transect_pts = sampling_params.n_transect_points
            pts = _get_transect_points(endpoints = endpts, n = n_transect_pts)
            self.pts = pts

        #create the point buffers, if either 'transect' or 'point' is the 
        #chosen sampling scheme
        self.pts_area = None
        if sampling_params.scheme in ['transect', 'point']:
            self.pts_area = _make_point_buffers(self.pts,
                                        sampling_params.radius)

        #get the 'include_landscape' param (defaults False)
        self.include_landscape = False
        if ('include_landscape' in sampling_params.keys() and
            type(sampling_params.include_landscape) is bool):
            self.include_landscape = sampling_params.include_landscape

        #get the 'include_fixed_sites' param (defaults False)
        self.include_fixed_sites = False
        if ('include_fixed_sites' in sampling_params.keys() and
            type(sampling_params.include_landscape) is bool):
            self.include_fixed_sites = sampling_params.include_fixed_sites

        #set the when attribute to the 'when' parameter
        self.when = sampling_params.when
        #check type- and value-validity of self.when, and update its value
        #as necessary
        assert type(self.when in (list, float, int, type(None)))
        #if it's a list, make sure no values are greater than final timestep
        if type(self.when) is list:
            assert np.array([n < self.T for n in self.when]).all(), (' Values '
                'provided for sampling times must be less than total '
                'model run-time.')
        #if it's a float, int, or None
        elif type(self.when) in (float, int, type(None)):
            #check value is less than or equal to last timestep (or None)
            assert self.when is None or self.when < self.T, ('Values '
            'provided for sampling times must be less than total '
            'model run-time.')
            #make it a list containing nothing (last timestep will be added),
            #if 0 or None
            if self.when in (0, None):
                self.when = []
            #make it a stepwise timestep list, if integer other than 0
            else:
                self.when = [*range(0, self.T, int(self.when))]
        #add the final timestep, if not already in self.when
        if (len(self.when) == 0 or self.when[-1] != self.T-1):
            self.when.append(self.T-1)
        #now turn the when attribute into an iterator
        self.when = iter(self.when)
        #and set next_t
        self._set_next_t()

        #grab the genetic data formats as a separate attribute
        self.gen_formats = format_params.gen_format
        #change the gen_formats attribute to a list if it came in as a string
        if type(self.gen_formats) == str:
            self.gen_formats = [self.gen_formats]
        #also grab the geographic data formats as a separate attribute
        self.geo_formats = [format_params.geo_vect_format]
        #and grab the raster format, if a raster is required
        #NOTE: added to a separate attribute because this is written per
        #timestep, not per species within timestep
        self.rast_format = None
        if (sampling_params.include_landscape
            and 'geo_rast_format' in format_params.keys()):
            self.rast_format = format_params.geo_rast_format

        # grab the format for CSV files containing the non-neutral loci
        # NOTE: None indicates that these files should not be produced
        self.nonneut_loc_format = format_params.nonneut_loc_format
        assert self.nonneut_loc_format in ['csv', None], ("the "
            "'nonneut_loc_format' parameter must be either 'csv' or None.")

    #method to set self.next_t
    def _set_next_t(self):
        try:
            self.next_t = next(self.when)
        #if there are no further timestep values left then this is the last
        #timestep, so set self.next_t to None
        except StopIteration:
            assert self.next_t == self.T -1, ("Model._set_next_t() threw a "
            "StopIteration error, but the current value of Model.next_t is "
            "not the final timestep (instead, it is %i).\n\n") % self.next_t
            self.next_t = None

    #method to create filenames for genetic and geographic datafiles
    def _make_filenames(self, iteration, spp_name):
        filenames = []
        for att_name in ['gen_formats', 'geo_formats']:
            filenames.append(['mod-%s_it-%i_t-%i_spp-%s.%s' % (self.model_name,
              iteration, self.next_t, spp_name, self.file_extension_dict[fmt])
                        for fmt in getattr(self, att_name)])
        return(filenames)

    #a method to be called each timestep, which will collect needed
    #data and then write the data (if write_intermittent == True) if it's
    #the right timestep
    #TODO: CONSIDER NOMENCLATURE CHANGE HERE AND IN CLASS NAME!
    def _write_data(self, community, land, iteration):

        #if this timestep is scheduled for sampling
        if community.t == self.next_t:

            #for each species
            for spp in community.values():

                #TODO: Probably get rid of this conditional; should be
                #unnecessary, and is not a direct check of the in-sync
                #assumption at any rate
                #double-check that each species' timestep is in sync with
                #comm.t and is scheduled for sampling 
                if spp.t == self.next_t:

                    #get the data directory name for this timestep
                    dirname = os.path.join(os.getcwd(),
                          'GNX_mod-%s' % self.model_name,
                                    'it-%i' % iteration)

                    #get the subdirectory for this species
                    subdirname = os.path.join(dirname, 'spp-%s' % spp.name)

                    #and create (and its parent data directory, if needed)
                    os.makedirs(subdirname, exist_ok = True)

                    #get filenames
                    gen_files, geo_files = self._make_filenames(
                                    iteration = iteration, spp_name = spp.name)

                    #sample individuals according to the scheme defined 
                    sample = self._get_sample(spp)

                    #write files, if sample length > 0 
                    #(NOTE: otherwise, an empty file with "ZERO_SAMPLE" in the 
                    #name will be written, below)
                    if len(sample) > 0:

                        #save genetic data, if the spp has a
                        #genomic architeecture
                        if spp.gen_arch is not None:
                            #for each genetic data format to be written
                            for n, data_format in enumerate(self.gen_formats):

                                #format the data accordingly
                                data = self._format_gen_data(
                                    data_format = data_format, sample = sample,
                                                                    spp = spp)

                                #then write it to disk
                                gen_filepath = os.path.join(subdirname,
                                                            gen_files[n])
                                self._write_gendata(filepath = gen_filepath,
                                                            gen_data = data)

                        #also write the geodata for this spp
                        for n, data_format in enumerate(self.geo_formats):
                            #write the geodata to disk
                            geo_filepath = os.path.join(subdirname,
                                                            geo_files[n])
                            self._write_geodata(filepath = geo_filepath,
                                               data_format = data_format,
                                               sample = sample)

                    #if sample was empty, write a placeholder file with name
                    #"<base_filename>_ZERO_SAMPLE"
                    else:
                        filenames = [gen_files, geo_files]
                        filenames = [*chain.from_iterable(filenames)]
                        filename = os.path.splitext(filenames[0])[0]
                        filename = filename + '_ZERO_SAMPLE'
                        filepath = os.path.join(subdirname, filename)
                        self._write_gendata(filepath, '')

                # write the non-neutral loci CSV file, if required,
                # and if this is the final timestep
                if self.nonneut_loc_format is not None:
                    self._write_nonneut_loc_file(spp, subdirname, iteration)

            #write the raster, if necessary
            if self.rast_format is not None:
                #for each Layer
                for lyr in land.values():
                    #get the raster filename
                    filename = 'mod-%s_it-%i_t-%i_lyr-%s.%s' % (
                        self.model_name, iteration, self.next_t, lyr.name,
                                    self.file_extension_dict[self.rast_format])
                    filepath = os.path.join(dirname, filename)
                    #and write it to disk
                    lyr._write_raster(filepath, self.rast_format)

            #update self.next_t to the next timestep to be sampled
            self._set_next_t()


    def _get_random_sample(self, individuals):
        if len(individuals) > self.n:
            sample = r.choice(individuals, size = self.n, replace = False)
        else:
            sample = individuals
        return(sample)


    def _get_point_sample(self, spp):
        #TODO: check if this should be sped up any more
        sample = [i for i,v in spp.items() if self.pts_area.contains(
                                                        Point(v.x, v.y))]
        if len(sample) > self.n:
            sample = self._get_random_sample(individuals = sample)
        return(sample)


    def _get_sample(self, spp):
        #get a set of indices for the individuals in the sample
        sample = set()
        #take the whole species, if scheme is 'all'
        if self.scheme == 'all':
            sample.update([*spp])
        #or take a random sample, if scheme is 'random'
        elif self.scheme == 'random':
            inds = self._get_random_sample([*spp])
            sample.update(inds)
        #or take individuals within a given radius of self.buffs (which could
        #have come from a set of input points or from a calculated set of
        #transect points), if scheme is 'point' or 'transect'
        elif self.scheme in ['point', 'transect']:
            inds = self._get_point_sample(spp)
            sample.update(inds)
        #sort sample, then convert to a dict of individuals
        sample = sorted([*sample])
        sample = {i:spp[i] for i in sample}
        return(sample)


    def _format_gen_data(self, data_format, sample, spp):
        '''<data_format> can be:
                            'fasta'
                            'vcf'
        '''
        genotypes = spp._get_genotypes(individs=[*sample], as_dict=True)
        if data_format == 'fasta':
            formatted_data = _format_fasta(sample, genotypes)
        elif data_format == 'vcf':
            formatted_data = _format_vcf(sample, genotypes, spp.gen_arch,
                    include_fixed_sites = self.include_fixed_sites)
        return(formatted_data)


    def _write_gendata(self, filepath, gen_data):
        _write_file(filepath, gen_data)


    def _write_geodata(self, filepath, data_format, sample):
        write_fn = self.write_geodata_fn_dict[data_format]
        write_fn(filepath = filepath, individuals = sample)


    def _write_nonneut_loc_file(self, spp, subdir, iteration,
                                format='csv'):
        """
        NOTE: this is agnostic of the fact that the number of non-neutral
              loci can grow throughout a simulation because of mutation.
              the user would need to compare across nonneutral locus files,
              and/or compare with the mutation log,
              in order to determine when loci became non-neutral.
        """
        # get all the non-neutral loci for each trait
        locs_dict = {trt.name:
                     [*trt.loci] for trt in spp.gen_arch.traits.values()}
        # make the number of rows even across columns
        max_nrow = np.max([len(v) for v in locs_dict.values()])
        for k,v in locs_dict.items():
            if len(v) < max_nrow:
                locs_dict[k] = np.array([*v] + [np.nan]*(max_nrow-len(v)))
        # recast as a DataFrame
        locs_df = pd.DataFrame.from_dict(locs_dict)
        # make the filename and filepath
        filename = 'mod-%s_it-%i_t-%i_spp-%s_NONNEUTS.%s' % (self.model_name,
                                                             iteration,
                                                             self.next_t,
                                                             spp.name,
                                                             format)
        filepath = os.path.join(subdir, filename)
        # write the file, with trait names as column names
        # and locus numbers down the non-fixed-length columns
        locs_df.to_csv(filepath, index=False)



######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#a function to get a set of n evenly spaced points between endpoints
def _get_transect_points(endpoints, n):
    x_pts = np.linspace(endpoints[0][0] , endpoints[1][0], n)
    y_pts = np.linspace(endpoints[0][1] , endpoints[1][1], n)
    return(list(zip(x_pts, y_pts)))


#a function to make shapely geometry buffers around a set of points
def _make_point_buffers(points, radius):
    pts = [Point(*p) for p in points]
    buffs = [p.buffer(radius) for p in pts]
    cu = cascaded_union(buffs)
    return(cu)


# a function that returns an ad hoc sample of size n
# from the current population of the given spp,
# in the format index-keyed dict format needed by the data-writing functions
# (to be used by the convenience methods of the Model object)
def _get_adhoc_sample(spp, n):
    #get a set of indices for the individuals in the sample
    sample = set()
    #take the whole species, if scheme is 'all'
    if n is None:
        sample.update([*spp])
    #or take a random sample, if scheme is 'random'
    else:
        if len(spp) > n:
            inds = r.choice([*spp], size = n, replace = False)
        else:
            inds = [*spp]
        sample.update(inds)
    #sort sample, then convert to a dict of individuals
    sample = sorted([*sample])
    sample = {i:spp[i] for i in sample}
    return(sample)


def _format_fasta(sample, genotypes):
    '''
    FASTA FORMAT:

    >idx:haploid_num|x_location|y_location|phenotype0;phenotype1;...;
                    phenotypeN|env_var0;env_var1;...;env_varN|age|sex
    001110101010101010010101011101010110.....01011110

    '''
    # ensure that sample and genotypes have identical lenghts and orders
    assert len(sample) == len(genotypes), ("'sample', and 'genotypes' have"
                                           "different lenghts!")
    assert np.all([*sample] == [*genotypes]), ("'sample' and 'genotypes' do not"
                                               " have identical orders!")

    row1 = '>%s:HAP;%s;%s;%s;%s;%s;%s\n'
    file_text = ''

    for individ, genotype in zip(sample.values(), genotypes.values()):
        for hap in range(2):
            individ_row1 = re.sub('HAP', str(hap), row1)
            replace = tuple(map(lambda att: re.sub(',', '|', re.sub('[\[\] ]',
                '', str(getattr(individ, att)))), ['idx', 'x', 'y', 'age',
                                            'sex', 'z', 'e']))
            individ_row1 = individ_row1 % replace
            individ_row2 = ''.join([str(
                            base) for base in genotype[:,hap]]) + '\n'

            file_text = file_text + individ_row1 + individ_row2

    return(file_text)


def _format_vcf(sample, genotypes, gen_arch, include_fixed_sites=False):

    # ensure that sample and genotypes have identical lenghts and orders
    assert len(sample) == len(genotypes), ("'sample', and 'genotypes' have"
                                           "different lenghts!")
    assert np.all([*sample] == [*genotypes]), ("'sample' and 'genotypes' do not"
                                               " have identical orders!")

    #create a template header
        #NOTE: has 1 string slot for a date

        #TODO: DECIDE ON NECESSARY INFO AND FORMAT CONTENTS,
        #THEN ADD METADATA ROWS HERE

    header = '''##fileformat=VCFv4.2
##fileDate=%s
##source=Geonomics
'''

    #template column-header row
    col_header_row = ('#CHROM\tPOS\tID\tREF\tALT\tQUAL'
                      '\tFILTER\tINFO\tFORMAT\t%s\n')
        #NOTE: this has 1 string slot for a tab-separated
        #list of all individ ids

    #template data row
    #TODO: UPDATE/CHANGE THE INFO AND FORMAT PORTIONS OF THIS TEMPLATE,
    #AFTER I DECIDE ON THEIR CONTENTS (above)
    data_row = ('%i\t%i\t.\tA\tT\t1000\tPASS\t%s\tGT\t%s\n')
        #NOTE: this has 2 integer slots, then 1 string slot for:
            #- chrom number (NOTE: unpythonically, starts from 1)
            #- locus number (NOTE: reported cumulative from locus 0,
               #not from start of each chrom)
            #- a tab-separated list of individs' genotypes at this locus

    #create a col_header_row for this data
    inds = [*sample.keys()]
    ind_cols = '\t'.join([str(i) for i in inds])
    cols = col_header_row % (ind_cols)

    #get a list of the chromosome numbers
    #chroms = np.cumsum(gen_arch.l_c)

    #and get all individuals' genomic data in a 'samplome' object (a 3-d array)
    samplome = np.array([genotypes[i] for i in inds])

    #get all segregating sites
    max_val = 2 * len(sample)
    segs = np.where(samplome.sum(axis = 2).sum(axis = 0) > 0)[0]
    segs2 = np.where(samplome.sum(axis = 2).sum(axis = 0) < max_val)[0]
    segs = sorted(list(set(segs).intersection(set(segs2))))

    #if not_include_fixed_sites, include only the segregating sites in the VCF 
    if not include_fixed_sites:
        loci = segs
    #or else get all loci
    else:
        loci = range(gen_arch.L)

    #and get the sites' chrom nums
    #chroms = [list((locus - chroms) < 0).index(True) for locus in loci]

    #build all the VCF data rows
    rows = ''
    for n, locus in enumerate(loci):
        gts = samplome[:,locus,:]
        gts = '\t'.join(['%i|%i' % (gts[i,0],
                gts[i,1]) for i in range(np.shape(gts)[0])])
        #get an indicator for whether the site is segregating or fixed
        seg = locus in segs
        seg_dict = {True: 'SEG', False: 'FIX'}

        rows = rows + data_row % (0, locus, seg_dict[seg], gts)

    #get the date
    now = datetime.datetime.now()
    month = str(now.month).zfill(2)
    day = str(now.day).zfill(2)
    date = '%d%s%s' % (now.year, month, day)

    #paste all the VCF content together
    out_vcf = ''.join([header % date, cols, rows])

    #return it
    return(out_vcf)


def _tskit_table_to_pandas(table):
    """
    Takes a raw tskit.TableCollection table dict,
    as returned from tskit.TableCollection.<TABLENAME>.asdict(),
    and converts it into a pandas.DataFrame, for write-out
    """
    # get the max column length
    max_col_len = max([len(col) for col in table.values()])
    # fill in columns with np.nan
    filled_table = {}
    for k, col in table.items():
        if len(col) > 0:
            filled_table[k] = np.concatenate((col,
                                    np.nan * np.ones(max_col_len - len(col)))) 
        else:
            filled_table[k] = np.ones(max_col_len) * np.nan
    # convert to pandas.DataFrame
    pd_table = pd.DataFrame.from_dict(filled_table)
    return pd_table

