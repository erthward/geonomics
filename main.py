#!/usr/bin/python
# main.py

'''
##########################################

Module name:          main


Module contains:
                      - the Geonomics main module, containing the key, 
                        highest-level functions the common user would need


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           07-06-18
Documentation:        URL


##########################################
'''

#geonomics imports
from sim import model

#other imports
import re
import os, sys
import numpy as np
import pandas as pd

######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#function that takes a bunch of arguments and returns a template params file to
#be filled in; args should include:
    #num of scapes and whether they should change
    #num of populations, whether they should have genomes, move, and change
    #whether data and stats should be generated
#and they should have sane defaults
#TODO
def make_params_file(filename, n_scapes, n_pops, scape_change=None,
                     pop_change=None, movement=True, genomes=True, data=True,
                     custom_gen_arch=None):
    #define a function to generate an assertion error string for arguments that
    #can be None, a boolean, or a list of booleans of a certain length
    def make_boollist_assert_str(arg, object_ref):
        arg_name = [ k for k,v in locals().items() if id(v) == id(arg) and not
                    k.startswith('_') and re.search('^[a-z]', k)]
        string = ("The %s argument can be None (default), a boolean (which "
                  "will apply to all %ss), or a list of booleans of length "
                  "equal to the number of %ss.") 
        string = string % (arg_name, object_reference, object_reference)
        return(string)
    #define a function to run asserts on arguments that can be None, a boolean,
    #or a list of booleans of a certain length
    def assert_boollist(arg, target_object, target_len):
        assert type(arg) in (type(None), bool, list), make_boollist_assert_str(arg, target_object)
        if type(arg) is list:
            assert (np.all([type(x) is bool for x in arg]) 
            and len(arg) == target_len), make_boollist_assert_str(arg, target_len)
    #run assertions to check argument type validity
    assert type(filename) is str, "The filename must be provided as a string."
    assert type(n_scapes) is int, "The number of scapes must be provided as an integer."
    assert type(n_pops) is int, "The number of populations must be provided as an integer."
    assert_boollist(scape_change, 'scape', n_scapes)
    for arg in [pop_change, genomes, movement, custom_gen_arch]:
        assert_boollist(arg, 'population', n_pops)
    #coerce requisite arguments from bool to a list of bools
    if type(scape_change) is bool:
        scape_change = [scape_change] * n_scapes
    if type(pop_change) is bool:
        pop_change = [pop_change] * n_pops
    if type(movement) is bool:
        movement = [movement] * n_pops
    if type(genomes) is bool:
        genomes = [genomes] * n_pops
    if type(custom_gen_arch) is bool:
        custom_gen_arch = [custom_gen_arch] * n_pops
        #should only be made True if that population will have genomes
        for n, pop in enumerate(genomes):
            if not pop:
                custom_gen_arch[n] = False

    #generate a dictionary of params-file sections for each scape and pop
    scapes_dict = {}
    pops_dict = {}

    #generate a params file template, using the arguments above, and to save it
    #to the filename provided

    #and if required, generate a custom genomic-architecture file for each
    #applicable population, saving it as '<params_filename>_pop_<n>_gen_arch.csv'
    if custom_gen_arch is not None:
        if True in custom_gen_arch:
            for n, val in enumerate(custom_gen_arch):
                if val:
                    #create the dataframe for the CSV file
                    cols = ('locus', 'p', 'r', 'trait', 'alpha')
                    row1 = ([0], [np.nan], [0.5], [np.nan], [np.nan])
                    df_dict = dict(zip(cols, row1))
                    df = pd.DataFrame.from_dict(df_dict)
                    #write it to file, without the index
                    df.to_csv(os.path.splitext(filename)[0] + 
                              '_pop%i_gen_arch_.csv' % n, index = False)



def read_params_file(params_file):
    name = os.path.splitext(os.path.split(params_file)[-1])[0]
    params = sim.params.read(params_file)
    params['name'] = name
    return(params)


def make_model(params_file, verbose=False):
    mod = model.Model(name, params, verbose=verbose)
    return(mod)


