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
def make_params_file(args):
    pass

def read_params_file(params_file):
    name = os.path.splitext(os.path.split(params_file)[-1])[0]
    params = sim.params.read(params_file)
    params['name'] = name
    return(params)

def make_model(params_file, verbose=False):
    mod = model.Model(name, params, verbose=verbose)
    return(mod)


