#!/usr/bin/python
# _str_repr_.py

'''
##########################################

Module name:          utils.str_repr


Module contains:
                      - helper functions for defining class __str__ and
                      __repr__ special methods


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           10-13-18
Documentation:        URL


##########################################
'''


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#function to get the spaces to insert between an item and its value for a line
#of output in a class' __repr__ string
def _get_str_spacing(unformatted_str):
    tot_spacing = 48
    spacing_len = tot_spacing - (1 + len(unformatted_str.split(':')[0]))
    spacing = ' ' * spacing_len
    return spacing

