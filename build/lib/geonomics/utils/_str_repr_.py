#!/usr/bin/python
# _str_repr_.py

'''
Defines core functions to aid in formatting strings for STDOUT
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

