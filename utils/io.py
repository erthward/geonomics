#!/usr/bin/python
#io.py

'''
##########################################

Module name:          utils.io


Module contains:
                      - basic IO functions


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           07-18-18
Documentation:        URL


##########################################
'''

######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def read_raster(filepath):
    from osgeo import gdal
    rast = gdal.Open(filepath)
    array = rast.ReadAsArray()
    dim = array.shape
    res = [x for n,x in enumerate(rast.GetGeoTransform()) if n in [1,5]]
    ulc = [x for n,x in enumerate(rast.GetGeoTransform()) if n in [0,3]]
    return(array, dim, res, ulc)

