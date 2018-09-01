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

#genomics imports

#other imports
import numpy as np
import pandas as pd
import geopandas as gpd
from osgeo import gdal
import os
from shapely.geometry import Point



######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################


    ########
    # Read #
    ########

def read_raster(filepath):
    if os.path.splitext(filepath)[1].lower() == 'txt':
        rast = np.fromfile(filepath, sep = ' ')
        assert len(rast) == prod(dim), ('The raster read in from the .txt '
        'file provided does not have a size equal to the product of the '
        'dimensions provided, and thus cannot be coerced to an array of '
        'those dimensions.')
        rast = rast.reshape(dim)
        dim = dim
        res = (1,1)
        ulc = (0,0)
        prj = None
    else:
        rast_file = gdal.Open(filepath)
        rast = rast_file.ReadAsArray()
        dim = rast.shape
        res = [i for n,i in enumerate(rast_file.GetGeoTransform()) if n in [1,5]]
        ulc = [i for n,i in enumerate(rast_file.GetGeoTransform()) if n in [0,3]]
        #get the projection as WKT
        prj = rast_file.GetProjection()
    return(rast, dim, res, ulc, prj)

    
    #########
    # Write #
    #########

def write_file(dirname, filename, data):
    filepath = os.path.join(dirname, filename)
    with open(filepath, 'w') as f:
        f.write(data)

def write_geopandas(dirname, filename, individuals, driver):
    #get full path and filename
    filepath = os.path.join(dirname, filename)
    attributes = ['idx', 'phenotype', 'habitat', 'age', 'sex']
    #FIXME: replace the call to str() below with something more sophisticated
    #that will actually separate distinct phenotype and habitat values into
    #separate, labeled columns
    df_dict = {att: [str(getattr(ind, att)) for ind in individuals] for att in attributes}
    pts = [Point(ind.x, ind.y) for ind in individuals]
    df_dict['pt'] = pts
    df = pd.DataFrame.from_dict(df_dict)
    gdf = gpd.GeoDataFrame(df, geometry = 'pt')
    if driver == 'CSV':
        gdf['x'] = gdf.pt.x
        gdf['y'] = gdf.pt.y
        gdf = gdf.drop(labels = ('pt'), axis = 1)
        gdf.to_csv(filepath, index = False)
    else:
        gdf.to_file(filepath, driver = driver)


def write_shapefile(dirname, filename, individuals):
    write_geopandas(dirname, filename, individuals, driver='ESRI Shapefile')


def write_geojson(dirname, filename, individuals):
    write_geopandas(dirname, filename, individuals, driver='GeoJSON')


def write_csv(dirname, filename, individuals):
    write_geopandas(dirname, filename, individuals, driver = 'CSV')


def write_array(dirname, filename, scape):
    filepath = os.path.join(dirname, filename)
    np.savetxt(filepath, scape.rast, fmt = '%0.5f')


def write_geotiff(dirname, filename, scape):
    #TODO: this is a tweak on code taken from https://gis.stackexchange.com/
    #questions/58517/python-gdal-save-array-as-raster-with-projection-from-
    #other-file
    #get the driver
    driver = gdal.GetDriverByName('GTiff')
    #get values
    #number of pixels in x and y
    x_pixels = scape.dim[0]
    y_pixels = scape.dim[1]
    #resolution
    PIXEL_SIZE = scape.res[0]
    #x_min & y_max are the "top-left corner"
    x_min = scape.ulc[0]
    y_max = scape.ulc[1] + (scape.dim[1]*scape.res[1]) 
    #get the WKT projection
    wkt_projection = scape.prj
    if wkt_projection is None:
        #TODO: SET SOME DEFAULT, MEANINGLESS PROJECTION??
        pass
    dataset = driver.Create(
       filepath,
       x_pixels,
       y_pixels,
       1,
       gdal.GDT_Float32, )
    dataset.SetGeoTransform((
       x_min,        # 0
       PIXEL_SIZE,   # 1
       0,            # 2
       y_max,        # 3
       0,            # 4
       -PIXEL_SIZE)) # 5
    #set the projection
    dataset.SetProjection(wkt_projection)
    #and write to disk
    dataset.GetRasterBand(1).WriteArray(scape.rast)
    dataset.FlushCache()


