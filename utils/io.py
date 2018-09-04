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
# VARIABLES -------------------------#
# -----------------------------------#
######################################


#TODO: For now, defining a default projection to use for saving geotiff rasters
#if no projection is set, but I need to decide if it even makes sense to do
#this, rather than just disallowing geotiff as a raster-saving format if no
#projection is set (because of course this makes random landscapes project as
#rasters at lat:0, lon:0 in Web Mercator, which is at the Equator and on the
#prime meridian (basically due west of São Tomé and due south of Accra)
__DEFAULT_WEB_MERCATOR__ = '''PROJCS["WGS 84 / Pseudo-Mercator",
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.0174532925199433,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]],
    PROJECTION["Mercator_1SP"],
    PARAMETER["central_meridian",0],
    PARAMETER["scale_factor",1],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["X",EAST],
    AXIS["Y",NORTH],
    EXTENSION["PROJ4","+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"],
    AUTHORITY["EPSG","3857"]]'''

######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################


    ########
    # Read #
    ########

def read_raster(filepath, dim=None):
    if os.path.splitext(filepath)[1].lower() == '.txt':
        rast = np.fromfile(filepath, sep = ' ')
        assert len(rast) == np.prod(dim), ('The raster read in from the .txt '
        'file provided does not have a size equal to the product of the '
        'dimensions provided, and thus cannot be coerced to an array of '
        'those dimensions.')
        rast = np.float32(rast.reshape(dim))
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

#write a data (a block of text, as a string) to a file
def write_file(filepath, data):
    with open(filepath, 'w') as f:
        f.write(data)

#write a data file from a geopandas object created from an index-keyed dict 
#of individual.Individual objects
def write_geopandas(filepath, individuals, driver):
    #get full path and filename
    attributes = ['idx', 'phenotype', 'habitat', 'age', 'sex']
    #TODO: FIXME: replace the call to str() below with something more sophisticated
    #that will actually separate distinct phenotype and habitat values into
    #separate, labeled columns
    #TODO: also round numbers to 5 decimals or so?
    df_dict = {att: [str(getattr(ind, att)) for ind in individuals.values()] for att in attributes}
    pts = [Point(ind.x, ind.y) for ind in individuals.values()]
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


#write a shapefile from an index-keyed dict of individual.Individual objects
def write_shapefile(filepath, individuals):
    filepath = set_extension(filepath, 'shp')
    write_geopandas(filepath, individuals, driver='ESRI Shapefile')


#write a geojson from an index-keyed dict of individual.Individual objects
def write_geojson(filepath, individuals):
    filepath = set_extension(filepath, ['json', 'geojson'])
    write_geopandas(filepath, individuals, driver='GeoJSON')


#write a csv from an index-keyed dict of individual.Individual objects
def write_csv(filepath, individuals):
    filepath = set_extension(filepath, 'csv')
    write_geopandas(filepath, individuals, driver = 'CSV')


#write a txt array from a landscape.Scape object's numpy-array raster
def write_txt_array(filepath, scape):
    filepath = set_extension(filepath, 'txt')
    np.savetxt(filepath, scape.rast, fmt = '%0.5f')


#write a geotiff from a landscape.Scape object's numpy-array raster
def write_geotiff(filepath, scape):
    filepath = set_extension(filepath, ['tif', 'tiff'])
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
        #TODO: FOR NOW THIS DEFAULTS TO THE __DEFAULT_WEB_MERCATOR__
        #VARIABLE, BUT THINK ABOUT WHETHER IT MAKES ANY SENSE TO HAVE
        #SUCH A DEFAULT, AND IF SO, WHETHER OR NOT THERE'S A BETTER
        #OPTION
        wkt_projection = __DEFAULT_WEB_MERCATOR__
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


    #########
    # Other #
    #########

#add the correct extension to a filepath, if necessary
def set_extension(filepath, ft_ext):
    #make ft_ext iterable, if just a string is provided
    if type(ft_ext) is str:
        ft_ext = [ft_ext]

    #check if filepath has any extension
    any_ext =  len([i for i in os.path.splitext(filepath) if len(i) >0]) > 1
    #and check for an existing valid extension
    valid_ext = os.path.splitext(filepath.lower())[1].lstrip('.') in ft_ext

    #if the file has an extension not concordant with the filetype extension
    #(ft_ext), throw an error
    if any_ext and not valid_ext:
        raise ValueError(("File name already contains an extension "
            "('%s'), but it is incompatible with the filetype "
            "to be written (which requires one of the following "
            " extensions: %s).") % (os.path.splitext(filepath)[1],
            ','.join(["'.%s'" % ext for ext in ft_ext])))

    #else if it has no extension, append the extension
    elif not any_ext and not valid_ext:
        filepath = '.'.join([filepath, ft_ext[0]])
    #else, leave it be
    else:
        pass
    #return the resulting filepath
    return filepath
