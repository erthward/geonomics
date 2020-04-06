#!/usr/bin/python
# io.py

'''
Defines core IO functions
'''

# genomics imports

# other imports
import os
import csv
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import rasterio


######################################
# -----------------------------------#
# VARIABLES -------------------------#
# -----------------------------------#
######################################


# TODO: For now, defining a default projection for saving geotiff rasters
# when no projection is set, but I need to decide if it even makes sense to do
# this, rather than just disallowing geotiff as a raster-saving format if no
# projection is set (because of course this makes random layers project as
# rasters at lat:0, lon:0 in Web Mercator, which is at the Equator and on the
# prime meridian (basically due west of São Tomé and due south of Accra)
_DEFAULT_PROJ = '''PROJCS["WGS 84 / Pseudo-Mercator",
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


#    ########
#    # Read #
#    ########

def _read_raster(filepath, coord_prec, dim=None):
    if os.path.splitext(filepath)[1].lower() == '.txt':
        rast = np.fromfile(filepath, sep=' ')
        assert len(rast) == np.prod(dim), ('The raster read in from the .txt '
                                           'file provided does not have a '
                                           'size equal to the product of the '
                                           'dimensions provided, and thus '
                                           'cannot be coerced to an array of '
                                           'those dimensions.')
        rast = np.float32(rast.reshape(dim))
        dim = rast.shape[::-1]
        res = (1, 1)
        ulc = (0, 0)
        prj = None
    else:
        rast_file = rasterio.open(filepath)
        rast = rast_file.read()[0, :, :]
        dim = rast.shape[::-1]
        # NOTE: rasterio has switched to using the proper affine transform
        # matrix as the dataset.transform attribute (see, for example,
        # https://www.perrygeo.com/python-affine-transforms.html). However, the
        # get_transform() method still returns the affine transform's tuple
        # in GDAL format, which is specified (according to the
        # tutorial at https://www.gdal.org/gdal_tutorial.html) as:
        # adfGeoTransform[0] /* top left x */
        # adfGeoTransform[1] /* w-e pixel resolution */
        # adfGeoTransform[2] /* 0 */
        # adfGeoTransform[3] /* top left y */
        # adfGeoTransform[4] /* 0 */
        # adfGeoTransform[5] /* n-s pixel resolution (negative value) */
        res = tuple([i for n, i in enumerate(rast_file.get_transform(
                                                        )) if n in [1, 5]])
        res = np.round(res, coord_prec)
        ulc = tuple([i for n, i in enumerate(rast_file.get_transform(
                                                        )) if n in [0, 3]])
        ulc = np.round(ulc, coord_prec)
        # get the projection as a CRS object
        prj = rast_file.crs
    return(rast, dim, res, ulc, prj)


# read a txt file containing a stack of 2d arrays (such as is created to save 
# LD data)
def _read_array_stack(filepath):
    array = np.loadtxt(filepath)
    return(array)

    #########
    # Write #
    #########


# write data (a block of text, as a string) to a file
def _write_file(filepath, data):
    with open(filepath, 'w') as f:
        f.write(data)


# write a CSV of data from a dictionary of columns
def _write_dict_to_csv(filepath, array_1d_dict):
    # make a pandas DataFrame from the species' stats
    df = pd.DataFrame.from_dict(array_1d_dict)
    # add a timestep column
    df['t'] = [*df.index]
    # reorder the columns so that timestep is the first
    ordered_cols = ['t'] + [col for col in df.columns if col != 't']
    df = df[ordered_cols]
    # write to disk
    df.to_csv(path_or_buf=filepath, float_format='%0.5f', index=False)


# append a 2d array layer to an array stack (i.e. an np.txt file
# intended to eventually be read in as a 3D array)
def _append_array2d_to_array_stack(filepath, locuswise_array2d):
    tmp_file_dir = os.path.split(filepath)[0]
    tmp_filename = 'tmp_%s.txt' % str(np.random.randint(0, 1000)).zfill(4)
    tmp_filename = os.path.join(tmp_file_dir, tmp_filename)
    np.savetxt(tmp_filename, locuswise_array2d, fmt='%0.5f')
    with open(tmp_filename, 'r') as f:
        tmp_txt = f.read()
    with open(filepath, 'a') as f:
        f.write(tmp_txt)
    os.remove(tmp_filename)


# append a row of data to a CSV file
def _append_row_to_csv(filepath, locuswise_array1d, t):
    data_dict = dict(zip(range(locuswise_array1d.size), locuswise_array1d))
    write_header = not os.path.exists(filepath)
    with open(filepath, 'a') as f:
        dict_writer = csv.DictWriter(f, ['t'] + [*data_dict.keys()])
        if write_header:
            dict_writer.writeheader()
        dict_writer.writerow({'t': t, **data_dict})


# write a data file from a geopandas object created from an index-keyed dict
# of individual.Individual objects
def _write_geopandas(filepath, individuals, driver):
    # get full path and filename
    attributes = ['idx', 'z', 'e', 'age', 'sex']
    # TODO: FIXME: replace the call to str() below with something more
    # sophisticated that will actually separate distinct phenotype
    # and environment values into separate, labeled columns
    # TODO: also round numbers to 5 decimals or so?
    df_dict = {att: [str(getattr(
            ind, att)) for ind in individuals.values()] for att in attributes}
    pts = [Point(ind.x, ind.y) for ind in individuals.values()]
    df_dict['pt'] = pts
    df = pd.DataFrame.from_dict(df_dict)
    gdf = gpd.GeoDataFrame(df, geometry='pt')
    if driver == 'CSV':
        gdf['x'] = gdf.pt.x
        gdf['y'] = gdf.pt.y
        gdf = gdf.drop(labels=('pt'), axis=1)
        gdf.to_csv(filepath, index=False)
    else:
        gdf.to_file(filepath, driver=driver)


# write a shapefile from an index-keyed dict of individual.Individual objects
def _write_shapefile(filepath, individuals):
    filepath = _set_extension(filepath, 'shp')
    _write_geopandas(filepath, individuals, driver='ESRI Shapefile')


# write a geojson from an index-keyed dict of individual.Individual objects
def _write_geojson(filepath, individuals):
    filepath = _set_extension(filepath, ['json', 'geojson'])
    _write_geopandas(filepath, individuals, driver='GeoJSON')


# write a csv from an index-keyed dict of individual.Individual objects
def _write_csv(filepath, individuals):
    filepath = _set_extension(filepath, 'csv')
    _write_geopandas(filepath, individuals, driver='CSV')


# write a txt array from a landscape.Layer object's numpy-array raster
def _write_txt_array(filepath, lyr):
    filepath = _set_extension(filepath, 'txt')
    np.savetxt(filepath, lyr.rast, fmt='%0.5f')


# write a geotiff from a landscape.Layer object's numpy-array raster
def _write_geotiff(filepath, lyr):
    filepath = _set_extension(filepath, ['tif', 'tiff'])
    # create a dict that contains all the necessary outfile metadata
    out_meta = {'driver': 'GTiff',
                'dtype': 'float32',     # always save as continuous data
                'nodata': None,
                'width': lyr.dim[0],
                'height': lyr.dim[1],   
                'count': 1,             # n layers
                'tiled': False,
                'compress': 'lzw',
                'interleave': 'band'
               }

    # take the Layer's projection, if it has one, or else just assign
    # it unprojected lat-lon, WGS84
    if lyr.prj is None:
        out_meta['crs'] = rasterio.crs.CRS.from_dict(init='epsg:4326')
    else:
        out_meta['crs'] = lyr.prj

    # use the res and ulc info to create the transform object
    aff = rasterio.transform.Affine(res[0],
                                    0,
                                    ulc[0],
                                    0,
                                    res[1],
                                    ulc[1]
                                   )
    out_meta['transform'] = aff

    with rasterio.open(filepath, 'w', **out_meta) as dst:
        dst.write(lyr.rast, 1)
    return


#    #########
#    # Other #
#    #########

# add the correct extension to a filepath, if necessary
def _set_extension(filepath, ft_ext):
    # make ft_ext iterable, if just a string is provided
    if type(ft_ext) is str:
        ft_ext = [ft_ext]

    # check if filepath has any extension
    any_ext = len([i for i in os.path.splitext(filepath) if len(i) > 0]) > 1
    # and check for an existing valid extension
    valid_ext = os.path.splitext(filepath.lower())[1].lstrip('.') in ft_ext

    # if the file has an extension not concordant with the filetype extension
    # (ft_ext), throw an error
    if any_ext and not valid_ext:
        raise ValueError(("File name already contains an extension "
                          "('%s'), but it is incompatible with the filetype "
                          "to be written (which requires one of the following "
                          " extensions: %s).") % (os.path.splitext(
                            filepath)[1], ','.join(
                            ["'.%s'" % ext for ext in ft_ext])))

    # else if it has no extension, append the extension
    elif not any_ext and not valid_ext:
        filepath = '.'.join([filepath, ft_ext[0]])
    # else, leave it be
    else:
        pass
    # return the resulting filepath
    return filepath
