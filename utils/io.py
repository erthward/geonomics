#!/usr/bin/python
#io.py

def read_raster(filepath):
    from osgeo import gdal
    rast = gdal.Open(filepath)
    array = rast.ReadAsArray()
    dim = array.shape
    res = [x for n,x in enumerate(rast.GetGeoTransform()) if n in [1,5]]
    ulc = [x for n,x in enumerate(rast.GetGeoTransform()) if n in [0,3]]
    return(array, dim, res, ulc)


