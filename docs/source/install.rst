.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash

Installation
************

Some Geonomics dependencies (i.e. other Python packages that Geonomics uses)
depend in turn on a stack of open-source geospatial libraries
in C (`GEOS`_, `GDAL`_, and `PROJ`_). 

They can sometimes be tricky to install. Please make sure you have them installed first.
Then you can follow the `pip` installation instructions below.

Installing with `pip`
---------------------

If you already have `GEOS`_, `GDAL`_, and `PROJ`_ installed, then installing
with `pip`_ probably will not be an issue. 

The `pip` software can be installed by following
the installation instructions on the `pip`_ webpage.

Then, Geonomics can be installed using:

.. code-block:: bash

     pip install geonomics

The installation can then be tested by launching a Python prompt and running:

.. code-block:: python
    
    import geonomics as gnx
    gnx.run_default_model()

This will load the Geonomics package, create in your current working
directory a new Geonomics parameters file containing the default
parameter values, use that file to instantiate and run a Geonomics model,
then delete the parameters file (by default).


Dependencies
------------

Required dependencies:
......................

- `numpy <http://numpy.org/>`_

- `matplotlib <http://matplotlib.org/>`_

- `pandas <http://pandas.pydata.org/>`_ (version 0.23.4 or later)

- `geopandas <http://geopandas.org/>`_

- `scipy <http://www.scipy.org/scipylib/index.html>`_ (version 1.3.1 or later)

- `scikit-learn <http://scikit-learn.org/stable/>`_

- `statsmodels <http://www.statsmodels.org/stable/index.html>`_ (version
  0.9.0 or later)

- `shapely <http://shapely.readthedocs.io/en/stable/project.html>`_

- `bitarray <http://pypi.org/project/bitarray/>`_

- `rasterio <https://rasterio.readthedocs.io/en/latest/index.html>`_

Optional dependencies:
......................

- `nlmpy <http://pypi.org/project/nlmpy/>`_



------------------------------------------------------------

Troubleshooting
---------------

Here is a list of issues that we have seen come up during Geonomics
installation, each with some suggestions for troubleshooting:


Mac Segmentation Fault
......................

**Geonomics installs, but running the default model crashes Python with a 
segmentation fault**:
    
  We have only seen this issue on older Macs (ca. 2012-2015). The segmentation
  fault appears to be caused by functions within
  `scipy <http://www.scipy.org/scipylib/index.html>`_'s interpolate module
  (e.g. `scipy.interpolate.griddata`), such as discussed
  `here <https://stackoverflow.com/questions/59274750/segmentation-fault-when-running-scipy-interpolate>`_). The only fix we are aware of, as of now, is to use
  `conda`_ to set up a clean `conda` environment, then install Geonomics
  and its dependencies there. While we are not sure why, this appears to
  make things all better.


.. _GDAL: https://www.gdal.org/

.. _GEOS: https://geos.osgeo.org

.. _PROJ: https://proj.org/

.. _conda: https://docs.conda.io/en/latest/

.. _pip: https://pip.pypa.io/en/stable/
