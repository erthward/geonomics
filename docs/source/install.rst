.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash

Installation
************

Some Geonomics dependencies (i.e. other Python packages that Geonomics uses)
depend in turn on a stack of open-source geospatial libraries
in C (`GEOS`_, `GDAL`_, and `PROJ`_). They can sometimes be tricky to install.
Thus, if you are unsure if you already have a working install of those
libraries, then we recommend :ref:`Installing with \`conda\```.

Alternatively, you can try :ref:`Installing with \`pip\```. That will work fine,
as long as `pip` can successfully install Geonomics' dependencies as well.

You can find a list of Geonomics' :ref:`Dependencies` below.


Installing with `conda`
-----------------------

You may want to install with `conda`_, as it will automatically install
third-party dependencies (such as `GEOS`_, `GDAL`_, and `PROJ`_) for you
as well. 

The `conda` software can be installed by following
the installation instructions on the `conda`_ webpage.

Then, Geonomics can be installed using:

.. code-block:: bash
  
    conda config --env --add channels conda-forge
    conda install geonomics

If you prefer, you may first create a new `conda` environment, then install
Geonomics there. This will avoid potential dependency conflicts with other
packages using Geonomics' dependencies. You can do that by running:

.. code-block:: bash

  conda create -n gnx
  conda activate gnx
  conda config --env --add channels conda-forge
  conda install python=3 geonomics

Once installed, Geonomics can then be tested by launching a Python prompt and running:

.. code-block:: python
    
    import geonomics as gnx
    gnx.run_default_model()

This will load the Geonomics package, create in your current working
directory a new Geonomics parameters file containing the default
parameter values, use that file to instantiate and run a Geonomics model,
then delete the parameters file (by default).


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

- `pyvcf <http://pyvcf.readthedocs.io/en/latest/>`_

Optional dependencies:
......................

- `nlmpy <http://pypi.org/project/nlmpy/>`_

- `GDAL (Python package) <http://pypi.org/project/GDAL/>`_




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
  `conda`_ to set up a clean `conda` environment, the install Geonomics
  and its dependencies there. While we are not sure why, this appears to
  make things all better.

  In :ref:`Installing with \`conda\``` we have provided code for
  creating a new `conda` environment and installing Geonomics there.



.. _GDAL: https://www.gdal.org/

.. _GEOS: https://geos.osgeo.org

.. _PROJ: https://proj.org/

.. _conda: https://docs.conda.io/en/latest/

.. _pip: https://pip.pypa.io/en/stable/
