.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


############
Installation
############

Geonomics depends on a number of common Python packages. See the
[Dependencies] section below for details.

*********************
Installing with `pip`
*********************

To install Geonomics we recommend using
`pip <https://pypi.org/project/pip/>`. It can be installed by following
the installation instructions on the `pip <https://pypi.org/project/pip/>`
webpage. Geonomics can then be installed using:

.. code-block:: bash

     pip install geonomics

The installation can then be tested by launching a Python prompt and running:

.. code-block:: python
    
    import geonomics as gnx
    gnx.run_default_model()

This will load the Geonomics package, create in your current working
directory a new Geonomics parameters file containing the default
parameter values, then use that file to instantiate and run a Geonomics model.

************
Dependencies
************

Required dependencies:

`numpy <http://numpy.org/>`_

`matplotlib <http://matplotlib.org/>`_

`pandas <http://pandas.pydata.org/>`_ (version 0.23.4 or later)

`geopandas <http://geopandas.org/>`_

`scipy <http://www.scipy.org/scipylib/index.html>`_

`scikit-learn <http://scikit-learn.org/stable/>`_

`statsmodels <http://www.statsmodels.org/stable/index.html>`_ (version
0.9.0 or later)

`shapely <http://shapely.readthedocs.io/en/stable/project.html>`_

`bitarray <http://pypi.org/project/bitarray/>`_

`pyvcf <http://pyvcf.readthedocs.io/en/latest/>`_

Optional dependencies:

`nlmpy <http://pypi.org/project/nlmpy/>`_

`GDAL <http://pypi.org/project/GDAL/>`_
