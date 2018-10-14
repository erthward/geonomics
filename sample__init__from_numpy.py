#THIS IS THE __init__.py SCRIPT TAKEN FROM NUMPY AND TWEAKED A BIT. CONTINUE TO
#STUDY THIS TO FIGURE OUT HOW TO SET UP MY OWN!

"""
#########
Geonomics
#########

An individual-based, spatially explicit, forward-time, and highly
customizable landscape-genomics simulator written in Python.

The package provides:

  #. Land, GenomicArchitecture, Individual, Population, and Community
     classes
  #. A Model class that builds these objects according to the scenario
     stipulated in a parameters file, uses them to run numerous simulation
     iterations, and collects data and basic statistics from each iteration
  #. Functionality for planning arbitrarily complex landscape- and
     demographic-change scenarios
  #. Numerous visualization methods, to help users design, run, explore,
     present, and explain their models' behavior and results


*****************
The documentation
*****************
Documentation is available in two forms: as docstrings provided
with the code, and as a reference guide on the `geonomics webpage
<http://www.PAGE_HERE.com>`_.

Examples in the documentation assume `genomics` has been imported as `gnx`::

  >>> import geonomics as gnx

Code snippets in the examples are indicated by three right carets::

  >>> gnx.make_parameters_file()
  >>> gnx.make_model() 
  >>> model.run(verbose = True)

To view a function's docstring, use the Python `help` function::

  >>> help(gnx.make_model)

Or use the `?` shorthand, if you're working in iPython::

  >>> ?gnx.make_model


#####################################
#####################################
#####################################
#####################################
TODO: CHANGE THESE TO MY SUBPACKAGES!
#####################################
#####################################
#####################################
#####################################

***********
Subpackages
***********
doc
    Topical documentation on broadcasting, indexing, etc.
lib
    Basic functions used by several sub-packages.
random
    Core Random Tools
linalg
    Core Linear Algebra Tools
fft
    Core FFT routines
polynomial
    Polynomial tools
testing
    NumPy testing tools
f2py
    Fortran to Python Interface Generator.
distutils
    Enhancements to distutils with support for
    Fortran compilers support and more.

---------
Utilities
---------
test
    Run numpy unittests
show_config
    Show numpy build configuration
dual
    Overwrite certain functions with high-performance Scipy tools
matlib
    Make everything matrices.
__version__
    NumPy version string

"""

import sys
import warnings

from ._globals import ModuleDeprecationWarning, VisibleDeprecationWarning
from ._globals import _NoValue

# We first need to detect if we're being called as part of the numpy setup
# procedure itself in a reliable manner.
try:
    __NUMPY_SETUP__
except NameError:
    __NUMPY_SETUP__ = False

if __NUMPY_SETUP__:
    sys.stderr.write('Running from numpy source directory.\n')
else:
    try:
        from numpy.__config__ import show as show_config
    except ImportError:
        msg = """Error importing numpy: you should not try to import numpy from
        its source directory; please exit the numpy source tree, and relaunch
        your python interpreter from there."""
        raise ImportError(msg)

    from .version import git_revision as __git_revision__
    from .version import version as __version__

    from ._import_tools import PackageLoader

    def pkgload(*packages, **options):
        loader = PackageLoader(infunc=True)
        return loader(*packages, **options)

    from . import add_newdocs
    __all__ = ['add_newdocs',
               'ModuleDeprecationWarning',
               'VisibleDeprecationWarning']

    pkgload.__doc__ = PackageLoader.__call__.__doc__

    # Allow distributors to run custom init code
    from . import _distributor_init

    from . import core
    from .core import *
    from . import compat
    from . import lib
    from .lib import *
    from . import linalg
    from . import fft
    from . import polynomial
    from . import random
    from . import ctypeslib
    from . import ma
    from . import matrixlib as _mat
    from .matrixlib import *
    from .compat import long

    # Make these accessible from numpy name-space
    # but not imported in from numpy import *
    if sys.version_info[0] >= 3:
        from builtins import bool, int, float, complex, object, str
        unicode = str
    else:
        from __builtin__ import bool, int, float, complex, object, unicode, str

    from .core import round, abs, max, min

    __all__.extend(['__version__', 'pkgload', 'PackageLoader',
               'show_config'])
    __all__.extend(core.__all__)
    __all__.extend(_mat.__all__)
    __all__.extend(lib.__all__)
    __all__.extend(['linalg', 'fft', 'random', 'ctypeslib', 'ma'])

    # Filter out Cython harmless warnings
    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
    warnings.filterwarnings("ignore", message="numpy.ndarray size changed")

    # oldnumeric and numarray were removed in 1.9. In case some packages import
    # but do not use them, we define them here for backward compatibility.
    oldnumeric = 'removed'
    numarray = 'removed'

    # We don't actually use this ourselves anymore, but I'm not 100% sure that
    # no-one else in the world is using it (though I hope not)
    from .testing import Tester

    # Pytest testing
    from numpy.testing._private.pytesttester import PytestTester
    test = PytestTester(__name__)
    del PytestTester


    def _sanity_check():
        """
        Quick sanity checks for common bugs caused by environment.
        There are some cases e.g. with wrong BLAS ABI that cause wrong
        results under specific runtime conditions that are not necessarily
        achieved during test suite runs, and it is useful to catch those early.

        See https://github.com/numpy/numpy/issues/8577 and other
        similar bug reports.

        """
        try:
            x = ones(2, dtype=float32)
            if not abs(x.dot(x) - 2.0) < 1e-5:
                raise AssertionError()
        except AssertionError:
            msg = ("The current Numpy installation ({!r}) fails to "
                   "pass simple sanity checks. This can be caused for example "
                   "by incorrect BLAS library being linked in.")
            raise RuntimeError(msg.format(__file__))

    _sanity_check()
    del _sanity_check
