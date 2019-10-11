# flake8: noqa

from . import ops
from . import sim
from . import utils
from . import help
from .main import *

__docformat__ = 'restructuredtext'

#package-level docstring
__doc__ = """
geonomics - a Python package for running complex landscape genomics simulations
===============================================================================

"""
with open('./geonomics/README.rst', 'r') as f:
    __doc__ = __doc__ + f.read()
