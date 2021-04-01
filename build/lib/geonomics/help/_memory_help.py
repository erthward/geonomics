#!/usr/bin/python
#_memory_help.py

'''
Helper functions for assesing the memory usage of Geonomics objects.
'''

#------------------------------------
# CLASSES ---------------------------
#------------------------------------


#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------

#A nice function for recursively iterate and sum the size of an object
#Courtesy of Aaron Hall's response at
#http://stackoverflow.com/questions/449560/
#how-do-i-determine-the-size-of-an-object-in-python
def get_obj_size(obj_0):
    """
    Get the full size of a complex object

    Recursively iterate through the structure of some complex object,
    summing the size of all components, to return the full size of the object.

    Parameters
    ----------
    obj_0, object
        The object whose full size should be calculated.

    Returns
    -------
    out : float
        The full calculated size of the the input object.
    """

    import sys
    from numbers import Number
    from collections import Set, Mapping, deque

    try: # Python 2
        zero_depth_bases = (basestring, Number, xrange, bytearray)
        iteritems = 'iteritems'
    except NameError: # Python 3
        zero_depth_bases = (str, bytes, Number, range, bytearray)
        iteritems = 'items'

    def inner(obj, _seen_ids = set()):
        obj_id = id(obj)
        if obj_id in _seen_ids:
            return 0
        _seen_ids.add(obj_id)
        size = sys.getsizeof(obj)
        if isinstance(obj, zero_depth_bases):
            pass # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, Set, deque)):
            size += sum(inner(i) for i in obj)
        elif isinstance(obj, Mapping) or hasattr(obj, iteritems):
            size += sum(inner(k) + inner(v) for k, v in getattr(obj,
                iteritems)())
        # Check for custom object instances - may subclass above too
        if hasattr(obj, '__dict__'):
            size += inner(vars(obj))
        if hasattr(obj, '__slots__'): # can have __slots__ with __dict__
            size += sum(inner(getattr(obj,
                s)) for s in obj.__slots__ if hasattr(obj, s))
        return size
    return inner(obj_0)

