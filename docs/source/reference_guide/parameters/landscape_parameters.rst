.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


====================
Landscape parameters
====================

----
Main
----

------------------------------------------------------------------------------

**dim**

.. code-block:: python

              # x,y (a.k.a. j,i) dimensions of the Landscape
              'dim':                      (20,20),

:py:`tuple`

default: :py:`(20,20)`

reset: P
  
  This defines the y,x dimensions of the :py:`Landscape`,
  in units of cells. As you might imagine, these values are used 
  for a wide variety of basic operations throughout Geonomics. Change the
  default value to the dimensions of the landscape you wish to simulate on.


------------------------------------------------------------------------------

**res**

.. code-block:: python

              # x,y resolution of the Landscape
              'res':                      (1,1),

:py:`tuple`
  
default: :py:`(1,1)`

reset: N

  This defines the :py:`Landscape` resolution (or cell-size) in the y,x
  dimensions (matching the convention of the **dim** parameter).
  This information is only used if GIS rasters of :py:`Landscape` 
  layers are to be written out as GIS raster files (as parameterized in the
  'Data' parameters). Defaults to the meaningless value (1,1), and this value
  generally needn't be changed in your parameters file, because it will 
  be automatically updated to the resolution of any GIS rasters that 
  are read in for use as :py:`Layers` (assuming they all share the same
  resolution; otherwise, an Error is thrown). 


------------------------------------------------------------------------------

**ulc**

.. code-block:: python

              # x,y coords of upper-left corner of the Landscape
              'ulc':                      (0,0),

:py:`tuple`

default: :py:`(0,0)`

reset: N

  This defines the x,y upper-left corner (ULC) of the 
  :py:`Landscape` (in the units of
  some real-world coordinate reference system, e.g. decimal degrees, or
  meters). This information is only used if GIS rasters of 
  :py:`Landscape` layers are to be written out as GIS raster files. 
  Defaults to the meaningless value
  (0,0), and this value usually needn't be changed in your parameters file,
  because it will be automatically updated to match the ULC value 
  of any GIS rasters that are read in for use as :py:`Layers` (assuming 
  they all share the same ULC; otherwise, an Error is thrown).

        
------------------------------------------------------------------------------

**prj**

.. code-block:: python
              
              #projection of the Landscape
              'prj':                      None,

:py:`str`; (WKT projection string)

default: :py:`None`

reset: N

  This defines the projection of the :py:`Landscape`, as a
  string of Well Known Text (WKT). 
  This information is only used if GIS rasters of :py:`Landscape` layers are
  to be written out as GIS raster files. Defaults to :py:`None`, which is fine,
  because this value will be automatically updated to match the projection
  of any GIS rasters that are read in for us as :py:`Layers` (assuming they
  all share the same projection; otherwise, an Error is thrown)



------
Layers
------

------------------------------------------------------------------------------

**layer_<n>**

.. code-block:: python
     
      ################
      #### layers ####
      ################
          'layers': {
              #layer name (LAYER NAMES MUST BE UNIQUE!) 
              'layer_0': {

{:py:`str`, :py:`int`}

default: :py:`layer_<n>` 

reset? P

This parameter defines the name for each :py:`Layer`. (Note that unlike most
parameters, this parameter is a :py:`dict` key,
the value for which is a :py:`dict`
of parameters defining the :py:`Layer` being named.) As the capitalized
reminder in the parameters states, each :py:`Layer` must have a unique name
(so that a parameterized :py:`Layer` isn't overwritten in the
:py:`ParametersDict` by a second, identically-named :py:`Layer`; Geonomics
checks for unique names and throws an Error if this condition is not met.
:py:`Layer` names can, but needn't be, descriptive of what each 
:py:`Layer` represents. Example valid values include: 0, 0.1, 'layer_0', 1994,
'1994', 'mean_ann_tmp'. Names default to :py:`layer_<n>`,
where n is a series of integers starting from 0 and counting the number
of :py:`Layer`\s.



^^^^
Init
^^^^

There are four different types of :py:`Layers` that can be created. The
parameters for each are explained in the next four subsections.

""""""
random
""""""

------------------------------------------------------------------------------

**n_pts**

.. code-block:: python
    
                      #parameters for a 'random'-type Layer
                      'rand': {
                          #number of random points
                          'n_pts':                        500,

:py:`int`

default: 500

reset? P

This defines the number of randomly located, randomly valued points
from which the random :py:`Layer` will be interpolated. (Locations drawn
from uniform distributions between 0 and the :py:`Landscape` dimensions on
each axis. Values drawn from a uniform distribution between 0 and 1.)


------------------------------------------------------------------------------

**interp_method**

.. code-block:: python

                          #interpolation method ('linear', 'cubic', or 'nearest')
                          'interp_method':                'linear',
                          },

{:py:`'linear'`, :py:`'cubic'`, :py:`'nearest'`}

default: :py:`'linear'`

reset? N

This defines the method to use to interpolate random points to the array that
will serve as the :py:`Layer`'s raster. Whichever of the three valid values
is chosen (:py:`'linear'`, :py:`'cubic'`, or :py:`'nearest'`) will be passed
on as an argument to :py:`scipy.interpolate.griddata`. Note that the
:py:`'nearest'` method will generate a random categorical array, such as
might be used for modeling habitat types.


"""""""
defined
"""""""

------------------------------------------------------------------------------

**rast**

.. code-block:: python
   
                      #parameters for a 'defined'-type Layer 
                      'defined': {
                          #raster to use for the Layer
                          'rast':                    np.ones((100,100)),

nx2 :py:`np.ndarray`

default: :py:`np.ones((100,100))`

reset? Y

This defines the raster that will be used for this :py:`Layer`. Can be set to
:py:`None` if an array for the raster should instead be interpolated from a
set of valued points using the **pts**, **vals**, and **interp_method**
parameters. Dimensions of this array must match the dimensions of the
:py:`Landscape`.


------------------------------------------------------------------------------

**pts**

.. code-block:: python
   
                      #parameters for a 'defined'-type Layer 
                      'defined': {
                          #point coordinates
                          'pts':                    None,

nx2 :py:`np.ndarray`

default: :py:`None`

reset? Y

This defines the coordinates of the points to use to
interpolate this :py:`Layer`. Can be left as :py:`None` if the **rast**
parameter is given a :py:`numpy.ndarray`.


------------------------------------------------------------------------------

**vals**

.. code-block:: python

                           #point values
                           'vals':                  None,

{:py:`list`, 1xn :py:`np.ndarray`}

default: :py:`None`

reset? Y

This defines the values of the points to use to 
interpolate this :py:`Layer`. Can be left as :py:`None` if the **rast**
parameter is given a :py:`numpy.ndarray`.



------------------------------------------------------------------------------

**interp_method**

.. code-block:: python

                          #interpolation method {None, 'linear', 'cubic',
                          #'nearest'}
                          'interp_method':                None,
                          },

{:py:`'linear'`, :py:`'cubic'`, :py:`'nearest'`}

default: :py:`None`

reset? N

This defines the method to use to interpolate random points to the array that
will serve as the :py:`Layer`'s raster. Whichever of the valid string values
is chosen (:py:`'linear'`, :py:`'cubic'`, or :py:`'nearest'`) will be passed
on as an argument to :py:`scipy.interpolate.griddata`. Note that the
:py:`'nearest'` method will generate a random categorical array, such as
might be used for modeling habitat types. Can be left as :py:`None` if
the **rast** parameter is given a :py:`numpy.ndarray`.



""""
file
""""

------------------------------------------------------------------------------

**filepath**

.. code-block:: python
  
                      #parameters for a 'file'-type Layer 
                      'file': {
                          #</path/to/file>.<ext>
                          'filepath':                     '/PATH/TO/FILE.EXT',

:py:`str`

default: :py:`'/PATH/TO/FILE.EXT'`

reset? Y

This defines the location and name of the file that should be read in as the
raster-array for this :py:`Layer`. Valid file types include a '.txt' file
containing a 2d :py:`np.ndarray`, or any GIS raster file that can be read
by :py:`osgeo.gdal.Open`. In all cases, the raster-array read in from the
file must have dimensions equal to the stipulated dimensions of the
:py:`Landscape` (as defined in the **dims** parameter, above); otherwise,
Geonomics will throw an Error. Defaults to a dummy filename that must be
changed.


------------------------------------------------------------------------------

**scale_min_val**

.. code-block:: python

                          #minimum value to use to rescale the Layer to [0,1]
                          'scale_min_val':                None,

{:py:`float`, :py:`int`}

default: :py:`None`

reset? P

This defines the minimum value (in the units of the variable represented by
the file you are reading in) to use when rescaling the file's array to
values between 0 and 1. (This is done to satisfy the requirement that all
Geonomics :py:`Layer`\s have arrays in that interval). Defaults to :py:`None`
(in which case Geonomics will set it to the minimum value observed in this
file's array). But note that you should put good thought into
this parameter, because it *won't* necessarily be the minimum value
observed in the file; for example, if this file is being used
to create a :py:`Layer` that will undergo environmental change
in your `Model`, causing its real-world values to drop
below this file's minimum value, then you will probably want to set
this value to the minimum real-world value that will occur for this :py:`Layer`
during your :py:`Model` scenario, so that low values
that later arise on this `Layer` don't get truncated at 0.


------------------------------------------------------------------------------

**scale_max_val**

.. code-block:: python

                          #maximum value to use to rescale the Layer to [0,1]
                          'scale_max_val':                None,

{:py:`float`, :py:`int`}

default: :py:`None`

reset? P

This defines the maximum value (in the units of the variable represented by
the file you are reading in) to use when rescaling the file's array to
values between 0 and 1. (This is done to satisfy the requirement that all
Geonomics :py:`Layer`\s have arrays in that interval). Defaults to :py:`None`
(in which case Geonomics will set it to the maximum value observed in this
file's array). But note that you should put good thought into
this parameter, because it *won't* necessarily be the maximum value
observed in the file; for example, if this file is being used
to create a :py:`Layer` that will undergo environmental change
in your `Model`, causing its real-world values to increase
above this file's maximum value, then you will probably want to set
this value to the maximum real-world value that will occur for this 
:py:`Layer` during your :py:`Model` scenario, so that high values that 
later arise on this `Layer` don't get truncated at 1.


------------------------------------------------------------------------------

**coord_prec**

.. code-block:: python

                          #decimal-precision to use for coord-units (ulc & res)
                          'coord_prec':                5,

:py:`int`

default: 5

reset? P

This defines number of decimals to which to round upper-left corner
coordinates and resolution values read in from a raster file.
Because Geonomics requires equality of these values amongst all
input raster files, this allows the user to stipulate
the level of precision of their coordinate system, avoiding
false coordinate-system mismatch errors because of
arbitrary float imprecision.
(Note that for :py:`Layer`\s for which change rasters will be read in,
the same coordinate precision value will be used for all input rasters.)


------------------------------------------------------------------------------

**units**

.. code-block:: python

                          #units of this file's variable
                          'units':                       None,

{:py:`str`, :py:`None`}

default: None

reset? P

This is an optional parameter providing a string-representation
of the units in which a raster file's variable is expressed.
If provided, it will be used to label the colorbar on plots
of the raster's :py:`Layer`.


"""""
nlmpy
"""""

------------------------------------------------------------------------------

**function**

.. code-block:: python

                      #parameters for an 'nlmpy'-type Layer
                      'nlmpy': {
                          #nlmpy function to use the create this Layer
                          'function':                 'mpd',

:py:`str` that is the name of an :py:`nlmpy` function

default: :py:`'mpd'`

reset? P

This indicates the :py:`nlmpy` function that should be used to generate
this :py:`Layer`'s array. (:py:`nlmpy` is a Python package for
generating neutral landscape models; NLMs.) Defaults to :py:`'mpd'` (the
function for creating a midpoint-displacement NLM). Can be set to any other
:py:`str` that identifies a valid :py:`nlmpy` function, but then the
remaining parameters in this section must be changed to the parameters
that that function needs, and *only* those parameters 
(because they will be unpacked into this function,
i.e. passed on to it, at the time it is called.
(Visit the `Cheese Shop <https://pypi.org/project/nlmpy/>`_ for more 
information about the :py:`nlmpy` package and available functions).


------------------------------------------------------------------------------

**nRow**

.. code-block:: python

                          #number of rows (MUST EQUAL LAND DIMENSION y!)
                          'nRow':                     20,


:py:`int`

default: 20

reset? P

This defines the number of rows in the :py:`nlmpy` array that is created.
As the capitalized reminder in the parameters file mentions, this must be
equal to the y-dimension of the :py:`Landscape`; otherwise, an error
will be thrown. Note that this parameter (as for the remaining parameters in
this section, other than the **function** parameter) is valid for the
default :py:`nlmpy.mpd` function that is set by the
**function** parameter); if you are using a different :py:`nlmpy`
function to create this :py:`Layer` then this and the remaining parameters
must be changed to the parameters that that function needs, 
and *only* those parameters (because they will be unpacked into that function,
i.e. passed on to it, at the time it is called).


------------------------------------------------------------------------------

**nCol**

.. code-block:: python

                          #number of cols (MUST EQUAL LAND DIMENSION x!)
                          'nCol':                     20,


:py:`int`

default: 20

reset? P

This defines the number of columns in the :py:`nlmpy` array that is created.
As the capitalized reminder in the parameters file mentions, this must be
equal to the x-dimension of the :py:`Landscape`; otherwise, an error
will be thrown. Note that this parameter (as for the remaining parameters in
this section, other than the **function** parameter) is valid for the
default :py:`nlmpy.mpd` function that is set by the
**function** parameter); if you are using a different :py:`nlmpy`
function to create this :py:`Layer` then this and the remaining parameters
must be changed to the parameters that that function needs, 
and *only* those parameters (because they will be unpacked into that function,
i.e. passed on to it, at the time it is called).


------------------------------------------------------------------------------

**h**

.. code-block:: python

                          #level of spatial autocorrelation in element values
                          'h':                     1,


:py:`float`

default: 1

reset? P

This defines the level of spatial autocorrelation in the element values
of the :py:`nlmpy` array that is created.
Note that this parameter (and the remaining parameters in
this section, other than the **function** parameter) is valid for the
default :py:`nlmpy` function (:py:`nlmpy.mpd`, which is set by the
**function** parameter); but if you are using a different :py:`nlmpy`
function to create this :py:`Layer` then this and the remaining parameters
must be changed to the parameters that that function needs, 
and *only* those parameters (because they will be unpacked into that function,
i.e. passed on to it, at the time it is called).


^^^^^^
Change
^^^^^^

------------------------------------------------------------------------------

**change_rast**

.. code-block:: python

                  #land-change event for this Layer
                  'change': {
                      #array of file for final raster of event, or directory
                      #of files for each stepwise change in event
                      'change_rast':         '/PATH/TO/FILE.EXT',

{2d :py:`np.ndarray`, :py:`str`}

default: :py:`'/PATH/TO/FILE.EXT'`

reset? Y

This defines either the final raster of the :py:`Landscape` change event
(with valid values being a :py:`numpy.ndarray` or a string pointing
to a valid raster file, i.e. a file that can be read by :py:`osgeo.gdal.Open`);
or the stepwise series of changes to be made over the course of the
:py:`Landscape` change event (with the valid value being a string
pointing to a directory full of valid raster files).
Note that whether an array, a raster, or multiple rasters
are input, their dimensions must be equal to the dimensions of the :py:`Layer`
that is being changed (and hence to the :py:`Landscape` to which it belongs).
Also note that if a directory of stepwise-change rasters is provided, the
rasters' filenames must begin with the integer timesteps at which they
should be used during the change event, followed by underscores. (For example,
files with the filenames '50_mat_2001.tif', '60_mat_2011.tif',
'65_mat_2011.tif' would be used at timesteps 50, 60, and 65 during a model.)
Defaults to a dummy file name that must be changed.


------------------------------------------------------------------------------

**start_t**

.. code-block:: python

                   #starting timestep of event
                   'start_t':          50,

:py:`int`

default: 50

reset? P

This indicates the first timestep of the :py:`Landscape`-change event. 
Defaults to 50, but should be set to suit your specific scenario. 
If a directory of files is provided for the **change_rast** parameter,
then this must match the earliest timestep in that series of files
(as indicated by the integers at the beginning of the file names).


------------------------------------------------------------------------------

**end_t**

.. code-block:: python

                   #ending timestep of event
                   'end_t':          100,

:py:`int`

default: 100

reset? P

This indicates the last timestep of the
:py:`Landscape`-change event.
Defaults to 100, but should be set to suit your specific scenario.
If a directory of files is provided for the **change_rast** parameter,
then this must match the final timestep in that series of files
(as indicated by the integers at the beginning of the file names).


------------------------------------------------------------------------------

**n_steps**

.. code-block:: python

                   #number of stepwise changes in event
                   'n_steps':          5,

:py:`int`

default: 5

reset? P

This indicates the number of stepwise changes to use to model a
:py:`Landscape`-change event.
If the the **change_rast** parameter is a directory of files, 
then the value of this parameter must be the number of files in that directory. 
If the **change_rast** parameter is either an :py:`np.ndarray` or a file name,
then the changes during the :py:`Landscape`-change event
are linearly interpolated (cellwise for the whole :py:`Layer`) to this
number of discrete, instantaneous :py:`Landscape` changes between
the starting and ending rasters. Thus, the fewer the number of 
steps, the larger, magnitudinally, each change will be. So more
steps may be 'better', as it will better approximate change that is continuous
in time. However, there is a potenitally significant memory trade-off here:
The whole series of stepwise-changed arrays is computed when the
:py:`Model` is created, then saved and used at the appropriate timestep
during each :py:`Model` run (and if the :py:`Layer` that is changing is used
by any :py:`Species` as a :py:`_ConductanceSurface` then each 
intermediate :py:`_ConductanceSurface` is also calculated
when the :py:`Model` is first built, which can be much more memory-intensive
because these are 3-dimensional arrays).
These objects take up memory, which may be limiting for larger
:py:`Model`\s and/or :py:`Landscape` objects. This often will not be a
major issue, but depending on your use case it could pose a problem, so
is worth considering.
