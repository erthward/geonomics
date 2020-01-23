.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash

   
================
Other parameters
================

----
Main
----


-------------------------------------------------------------------------------

**T**

.. code-block:: python

          #total Model runtime (in timesteps)
          'T':            100,

:py:`int`

default: 100

reset? Y

This indicates the total number of timesteps for which the main portion of
a :py:`Model` (i.e. the portion after the burn-in has completed) will be run
during each iteration.


------------------------------------------------------------------------------

**burn_T**

.. code-block:: python

          #min burn-in runtime (in timesteps)
          'burn_T':       30,

:py:`int`

default: 30

reset? P

This indicates the minimum number of timesteps for which a :py:`Model`'s
burn-in will run. (Note this is only a minimum because the test for
burn-in completion includes a check that at least this many timesteps have
elapsed, but also includes two statistical checks of stationarity of the
size of each :py:`Species` in a :py:`Community`.)


------------------------------------------------------------------------------

**seed**

.. code-block:: python

          #seed number
          'num':          None,

{:py:`int`, :py:`None`}

default: :py:`None`

reset? P
       
This indicates whether or not to set the seeds of the random number
generators (by calling :py:`np.random.seed` and :py:`random.seed`)
before building and running a :py:`Model`. If value is an integer, the seeds
will be set to that value. If value is :py:`None`, seeds will not be set.


----------
Iterations
----------

-------------------------------------------------------------------------------

**num_iterations**

.. code-block:: python

              #num iterations
              'n_its': 2,

:py:`int`

default: 2

reset? Y

This indicates the number of iterations for which the :py:`Model`
should be run. (Note that for each iteration a separate subdirectory of
data and stats will be written, if your :py:`Model` has parameterized data
and stats to be collected.)


------------------------------------------------------------------------------

**rand_landscape**

.. code-block:: python

              #whether to randomize Landscape each iteration
              'rand_landscape':    False,

:py:`bool`

default: False

reset? P

This indicates whether the :py:`Landscape` should be randomized for each
iteration. If True, a new :py:`Landscape` will be generated at the start
of each iteration. If False, the :py:`Landscape` from iteration 0 will be
saved and reused for each subsequent iteration.


------------------------------------------------------------------------------

**rand_community**

.. code-block:: python

              #whether to randomize Community each iteration
              'rand_comm':    False,

:py:`bool`

default: False

reset? P

This indicates whether the :py:`Community` should be randomized for each
iteration. If True, a new :py:`Community` will be generated at the start
of each iteration. If False, the :py:`Community` from iteration 0 will be
saved and reused for each subsequent iteration (and whether that
:py:`Community` is saved before or after being burned in will depend on
the value provided to the **repeat_burn** parameter).


------------------------------------------------------------------------------

**repeat_burn**

.. code-block:: python

              #whether to burn in each iteration
              'repeat_burn':  False,

:py:`bool`

default: False

reset? P

This indicates whether a reused :py:`Community` should be burned in
separately for each iteration for which it is reused. If True, the
:py:`Community` from iteration 0 will be saved as soon as its instantiated,
but will have a new burn-in run for each iteration in which it is used. If
False, the :py:`Community` from iteration 0 will be saved after its burn-in
is complete, and then will only have the main portion of its :py:`Model` run
separately during each iteration. (Note that if **rand_community** is set to True then
the value of this parameter will not be used.)


^^^^
Data
^^^^

""""""""
Sampling
""""""""

-------------------------------------------------------------------------------

**scheme**

.. code-block:: python

                  #sampling scheme {'all', 'random', 'point', 'transect'}
                  'scheme':               'random',

{:py:`'all'`, :py:`'random'`, :py:`'point'`, :py:`'transect'`}

default: :py:`'random'`

reset? P

This indicates the sampling scheme to use when collecting data from a
:py:`Model`. Currently valid values include :py:`'all'`, 
:py:`'random'`, :py:`'point'`, and :py:`'transect'`.

With :py:`'all'`, data will be collected for all :py:`Individual`\s
at each sampling timestep. With :py:`'random'`, data will be collected from a
random sample of :py:`Individual`\s (of size indicated by parameter **n**) 
from anywhere on the :py:`Landscape`.
With :py:`'point'`, data will be collected from random samples of size **n**
within a certain distance (**radius**) of each of a set of particular
points (**points**). With :py:`'transect'`, a linear transect of some
number of points (**n_transect_points**) between some endpoints
(**transect_endpoints**) will be created, and then data will be collected
from random samples of size **n** with a certain distance (**radius**)
of each point along the transect.


------------------------------------------------------------------------------

**n**

.. code-block:: python

                  #sample size at each point, for point & transect sampling
                  'n':                    250,

:py:`int`

default: 250

reset? P

This indicates the total number of :py:`Individual`\s to sample each time
data is collected (if **scheme** is :py:`'random'`), or the number of 
:py:`Individual`\s to sample around each one of a set of points (if **scheme**
is :py:`'point'` or :py:`'transect'`). This parameter will only be used if
**scheme** is :py:`'random'`, :py:`'point'`, or :py:`'transect'`; otherwise
it may be set to :py:`None`.


------------------------------------------------------------------------------

**points**

.. code-block:: python

                  #coords of collection points, for point sampling
                  'points':               None,

{:py:`tuple` of 2-:py:`tuple`\s, :py:`None`}

default: :py:`None`

reset? P

This indicates the points around which to sample :py:`Individual`\s for data
collection. This parameter will only be used if **scheme** is :py:`'point'`;
otherwise it may be set to :py:`None`.


------------------------------------------------------------------------------

**transect_endpoints**

.. code-block:: python

                  #coords of transect endpoints, for transect sampling
                  'transect_endpoints':   None,

{2-:py:`tuple` of 2-:py:`tuple`\s, :py:`None`}

default: :py:`None`

reset? P

This indicates the endpoints between which to create a transect, along which
:py:`Individual`\s will be sampled for data collection. 
This parameter will only be used if **scheme** is :py:`'transect'`; 
otherwise it may be set to :py:`None`.


------------------------------------------------------------------------------

**n_transect_points**

.. code-block:: python

                  #num points along transect, for transect sampling
                  'n_transect_points':    None,

{:py:`int`, :py:`None`}

default: :py:`None`

reset? P

This indicates the number of points to create on the transect along which
:py:`Individual`\s will be sampled for data collection. 
This parameter will only be used if **scheme** is :py:`'transect'`; 
otherwise it may be set to :py:`None`.


------------------------------------------------------------------------------

**radius**

.. code-block:: python

                  #collection radius around points, for point & transect sampling
                  'radius':               None,

{:py:`float`, :py:`int`, :py:`None`}

default: :py:`None`

reset? P

This indicates the radius around sampling points within which
:py:`Individual`\s may be sampled for data collection. 
This parameter will only be used if **scheme** is :py:`'point'` or 
:py:`'transect'`; otherwise it may be set to :py:`None`.


------------------------------------------------------------------------------

**when**

.. code-block:: python

                  #when to collect data
                  'when':                 None,

{:py:`int`, :py:`list` of :py:`int`\s, :py:`None`}

default: :py:`None`

reset? P

This indicates the timesteps during main :py:`Model` iterations
at which data should be collected (in addition to after the final timestep
of each iteration, when data is always collected for any :py:`Model` for which
data collection is parameterized). If value is a non-zero :py:`int`,
it will be treated as a frequency at which data should be collected (e.g.
a value of 5 will cause data to be collected every 5 timesteps). If value
is a list of :py:`int`\s, they will be treated as the particular timesteps
at which data should be collected. If value is 0 or :py:`None`, 
data will be collected only after the final timestep.


------------------------------------------------------------------------------

**include_landscape**

.. code-block:: python

                  #whether to save current Layers when data is collected
                  'include_landscape':         False,

:py:`bool`

default: False

reset? P

This indicates whether to include the :py:`Landscape` :py:`Layer`\s among the
data that is collected. If True, each :py:`Layer` will be written to a raster
or array file (according to the format indicated by
**geo_rast_format**) each time data is collected.


------------------------------------------------------------------------------

**include_fixed_sites**

.. code-block:: python

                  #whether to include fixed loci in VCF files
                  'include_fixed_sites':  False,

:py:`bool`

default: False

reset? P

This indicates whether fixed sites (i.e. loci which are fixed for either the
0 or 1 allele) should be included in any VCF files that are written. Thus,
this parameter is only relevant if :py:`'vcf'` is one of the genetic data
formats indicated by **gen_format**.



""""""
Format
""""""

-------------------------------------------------------------------------------

**gen_format**

.. code-block:: python

                  #format for genetic data {'vcf', 'fasta'}
                  'gen_format':           ['vcf', 'fasta'],

{:py:`'vcf'`, :py:`'fasta'`, :py:`['vcf', 'fasta']`}

default: :py:`['vcf', 'fasta']`

reset? P

This indicates the format or formats to use for writing genetic data.
data. Currently valid formats include :py:`'vcf'` and :py:`'fasta'` formats.
Either or both formats may be specified; all formats that are specified will
be written each time data is collected.


------------------------------------------------------------------------------

**geo_vect_format**

.. code-block:: python

                  #format for vector geodata {'csv', 'shapefile', 'geojson'}
                  'geo_vect_format':      'csv',

{:py:`'csv'`, :py:`'shapefile'`, :py:`'geojson'`}

default: :py:`'csv'`

reset? P

This indicates the format to use for writing geographic vector data (i.e.
:py:`Individual`\s' point locations). 
Currently valid formats include :py:`'csv'`, :py:`'shapefile'`,
and :py:`'geojson'`. Any one format may be specified.


------------------------------------------------------------------------------

**geo_rast_format**

.. code-block:: python

                  #format for raster geodata {'geotiff', 'txt'}
                  'geo_rast_format':      'geotiff',

{:py:`'geotiff'`, :py:`'txt'`}

default: :py:`'geotiff'`

reset? P

This indicates the format to use for writing geographic raster data (i.e.
:py:`Layer` arrays). Currently valid formats include :py:`'geotiff'`
and :py:`'txt'`. Either format may be specified. Note that this parameter
will only be used if the **include_landscape** parameter is set to True.


^^^^^
Stats
^^^^^

-------------------------------------------------------------------------------

The stats parameters section has subsection for each statistc that Geonomics
can calculate. (Currently valid statistics include:
  - *'Nt'*: number of individuals at timestep t
  - *'het'*: heterozygosity
  - *'maf'*: minor allele frequency
  - *'mean_fit'*: mean fitness of a :py:`Species`
  - *'ld'*: linkage disequilibrium
  
  There are only a few parameters, which are shared across all
of those subsections, and each parameter always means the same thing despite
which statistic it is parameterizing. Thus, hereafter we provide a single of
each of those parameters are how it works, regardless of the statistic for
which it used:


------------------------------------------------------------------------------

**calc**

.. code-block:: python

                #whether to calculate
                'calc':     True,


:py:`bool`

default: (varies by statistic)

reset? P

This indicates whether or not a given statistic should be calculated. Thus,
only those statistics whose **calc** parameters are set to True will be
calculated and saved when their :py:`Model` is run.


------------------------------------------------------------------------------

**freq**

.. code-block:: python

                #calculation frequency (in timesteps)
                'freq':     5,

:py:`int`

default: (varies by statistic)

reset? P

This indicates the frequency with which a given statistic should be calculated
during each iteration (in timesteps). If set to 0, Geonomics will calculate
and save this statistic for only the first and last timesteps
of each iteration.


------------------------------------------------------------------------------

**mean**

.. code-block:: python

                #whether to mean across sampled individs
                'mean': False,

:py:`bool`

default: (varies by statistic, and only valid for certain statistics)

reset? P

For some statistics that produce a vector of values each timestep
when they are collected (containing one value per :py:`Individual`),
such as heterozygosity, this indicates
whether those values should instead be meaned and saved as a
single value for each timestep.
