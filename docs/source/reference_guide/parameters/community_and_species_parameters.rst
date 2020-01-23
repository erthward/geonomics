.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


================================
Community and Species parameters
================================

-------
Species
-------


------------------------------------------------------------------------------

**spp_<n>**

.. code-block:: python
 
              #spp name (SPECIES NAMES MUST BE UNIQUE!) 
              'spp_0' :   {

{:py:`str`, :py:`int`}

default: :py:`spp_<n>` 

reset? P

This parameter defines the name for each :py:`Species`.
(Note that unlike most parameters, this parameter is 
a :py:`dict` key, the value for which is a :py:`dict`
of parameters defining the :py:`Species` being named.) As the capitalized
reminder in the parameters states, each :py:`Species`
must have a unique name (so that a parameterized 
:py:`Species` isn't overwritten in the :py:`ParametersDict` by a
second, identically-named :py:`Species`; Geonomics
checks for unique names and throws an Error if this condition is not met.
:py:`Species` names can, but needn't be, descriptive of what each 
:py:`Species` represents. Example valid values include: 0, 'spp0',
'high-dispersal', 'C. fasciata'. Names default to 
:py:`spp_<n>`, where n is a series of
integers starting from 0 and counting the number of :py:`Species`.

^^^^
Init
^^^^

------------------------------------------------------------------------------

**N**

.. code-block:: python
  
                  'init': {
                      #starting number of individs
                      'N':                250,

:py:`int`

default: 250

reset? P

This defines the starting size of this :py:`Species`. Importantly, this
may or may not be near the stationary size of the :py:`Species` after
the :py:`Model` has burned in, because that size will depend on the
carrying-capacity raster (set by the **K** parameter), and on
the dynamics of specific a :py:`Model` (because of the interaction of
its various parameters).


------------------------------------------------------------------------------

**K_layer**

.. code-block:: python

                      #name of the carrying-capacity Layer
                      'K_layer':         'layer_0',

:py:`str`

default: 'layer_0'

reset? P

This indicates, by name, the :py:`Layer` to be used as the
carrying-capacity raster for a :py:`Species`. The values of this
:py:`Layer`, multiplied by **K_factor**, should express
the carrying capacity at each cell, in number
of :py:`Individual`\s. Note that the sum of the values of the product of
this :py:`Layer` and **K_factor**
can serve as a rough estimate of the expected stationary 
number of individuals of a :py:`Species`; 
however, observed stationary size could vary
substantially depending on various other :py:`Model` parameters (e.g. birth
and death rates and mean number of offspring per mating event) as well
as on stochastic events (e.g. failure to colonize, or survive in, all
habitable portions of the :py:`Landscape`).


------------------------------------------------------------------------------

**K_factor**

.. code-block:: python

                      #multiplicative factor for carrying-capacity layer
                      'K_factor':         1,

{:py:`int`, :py:`float`}

default: 1

reset? P

This defines the factor by which the raster of the :py:`Layer` indicated
by **K_layer** will be multiplied to create a :py:`Species`' carrying-
capacity raster. Because :py:`Layer`\s' rasters are constrained to [0,1],
this allows the user to stipulate that cells have carrying capacities in
excess of 1.


^^^^^^
Mating
^^^^^^

------------------------------------------------------------------------------

**repro_age**

.. code-block:: python

                  'mating'    : {
                      #age(s) at sexual maturity (if tuple, female first)
                      'repro_age':            0,

{:py:`int`, :py:`(int, int)`, :py:`None`}

default: 0

reset? P

This defines the age at which :py:`Individual`\s in the :py:`Species`
can begin to reproduce. If the value provided is a 2-tuple of different
numbers (and the :py:`Species` uses separate sexes), then the first
number will be used as females' reproductive age, the second as males'.
If the value is 0, or :py:`None`, :py:`Individual`\s are capable
of reproduction from time of time.


------------------------------------------------------------------------------

**sex**

.. code-block:: python
        
                      #whether to assign sexes
                      'sex':                  False,

:py:`bool`

default: False

reset? P

This determines whether :py:`Individual`\s will be assigned separate sexes
that are used to ensure only male-female mating events.


------------------------------------------------------------------------------

**sex_ratio**

.. code-block:: python
                        
                      #ratio of males to females
                      'sex_ratio':            1/1,


{:py:`float`, :py:`int`}

default: 1/1

reset? P

This defines the ratio of males to females (i.e. it will be converted to
a probability that an offspring is a male, which is used as the probability
of a Bernoulli draw of that offspring's sex). 


------------------------------------------------------------------------------

**distweighted_birth**

.. code-block:: python

                      #whether P(birth) should be weighted by parental dist
                      'distweighted_birth':  False,


#NOTE: I WILL PROBABLY GET RID OF THIS PARAMETER...


------------------------------------------------------------------------------

**R**

.. code-block:: python

                      #intrinsic growth rate
                      'R':                    0.5,

:py:`float`

default: 0.5

reset? P

This defines a :py:`Species`' intrinsic growth rate, which is used
as the 'R' value in the spatialized logistic growth equation that
regulates population density (:math:`\frac{\mathrm{d}
N_{x,y}}{\mathrm{d}t}=rN_{x,y}(1-\frac{N_{x,y}}{K_{x,y}})`).


------------------------------------------------------------------------------

**b**

.. code-block:: python
                       
                      #intrinsic birth rate (MUST BE 0<=b<=1)
                      'b':                    0.2,

:py:`float` in interval [0, 1]

default: 0.2

reset? P

This defines a :py:`Species`' intrinsic birth rate, which is
implemented as the probability that an identified potential mating
pair successfully produces offspring. Because this is a probability, as
the capitalized reminder in the parameters file mentions, this value must
be in the inclusive interval [0, 1].

NOTE: this may later need to be re-implemented to allow for spatial
variation in intrinsic rate (i.e.. expression of a birth-rate raster),
and/or for density-dependent birth as well as mortality


------------------------------------------------------------------------------

**n_births_dist_lambda**

.. code-block:: python

                      #expectation of distr of n offspring per mating pair
                      'n_births_distr_lambda':      1,

{:py:`float`, :py:`int`}

default: 1

reset? P

This defines the lambda parameter for the Poisson distribution from 
which a mating pair's number of offspring is drawn (unless **n_births_fixed**
is set to True, in which case it defines the number of offspring 
produced by each successful mating event). Hence, this is either the
expected  or exact value for the number of offspring born in a
successful mating event (depending on how **n_births_fixed** is set).


------------------------------------------------------------------------------

**n_births_fixed**

.. code-block:: python

                      #whether n births should be fixed at n_births_dist_lambda
                      'n_births_fixed':           True,

:py:`bool`

default: True

reset? P

This determines whether or not the number of births for each mating event will
be fixed. If set to true, each successful mating event will produce
**n_births_distr_lambda** new offspring.


------------------------------------------------------------------------------

**mating_radius**

.. code-block:: python

                      #radius of mate-search area
                      'mating_radius':        1

{:py:`float`, :py:`int`}

default: 1

reset? Y

This defines the radius within which an :py:`Indvidual` can find a mate.
This radius is provided to queries run on the :py:`_KDTree` object.


^^^^^^^^^
Mortality
^^^^^^^^^

------------------------------------------------------------------------------

**max_age**

.. code-block:: python
                        
                      #maximum age
                      'max_age':              1,

{:py:`int`, :py:`None`}

default: 1

reset? P

This defines the maximum age an individual can achieve before being
forcibly culled from the :py:`Species`. Defaults to 1 (which will create
a Wright-Fisher-like simulation, with discrete generations). Can be set
to any other age, or can be set to :py:`None` (in which case no maxmimum
age is enforced).


------------------------------------------------------------------------------

**d_min**

.. code-block:: python
        
                      #min P(death) (MUST BE 0<=d_min<=1)
                      'd_min':                     0,

:py:`float` in interval [0, 1]

default: 0

reset? N

This defines the minimum probabilty of death that an :py:`Individual`
can face each time its Bernoulli death-decision is drawn. Because this 
is a probability, as the capitalized reminder in 
the parameters file mentions, this value must be in the 
inclusive interval [0, 1].

------------------------------------------------------------------------------

**d_max**

.. code-block:: python

                      #max P(death) (MUST BE 0<=d_max<=1)
                      'd_max':                    1,

:py:`float` in interval [0, 1]

default: 1

reset? N

This defines the minimum probabilty of death that an :py:`Individual`
can face each time its Bernoulli death-decision is drawn. Because this 
is a probability, as the capitalized reminder in 
the parameters file mentions, this value must be in the 
inclusive interval [0, 1].


------------------------------------------------------------------------------

**density_grid_window_width**


.. code-block:: python

                  'mortality'     : {
                      #width of window used to estimate local pop density
                      'dens_grid_window_width':   None,

{:py:`float`, :py:`int`, :py:`None`}

default: None

reset? N

This defines the width of the window used by the :py:`_DensityGridStack`
to estimate a raster of local :py:`Species` densities. The user should
feel free to set different values for this parameter (which could be
especially helpful when calling :py:`Model.plot_density` to inspect the
resulting surfaces calculated at different window widths, if trying
to heuristically choose a reasonable value to set for a
particular simulation scenario). But be aware that choosing particularly
small window widths (in our experience, windows smaller than ~1/20th of
the larger :py:`Landscape` dimension) will cause dramatic increases in the 
run-time of the density calculation (which runs twice per timestep).
Defaults to :py:`None`, which will internally be set to the integer
nearest to 1/10th of the larger :py:`Landscape` dimension; 
for many purposes this will work, but in some cases
the user may wish to control this.


^^^^^^^^
Movement
^^^^^^^^

------------------------------------------------------------------------------

**move**

.. code-block:: python

                     #whether or not the species is mobile
                     'move':                    True,

:py: `bool`

default: True

reset? P

This determines whether the :py: `Species` being parameterized is mobile
(i.e. whether its individuals should move). A :py:`Species` without movement
will still undergo dispersal of offspring, but after dispersing
those offspring will remain fixed in location until death.


------------------------------------------------------------------------------

**direction_distr_mu**

.. code-block:: python
 
                'movement': {
                     #mode of distr of movement direction
                     'direction_distr_mu':      1,

{:py:`int`, :py;`float`}

default: 1

reset? N

This is the :math:`\mu` parameter of the VonMises distribution
(a circularized normal distribution) from which
movement directions are chosen when movement is random and isotropic 
(rather than
being determined by a :py:`_ConductanceSurface`;
if a :py:`_ConductanceSurface`
is being usen this parameter is ignored). The :math:`\kappa` value
that is fed into this same distribution (**direction_distr_kappa**)
causes it to be very dispersed,
such that the distribution is effectively a uniform distribution on 
the unit circle (i.e. all directions are effectively equally probable).
For this reason, changing this parameter without changing the 
**direction_distr_kappa** value also, will make no change in the directions
drawn for movement.  If random, isotropic
movement is what you aim to model then there is probably little reason 
to change these parameters.


------------------------------------------------------------------------------

**direction_distr_kappa**

.. code-block:: python

                     #concentration of distr of movement direction
                     'direction_distr_kappa':  0,

{:py:`int`, :py:`float`}

default: 0

reset? N

This is the :math:`\kappa` parameter of the VonMises distribution
(a circularized normal distribution) from which
movement directions are chosen when movement is random and isotropic 
(rather than
being determined by a :py:`_ConductanceSurface`;
if a :py:`_ConductanceSurface`
is being usen this parameter is ignored). The default value of 0 will  
cause this distribution to be very dispersed, approximating a uniform
distribution on the unit circle and rendering the :math:`\mu`
value (**direction_distr_mu**) effectively meaningless. However, as this
parameter's value increases the resulting circular distributions will become
more concentrated around :math:`\mu`, making the value fed to
**direction_distr_mu** influential. If random, isotropic
movement is what you aim to model then there is probably little reason 
to change these parameters.

------------------------------------------------------------------------------


**distance_distr_mu**

.. code-block:: python

                     #mean of distr of movement distance
                     'distance_distr_mu':      0.5,

{:py:`int`, :py:`float`}

default: 0.5

reset? Y

This is the :math:`\mu` parameter of the Wald distribution used to draw
movement distances, expressed in units of raster cell widths
(or the wider of the two dimensions of a cell, in the case of a
non-square-resolution raster). 
This parameter and **distance_distr_sigma**
(the Wald distribution's :math:`sigma`) should be
set to reflect a distribution of movement distances that is appropriate
for your scenario.

------------------------------------------------------------------------------


**distance_distr_sigma**

.. code-block:: python

                     #variance of distr of movement distance
                     'distance_distr_sigma':   0.5,

{:py:`int`, :py:`float`}

default: 0.5 

reset? Y

This is the :math:`\sigma` parameter of the Wald distribution used to draw
movement distances, expressed in units of raster-cell widths
(or the wider of the two dimensions of a cell, in the case of a
non-square-resolution raster). 
This parameter and **distance_distr_mu**
(the Wald distribution's :math:`mu`) should be
set to reflect a distribution of movement distances that is appropriate
for your scenario.

------------------------------------------------------------------------------


**dispersal_distr_mu**

.. code-block:: python

                     #mean of distr of dispersal distance
                     'dispersal_distr_mu':     0.5,

{:py:`int`, :py:`float`}

default: 0.5

reset? Y

This is the :math:`\mu` parameter of the Wald distribution used to draw
dispersal distances, expressed in units of raster-cell widths
(or the wider of the two dimensions of a cell, in the case of a
non-square-resolution raster). 
This paramter and **distance_distr_sigma**
(the Wald distribution's :math:`sigma`) should be
set to reflect a distribution of dispersal distances that is appropriate
for your scenario.


------------------------------------------------------------------------------

**dispersal_distr_sigma**

.. code-block:: python

                     #variance of distr of dispersal distance
                     'dispersal_distr_sigma':  0.5,
                 
{:py:`int`, :py:`float`}

default: 0.5

reset? Y

This is the :math:`\sigma` parameter of the Wald distribution used to draw
dispersal distances, expressed in units of raster-cell widths
(or the wider of the two dimensions of a cell, in the case of a
non-square-resolution raster). 
This paramter and **distance_distr_mu**
(the Wald distribution's :math:`mu`) should be
set to reflect a distribution of dispersal distances that is appropriate
for your scenario.


"""""""""""""""""""""""""""""""""""""""""""
Movement and Dispersal _ConductanceSurfaces
"""""""""""""""""""""""""""""""""""""""""""

------------------------------------------------------------------------------

**layer**

.. code-block:: python

                     'move_surf'     : {
                         #move-surf Layer name
                         'layer':                'layer_0',

:py:`str`

default: :py:`'layer_0'`

reset? P

This indicates, by name, the :py:`Layer` to be used as to construct the
:py:`_ConductanceSurface` for a :py:`Species`. Note that this can also
be thought of as the :py:`Layer` that should serve as a
:py:`Species`' permeability raster (because :py:`Individual`\s moving
on this :py:`_ConductanceSurface` toward the higher
(if mixture distributions are used) or highest
(if unimodl distributions are used) values in their neighborhoods). 


------------------------------------------------------------------------------

**mixture**

.. code-block:: python

                         #whether to use mixture distrs
                         'mixture':              True,

:py:`bool`

default: True

reset? P

This indicates whether the :py:`_ConductanceSurface` should be built using
VonMises mixture distributions or unimodal VonMises distributions. 
If True, each cell in the :py:`_ConductanceSurface` will have an approximate
circular distribution that is a
weighted sum of 8 unimodal VonMises distributions (one per cell in the 8-cell
neighborhood); each of those summed unimodal distributions will have as its 
mode the direction of the neighboring cell on which it is based and as its 
weight the relative permeability of the cell on which it is based 
(relative to the full neighborhood). If False, each cell in the
:py:`_ConductanceSurface` will have an approximated circular distribution 
that is a single
VonMises distribution with its mode being the direction of the maximum-valued
cell in the 8-cell neighborhood and its concentration determined by
**vm_distr_kappa**.


------------------------------------------------------------------------------

**vm_distr_kappa**

.. code-block:: python

                         #concentration of distrs
                         'vm_distr_kappa':       12,

{:py:`int`, :py:`float`}

default: 12 

reset? N

This sets the concentration of the VonMises distributions used to build
the approximated circular distributions in the :py:`_ConductanceSurface`.
The default value was chosen heuristically as one that provides a reasonable
concentration in the direction of a unimodal VonMises distribution's mode 
without causing VonMises mixture distributions built from an 
evenly weighted sum of distributions pointing toward the 
8-cell-neighborhood directions to have 8 pronounced modes. 
There will probably be little need to change the default value, but if
interested then the user could create :py:`Model`\s with various values
of this parameter and then use the :py:`Model.plot_movement_surface`
method to explore the influence of the parameter on the resulting
:py:`_ConductanceSurface`\s.


------------------------------------------------------------------------------

**approx_len**

.. code-block:: python

                         #length of approximation vectors for distrs
                         'approx_len':       5000,

{:py:`int`}

default: 5000 

reset? P

This determines the length of the vector of values used to approximate each
distribution on the :py:`_ConductanceSurface` (i.e. the size of the z-axis
of the :py:`np.ndarray` used to hold all the distribution-approximations, where
the y and x axes have the same dimensions as the :py:`Landscape`). The default
value of 5000 is fine for many cases, but may need to be
reduced depending on the :py:`Landscape` dimensions (because for a larger
:py:`Landscape`, say 1000x1000 cells, it would create a 
:py:`_ConductanceSurface` that is roughly 4Gb,
and if the :py:`Layer` on which the :py:`_ConductanceSurface` is based will be
undergoing landscape changes then numerous versions of an object of this size
would need to be generated when the :py:`Model` is built and held in memory).
The value to use for this parameter will depend on the size of the
:py:`Landscape`, the exact scenario being simulated, and the memory of the
machine on which the :py:`Model` is to be run.

                   

^^^^^^^^^^^^^^^^^^^^
_GenomicArchitecture
^^^^^^^^^^^^^^^^^^^^

------------------------------------------------------------------------------

**gen_arch_file**

.. code-block:: python

                  'gen_arch': {
                      #file defining custom genomic arch
                      'gen_arch_file':            None,

{:py:`str`, :py:`None`}

default: {:py:`None`, :py:`'<your_model_name>_spp-<n>_gen_arch.csv'`

reset? P

This arguments indicates whether a custom genomic architecture file should
be used to create a :py:`Species`' :py:`GenomicArchitecture`, and if so,
where that file is located. If the value is :py:`None`, no file will be
used and the values of this :py:`Species`' other genomic
architecture parameters in the parameters file will be used to create
the :py:`GenomicArchitecture`. If the value is a :py:`str` pointing to a
custom genomic-architecture file 
(i.e. a CSV file with loci as rows and 'locus_num',
'p', 'dom', 'r', 'trait', and 'alpha' as columns stipulating the starting
allele frequencies, dominance values, inter-locus recombination rates,
trait names, and effect sizes of all loci; values can be left blank if not applicable).
Geonomics will create an empty
file of this format for each :py:`Species` for which the 
'genomes' argument is given the value 'custom' when
:py:`gnx.make_parameters_file` is called (which will be saved as
'<your_model_name>_spp-<n>_gen_arch.csv'). 

Note that when Geonomics reads in a custom genomic architecture file
to create a :py:`Model`, it will check
that the length (i.e. number of rows) in this file is equal to the length
stipulated by the **L** parameter, and will also check that the first value
at the top of the 'r' column is 0.5 (which is used to implement independent
assortment during gametogenesis). If either of these checks fails,
Geonomics throws an Error.


------------------------------------------------------------------------------

**L**

.. code-block:: python
 
                      #num of loci
                      'L':                        1000,

:py:`int`

default: 1000

reset? P

This defines the total number of loci in the genomes in a
:py:`Species`.


------------------------------------------------------------------------------

**l_c**

.. code-block:: python
                        
                      #num of chromosomes
                      'l_c':                      [100],

:py:`list` of :py:`int`\s

default: :py:`[100]`

reset? P

This defines the lengths (in number of loci) of each of the chromosomes 
in the genomes in a :py:`Species`.  Note that the sum of this :py:`list`
must equal **L**, otherwise Geonomics will throw an Error. 
Also note that Geonomics models genomes as single **L** x 2
arrays, where separate chromosomes are delineated by points along
the genome where the recombination rate is 0.5;
thus, for a model where recombination rates are often at or near 0.5, this
parameter will have little meaning.


------------------------------------------------------------------------------

**start_p_fixed**

.. code-block:: python
                        
                    #whether starting allele frequencies should be fixed at 0.5
                    'start_p_fixed':                      True,

:py:`bool`

default: True

reset? P

This indicates whether the starting 1-allele frequencies at all loci
should be set fixed at 0.5. Defaults to True.


------------------------------------------------------------------------------

**mu_neut**

.. code-block:: python

                      #genome-wide per-base neutral mut rate (0 to disable)
                      'mu_neut':                  0,

:py:`float`

default: 1e-9

reset? P

This defines the genome-wide per-base neutral mutation rate.
This value can be set to 0 to disable neutral mutation.


------------------------------------------------------------------------------

**mu_delet**

.. code-block:: python

                      #genome-wide per-base deleterious mut rate (0 to disable)
                      'mu_delet':                 0,

:py:`float`

default: 0

reset? P

This defines the genome-wide per-base deleterious mutation rate.
This value can be set to 0 to disable deleterious mutation. Note that all
deleterious mutation will fall outside the loci that affect any :py:`Trait`\s
a :py:`Species` may have, and will behave simply as globally
deleterious mutations (i.e. mutations that reduce the mutated
:py:`Individual`'s fitness regardless of that :py:`Individual`'s
spatial location).


------------------------------------------------------------------------------

**delet_alpha_distr_shape**

.. code-block:: python

                      #shape of distr of deleterious effect sizes
                      'delet_alpha_distr_shape':      0.2,

:py:`float`

default: 0.2

reset? P

This defines the shape parameter of the gamma distribution from which
the effect sizes of deleterious loci are drawn. (Values drawn will be
truncated to the interval [0,1].)


------------------------------------------------------------------------------

**delet_alpha_distr_scale**

.. code-block:: python

                      #scale of distr of deleterious effect sizes
                      'delet_alpha_distr_scale':      0.2,

:py:`float`

default: 0.2

reset? P

This defines the scale parameter of the gamma distribution from which
the effect sizes of deleterious loci are drawn. (Values drawn will be
truncated to the interval [0,1].)


------------------------------------------------------------------------------

**r_distr_alpha**

.. code-block:: python

                      #alpha of distr of recomb rates
                      'r_distr_alpha':            None,

{:py:`float`, :py:`None`}

default: None

reset? P

This defines the alpha parameter of the beta distribution from which
interlocus recombination rates are drawn. (Values drawn will be truncated to
the interval [0, 0.5].) Defaults to None, which will coerce all recombination
rates to 0.5 (i.e. will make all loci independent).


------------------------------------------------------------------------------

**r_distr_beta**

.. code-block:: python

                      #beta of distr of recomb rates
                      'r_distr_beta':            None,

{:py:`float`, :py:`None`}

default: None,

reset? P

This defines the beta parameter of the beta distribution from which
interlocus recombination rates are drawn. (Values drawn will be truncated to
the interval [0, 0.5].) Defaults to None, which will coerce all recombination
rates to 0.5 (i.e. will make all loci independent).


------------------------------------------------------------------------------

**dom**

.. code-block:: python

                      #whether loci should be dominant (for allele '1')
                      'dom':                      False,

:py:`bool`

default: False

reset? P

This indicates whether loci should be treated as dominant (if True) 
for the '1' allele  or as codominant (if False). Codominance is the default
behavior, because it is assumed that Geonomics will often be used
to model quantitative traits, for which this is a reasonable assumption.


------------------------------------------------------------------------------

**pleiotropy**

.. code-block:: python

                      #whether to allow pleiotropy
                      'pleiotropy':               False,

:py:`bool`

default: False

reset? P

This indicates whether pleiotropy should be allowed. If True, loci will be
permitted to contribute to more than one :py:`Trait`.


------------------------------------------------------------------------------

**recomb_rate_custom_fn**

.. code-block:: python

                      #custom fn for drawing recomb rates
                      'recomb_rate_custom_fn':    None,

{:py:`function`, :py:`None`}

default: :py:`None`

reset? P

This parameter allows the user to provide a custom function according to which
interlocus recombination rates will be assigned. If set to :py:`None`, the
default behavior (i.e. recombination rates chosen from a beta distribution
using **r_distr_alpha** and **r_distr_beta**) will be used.


------------------------------------------------------------------------------

**n_recomb_paths_mem**

.. code-block:: python

                      #number of recomb paths to hold in memory
                      'n_recomb_paths_mem': int(1e4),

:py:`int`

default: :py:`int(1e4)`

reset? P

This defines the maximum number of recombination paths for Genomics to hold in
memory at one time. Geonomics models recombination by using the interlocus
recombination rates to draw a large number of recombination 'paths'
along the Lx2 genome array (when the :py:`Model` is first built), and
then shuffling and cycling through those recombination paths as 
needed during :py:`Model` runs. Of the total number of paths created, some
subset will be held in memory (the number of these is defined by
this parameter), while the remainder will live in a temporary
file (which is occasionally read in whenever the paths in memory are close to
being used up). Thus, to avoid problems, the number provided to this parameter
should be comfortably larger than the largest anticipated number of
recombination paths that will be needed during a single mating event (i.e.
larger than two times the largest antipicated number of offspring to be born
to the :py:`Species` during one timestep).


------------------------------------------------------------------------------

**n_recomb_paths_tot**

.. code-block:: python

                      #total number of recomb paths to simulate
                      'n_recomb_paths':           int(1e5),

This defines the total number of recombination paths that Geonomics will
generate. Geonomics models recombination by using the interlocus
recombination rates to draw a large number of recombination 'paths'
along the Lx2 genome array (when the :py:`Model` is first built), and
then shuffling and cycling through those recombination paths as 
needed during :py:`Model` runs. The larger the total number of these paths
that is created, the more closely Geonomics will model truly
free recombination and the more prceisely it will model the exact
interlocus recombination rates defined in a :py:`Species`'
:py:`GenomicArchitecture`.


------------------------------------------------------------------------------

**allow_ad_hoc_recomb**

.. code-block:: python

                      #whether to generate recombination paths at each timestep
                      'allow_ad_hoc_recomb': False,

:py:`bool`

default: False

reset? P

This determines whether or not recombinants should be drawn each timestep
(rather than recombination paths being drawn and stored when a model is first
built, then used randomly throught the model run).
This is advantageous because it models recombination exactly (rather than
approximating recombination by drawing some number of fixed recombination paths
that get repeatedly used), and for combinations of larger genome sizes (L) and
larger mean population sizes (N) it avoids the memory used by storing so many
recombination paths drawn at model creation, thus making these
parameterizations feasible on computers with memory limitations).
It is disadvantageous, however, because it runs somewhat slower than the
default approach (recombinants drawn at model creation) for a range of
L and N values, and also because it is only available for parameterizations
with homogeneous recombination across the genome.


------------------------------------------------------------------------------

**mut_log**

.. code-block:: python

                      #whether to save mutation logs
                      'mut_log':                  None,


{:py:`str`, :py:`None`}

default: :py:`None`

reset? P

This indicates the location of the mutation-log file where  Geonomics should
save a record of each mutation that occurs for a :py:`Species`
:py:`Species`, for each iteration. If :py:`None`, no mutation log
will be created and written to.



""""""
Traits
""""""

------------------------------------------------------------------------------

**trait_<n>**

.. code-block:: python
 
              #trait name (TRAIT NAMES MUST BE UNIQUE!) 
              'trait_0' :   {

{:py:`str`, :py:`int`}

default: :py:`trait_<n>` 

reset? P

This parameter defines the name for each :py:`Trait`.
(Note that unlike most parameters, this parameter is a :py:`dict` key, 
the value for which is a :py:`dict`
of parameters defining the :py:`Trait` being named.) As the capitalized
reminder in the parameters states, each :py:`Trait`
must have a unique name (so that a parameterized 
:py:`Trait` isn't overwritten in the :py:`ParametersDict` by a
second, identically-named :py:`Trait`; Geonomics
checks for unique names and throws an Error if this condition is not met.
:py:`Trait` names can, but needn't be, descriptive of what each 
:py:`Trait` represents. Example valid values include: 0, 'trait0',
'tmp_trait', 'bill length'. Names default to :py:`trait_<n>`,
where n is a series of integers starting from 0 and counting the
number of :py:`Trait`\s for this :py:`Species`.


------------------------------------------------------------------------------

**layer**

.. code-block:: python

                              #trait-selection Layer name
                              'layer':                'layer_0',

:py:`str`

default: :py:`'layer_0'`

reset? P

This indicates, by name, the :py:`Layer` that serves as the selective force
acting on this :py:`Trait`. (For example, if this Trait is selected upon by
annual mean temperature, then the name of the :py:`Layer` 
representing annual mean temperature should be provided here.)


------------------------------------------------------------------------------

**phi**

.. code-block:: python

                              #polygenic selection coefficient
                              'phi':                  0.05,

{:py:`float`, :py:`np.ndarray` of :py:`float`\s}

default: 0.05

reset? P

This defines the polygenic selection coefficient on this :py:`Trait` (i.e
the selection coefficient acting on the phenotypes, rather than the genotypes,
of this :py:`Trait`). The effect of this value can be thought of as the
reduction (from 1) in an :py:`Individual`'s survival probability when that
:py:`Individual` is maximally unfit (i.e. when that :py:`Individual` has a
phenotypic value of 1.0 but is located in a location with an environmental
value of 0.0, or vice versa). When the value is a :py:`float` then the
strength of selection will be the same for all locations on the
:py:`Landscape`. When the value is an :py:`np.ndarray` of
equal dimensions to the :py:`Landscape` then the strength of
selection will vary across space, as indicated by the values in this array
(what Geonomics refers to as a "spatially contingent" selection regime).


------------------------------------------------------------------------------

**n_loci**

.. code-block:: python

                              #number of loci underlying trait
                              'n_loci':               1,

:py:`int`

default: 10

reset? P

This defines the number of loci that should contribute to the phenotypes
of this :py:`Trait`. These loci will be randomly drawn from across the
genome.


------------------------------------------------------------------------------

**mu**

.. code-block:: python

                              #mutation rate at loci underlying trait
                              'mu':                   1e-9,

:py:`float`

default: 1e-9

reset? P

This defines the mutation rate for this :py:`Trait` (i.e. the rate at which
mutations that affect the phenotypes of this :py:`Trait` will arise). Set to
0 to disable mutation for this :py:`Trait`.


------------------------------------------------------------------------------

**alpha_distr_mu**

.. code-block:: python

                              #mean of distr of effect sizes
                              'alpha_distr_mu' :      0.1,
:py:`float`

default: 0.1

reset? N

This defines the mean of the normal distribution from which a :py:`Trait`'s
initially parameterized loci and new mutations' effect sizes are drawn (with
the exception of monogenic traits, whose starting locus always has an alpha
value of 0.5, but whose later mutations are influenced by this parameter).
For effect sizes drawn from a distribution, it is recommended
to set this value set to 0 and adjust **alpha_distr_sigma**.
For fixed effect sizes, set this value to the fixed
effect size and set **alpha_distr_sigma** to 0; effects will alternate
between positive and negative when they are assigned to loci.
In either case, new mutations in a :py:`Trait`
will then be equally likely to decrease or increase :py:`Individual`\s'
phenotypes from the multigenic baseline phenotype of 0.5 (which is also
the central value on a Geonomics :py:`Landscape`).
It is also recmmended that the user consider the number of loci for a trait
when setting the fixed or distributed effect sizes; for example, for a trait
with 10 underlying loci, an average or fixed absolute effect size of 0.1
will enable phenotypes that cover the range of values on a
Geonomics :py:`Landscape` (i.e. phenotypes 0 <= z <= 1), whereas
0.05 will likely not enable that full range of phenotypes, and 0.5 will
generate many phenotypes that fall outside that range and will be selected
against at all locations on the :py:`Landscape`.


------------------------------------------------------------------------------

**alpha_distr_sigma**

.. code-block:: python

                              #variance of distr of effect size
                              'alpha_distr_sigma':    0,

:py:`float`

default: 0

reset? P

This defines the standard deviation of the normal distribution from which
a :py:`Trait`'s new mutations' effect sizes are drawn. 
For effect sizes drawn from a distribution, it is recommended
to set this value set to some nonzero number
and set **alpha_distr_mu** to 0. For fixed effect sizes,
set this value to 0 and set **alpha_distr_mu** to the fixed effect size;
effects will alternate between positive and negative when they are
assigned to loci. In either case, new mutations in a :py:`Trait`
will then be equally likely to decrease or increase :py:`Individual`\s'
phenotypes from the multigenic baseline phenotype of 0.5 (which is also
the central value on a Geonomics :py:`Landscape`).


------------------------------------------------------------------------------

**max_alpha_mag**

.. code-block:: python

                              #max allowed magnitude for an alpha value
                              'max_alpha':            None,

{:py:`float`}

default: None

reset? P

This defines the maximum value that can be drawn for a locus' effect size
(i.e. alpha). Defaults to None, but the user may want to set this to some
reasonable value, to prevent chance creation of loci with extreme effects.


------------------------------------------------------------------------------

**gamma**

.. code-block:: python

                              #curvature of fitness function
                              'gamma':                1,

{:py:`int`, :py:`float`}

default: 1

reset? N

This defines the curvature of the fitness function (i.e.
how fitness decreases as the absolute difference between an 
:py:`Individual`'s optimal and actual phenotypes increases). The user
will probably have no need to change this from the default value of 1
(which causes fitness to decrease linearly around the optimal
phenotypic value). Values < 1 will cause the fitness function to be
concave up; values > 1 will cause it to be concave down.


------------------------------------------------------------------------------

**univ_adv**

.. code-block:: python

                              #whether the trait is universally advantageous
                              'univ_adv':             False

:py:`bool`

default: False

reset? P

This indicates whether whether selection on a :py:`Trait` should be
universal (i.e. whether a phenotype of 1 should be optimal everywhere
on the :py:`Landscape`). When set to True, selection of the :py:`Trait`
will be directional on the entire :py:`Species`, regardless 
of :py:`Individual`\s' spatial contexts. 


^^^^^^^^^^^^^^
Species change
^^^^^^^^^^^^^^


""""""""""""""""""
Demographic change
""""""""""""""""""

------------------------------------------------------------------------------

**kind**

.. code-block:: python

                          #kind of event {'monotonic', 'stochastic',
                          #'cyclical', 'custom'}
                          'kind':             'monotonic',

{:py:`'monotonic'`, :py:`'stochastic'`, :py:`'cyclical'`, :py:`'custom'`}

default: :py:`'monotonic'`

reset? P

This indicates what type of demographic change is being parameterized.
Each event has a certain length (in timesteps; defined by the **start** and
**end** parameters). Note that of the other parameters in this section, only
those that are necessary to parameterize the type of change event indicated
here will be used.

In :py:`'monotonic'` change events, a :py:`Species`' 
carrying capacity raster (K) is multiplied by a constant factor 
(**rate**) at each timestep during the event. 
In :py:`'stochastic'` change events, K fluctuates
around the baseline value (i.e. the K-raster at the time that the change event
begins) at each required timestep during the event (where the sizes of the
fluctuations are drawn from the distribution indicated by
**distr**, the floor and ceiling on those sizes are set by
**size_range**, and the required timesteps are determined by **interval**). 
In :py:`'cyclical'` change events, K undergoes a number (indicated
by **n_cycles**) of sinusoidal cycles between some minimum and maximum
values (indicated by **size_range**). 
In :py:`'custom'` change events, the baseline K is multiplied by a series
of particular factors (defined by **sizes**) at a series of particular
timesteps (defined by **timesteps**).


------------------------------------------------------------------------------

**start_t**

.. code-block:: python

                          #starting timestep
                          'start_t':            50,

:py:`int`

default: 50

reset? P

This indicates the timestep at which the demographic change event
should start.


------------------------------------------------------------------------------

**end_t**

.. code-block:: python

                          #ending timestep
                          'end_t':            100,

:py:`int`

default: 100

reset? P

This indicates the last timestep of the change event.


------------------------------------------------------------------------------

**rate**

.. code-block:: python

                          #rate, for monotonic change
                          'rate':             1.02,

:py:`float`

default: 1.02

reset? P

This indicates the rate at which a :py:`'monotonic'` change event should occur.
At each timestep during the event, a new carrying capacity raster (K)
will be calculated by multiplying the previous step's K by this factor.
Thus, values should be expressed relative to 1.0 indicating no change.


------------------------------------------------------------------------------

**interval**

.. code-block:: python

                          #interval of changes, for stochastic change
                          'interval':         1,

:py:`int`

default: 1

reset? P

This indicates the interval at which fluctutations should occur during a
:py:`'stochastic'` change event (i.e. the number of timesteps to wait
between fluctuations).


------------------------------------------------------------------------------

**distr**

.. code-block:: python

                          #distr, for stochastic change {'uniform', 'normal'}
                          'distr':            'uniform',

{:py:`'uniform'`, :py:`'normal'`}

default: :py:`'uniform'`

reset? P

This indicates the distribution from which to draw the sizes of
fluctuations in a :py:`'stochastic'` change event. Valid options are
`'uniform'` and `'normal'`.


------------------------------------------------------------------------------

**n_cycles**

.. code-block:: python

                          #num cycles, for cyclical change
                          'n_cycles':         10,

:py:`int`

default: 10

reset? P

This indicates the number of cyclical fluctuations that should occur during
a :py:`'cyclical'` change event.


------------------------------------------------------------------------------

**size_range**

.. code-block:: python

                          #min & max sizes, for stochastic & cyclical change
                          'size_range':       (0.5, 1.5),

:py:`tuple` of :py:`float`\s

default: :py:`(0.5, 1.5)`

reset? P

This defines the minimum and maximum sizes of fluctuations that can occur
during :py:`'stochastic'` and :py:`'cyclical'` change events.


------------------------------------------------------------------------------

**timesteps**

.. code-block:: python

                          #list of timesteps, for custom change
                          'timesteps':        [50, 90, 95],

:py:`list` of :py:`int`\s

default: [50, 90, 95]

reset? P

This defines the series of particular timesteps at which fluctutations should
occur during a :py:`'custom'` change event.


------------------------------------------------------------------------------

**sizes**

.. code-block:: python

                          #list of sizes, for custom change
                          'sizes':        [2, 5, 0.5],

:py:`list` of :py:`float`\s

default: [2, 5, 0.5]

reset? P

This defines the series of particular fluctutations that should occur 
during a :py:`'custom'` change event.



"""""""""""""""""""
Life-history change
"""""""""""""""""""

------------------------------------------------------------------------------

**<life_hist_param>**

.. code-block:: python

                          #life-history parameter to change
                          '<life_hist_param>': {


:py:`str`

default: :py:`'<life_hist_param>'`

reset? P

This indicates the life-history parameter to be changed by this life-history
change event. (Note that unlike most parameters, this parameter is 
a :py:`dict` key, the value for which is a :py:`dict`
of parameters controlling how the life-history parameter that is named
will change.)


------------------------------------------------------------------------------

**timesteps**

.. code-block:: python

                          #list of timesteps
                          'timesteps':        [],

:py:`list` of :py:`int`\s

default: :py:`[]`

reset? P

This indicates the timesteps at which the life-history parameter being changed
should change (to the values indicated by **vals**).


------------------------------------------------------------------------------

**vals**

.. code-block:: python

                          #list of values
                          'vals':        [],

:py:`list` of :py:`float`\s

default: :py:`[]`

reset? P

This indicates the values to which the life-history parameter being changed
should change (at the timesteps indicated by **timesteps**).
