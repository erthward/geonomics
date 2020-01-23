.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


***************
Data structures
***************

The following sections discuss the structure and function of the key
Geonomics classes. Users will interface with these classes more or less
directly when running Geonomics models, so a fundamental understanding of how 
they're organized and how they work will be useful.

As you're reading, this conceptual diagram should be a useful reference.
Most of the core data structures are depicted in the central image.

.. image:: ../conceptual_diagram.png



=======================================
:py:`Landscape` and :py:`Layer` objects
=======================================

One of the core components of a Geonomics model is the land. The land is
modeled by the :py:`Landscape` class. This class is an 
integer-keyed :py:`dict` composed of numerous instances of the
class :py:`Layer`. Each :py:`Layer` represents a separate 
environmental variable (or 'layer', in GIS terminology),
which is modeled a 2d Numpy array (or raster; in
attribute 'rast'), of identical dimensions to each 
other :py:`Layer` in the :py:`Landscape`
object, and with the values of its environmental variable 'e' constrained to
the interval [0 <= e <= 1]. Each :py:`Layer` can be initialized from its own
parameters subsection within the 'land' parameters section of a Geonomics
parameters file. 

For each :py:`Species` (see section ':py:`Individuals`
and :py:`Species`', below), the different :py:`Layer`
layers in the :py:`Landscape` can be used to model habitat 
viability, habitat connectivity, or variables imposing spatially varying
natural selection. :py:`Landscape` and :py:`Layer` objects
also contain some metatdata (as public attributes), including
the resolution (attribute 'res'), upper-left corner ('ulc'),
and projection ('prj'), which default to 1, (0,0), and None but
will be set otherwise if some or all of the :py:`Layer` layers are read in from
real-world GIS rasters.


-------------------------------------------------------------------------------

===========================================================
Genomes, :py:`GenomicArchitecture`, and :py:`Trait` objects
===========================================================

:py:`Individual` objects (see section ':py:`Individuals`
and :py:`Species`', below) can optionally be assigned genomes.
If they are, each :py:`Individual`'s genome is modeled as a 
2-by-L Numpy array (where 2 is the ploidy, currently fixed at
diploidy, and L is genome length) containing 0s and 1s (because
Geonomics strictly models biallelic SNPs, i.e SNPs with '0'- and '1'-alleles). 

The parameter L, as well as numerous other genomic parameters (including 
locus-wise starting frequencies of the 1 alleles; locus-wise dominance effects;
locus-wise recombination rates; and genome-wide mutation rates for neutral, 
globally deleterious, and adaptive loci), are controlled by the 
:py:`GenomicArchitecture` object pertaining to the :py:`Species` to which an 
:py:`Individual` belongs. (For the full and detailed list of attributes in a 
:py:`GenomicArchitecture` object, see its class documentation, below.)
The genomes of the initial :py:`Individual`\s 
in a simulation, as well as those of 
:py:`Individual`\s in subsequent generations, are either drawn
or recombined, and are mutated, according to the values stipulated 
by the :py:`GenomicArchitecture` of
their :py:`Species`. The user can create a species with a 
:py:`GenomicArchitecture` and with corresponding
genomes by including a 'genome' subsection in that
species' section of the Geonomics parameters file (and 
setting the section's various parameters to their desired values). 

Geonomics can model :py:`Individual`\s' phenotypes.
It does this by allowing the 
user to create an arbitrary number of distinct :py:`Trait`\s
for each :py:`Species`. Each trait is
represented by a :py:`Trait` object, which 
maps genomic loci onto that trait, maps effect sizes ('alpha') onto those loci,
and sets the trait's polygenic selection
coefficient ('phi'). An :py:`Individual`'s
phenotype for a given trait is calculated as the 'null phenotype' plus a 
weighted sum of the products of its 'effective genotypes' at all loci 
underlying that :py:`Trait` and the effect sizes (i.e. 'alpha') of those loci:

.. math::

   z_{i,t} = null\_genotype + \sum_{l = 0}^{n} \alpha_{t,l} g_{i,l}

where :math:`z_{i,t}` is the phenotype of :py:`Individual` i for trait t, 
:math:`g_{i, l}` is the genotype of the :py:`Individual` at that locus, and 
:math:`\alpha_{t,l}` is the effect size of that locus for that trait.

The 'null phenotype' refers determines what would be the phenotypic value that
an :py:`Individual` who is homozygyous for
the 0 allele at all loci for a trait.
For monogenic traits the null phenotype is 0 and the effect size is fixed at 
0.5 (such that individuals can have phenotypes of 0, 0.5, or 1); 
for polygenic traits the null phenotype is 0.5 and effect sizes can be fixed 
at or distributed around a mean value (which is controlled in the 
parameters file).

The 'effective genotype' refers to how the genotype is calculated based on the 
dominance at a locus, as indicated by the following table of genotypes:

+--------------------+------------------+------------------+
| Biallelic genotype |   Codominant     |     Dominant     |
+====================+==================+==================+
|      0 : 0         |        0         |        0         |
+--------------------+------------------+------------------+
|      0 : 1         |       0.5        |        1         |
+--------------------+------------------+------------------+
|      1 : 1         |        1         |        1         |
+--------------------+------------------+------------------+

(For the full and detailed list of attributes in a :py:`Trait` object, 
see its class documentation, below.)

Note that for maximal control over the :py:`GenomicArchitecture`
of a :py:`Species`, the user can set the value of the 'gen_arch_file' 
parameter in the parameters file to the name of a separate CSV file 
stipulating the locus numbers, starting 1-allele frequencies, dominance 
effects, traits, and inter-locus recombination rates (as columns) of 
all loci (rows) in the :py:`GenomicArchitecture`;
these values will override any other values provided in the 'genome' 
subsection of the species' parameters.


-------------------------------------------------------------------------------

============================================================
:py:`Individual`, :py:`Species`, and :py:`Community` objects
============================================================

Being that Geonomics is an individual-based model, individuals serve as 
the fundamental units (or agents) of all simulations. They are represented by
objects of the :py:`Individual` class.
Each :py:`Individual` has an index (saved 
as attribute 'idx'), a sex (attribute 'sex'), an age (attribute 'age'), 
an x,y position (in continuous space; attributes 'x' and 'y'), and a 
:py:`list` of environment values (attribute 'e'), extracted from the 
:py:`Individual`'s current cell on each :py:`Layer`
of the :py:`Landscape` on which the :py:`Individual` lives.

The :py:`Species` class is an :py:`OrderedDict`
(defined by the :py:`collections` 
package) containing all :py:`Individaul`\s, (with 
their 'idx' attributes as keys). If a :py:`Species`
has a :py:`GenomicArchitecture` then the :py:`Individual`\s
in the :py:`Species` will also each have genomes (attribute 'g'),
and the :py:`GenomicArchitecture` includes :py:`Trait`\s
then each individual will also have a :py:`list` of 
phenotype values (one per :py:`Trait`; attribute 'z') and a 
single fitness value (attribute 'fit'). (These attributes all otherwise 
default to :py:`None`.)

Each :py:`Species` also has a number of other attributes of interest. Some 
of these are universal (i.e. they are created regardless of the 
parameterization of the :py:`Model` to which a :py:`Species` inheres). These 
include: the :py:`Species`' name (attribute 'name'); its current density 
raster (a Numpy array attribute called 'N'); and the number of births,
number of deaths, and terminal population size (i.e. total number of
individuals in the :py:`Species`) of each timestep (which are 
:py:`list` attributes called 'n_births', 'n_deaths', and 'Nt'). If the 
:py:`Species` was parameterized with a
:py:`GenomicArchitecture` then that will 
be created as the 'gen_arch' attribute (otherwise this attribute will be 
:py:`None`).

All of the :py:`Species` in a :py:`Model`
are collected in the :py:`Model`'s 
:py:`Community` object. The :py:`Community` class
is simply an integer-keyed :py:`dict` 
of :py:`Species`. For the time being, the :py:`Community` object allows a 
Geonomics :py:`Model` to simulate multiple :py:`Species` simultaneously on 
the same :py:`Landscape`, but otherwise affords no additional functionality
of interest. However, its implementation will facilitate the potential 
future development of methods for interaction between :py:`Species`. 
(e.g. to simulate coevolutionary, speciation, or hybridization scenarios).


-------------------------------------------------------------------------------

===================
:py:`Model` Objects
===================

Objects of the :py:`Model` class serve as the main interface between the user 
and the Geonomics program. (While it is certainly possible for a user 
to work directly with the :py:`Landscape`
and :py:`Species` or :py:`Community` objects to 
script their own custom models, the typical user should find that the 
:py:`Model` object allows them accomplish their goals with minimal toil.)
The main affordance of a :py:`Model` object is the :py:`Model.run` method, 
which, as one could guess, will run the :py:`Model`. The typical workflow 
for creating and running a  :py:`Model` object is as follows:

  1. Create a template paramters file containing the desired sections, 
     by calling :py:`gnx.make_parameters_file` with all revelant arguments;
  2. Define the scenario to be simulated, by opening and editing that 
     parameters file (and optionally, creating/editing corresponding 
     files, e.g. genomic-architecture CSV files;
     or raster or numpy-array files to be used as :py:`Layer`\s);
  3. Instantiate a :py:`Model` object from that parameters file, by calling 
     :py:`mod = gnx.make_model('/path/to/params_filename.py')`;
  4. Run the :py:`Model`, by calling :py:`mod.run()`.

For detailed information on usage of these functions, see their docstrings.
When a :py:`Model` is run, it will:

  1. Run the burn-in (until the mininmal burn-in length stipulated in the 
     parameters file and the built-in stationarity statistics 
     determine that the burn-in is complete);
  2. Run the main model for the stipulated number of timesteps;
  3. Repeat this for the stipulated number of iterations (retaining or 
     refreshing the first run's initial :py:`Landscape` and :py:`Species` 
     objects as stipulated).

The :py:`Model` object offers one other method, however, :py:`Model.walk`, 
which allows the user to run a model, in either 'burn' or 'main' mode, 
for an arbitrary number of timesteps within a single iteration (see its 
docstring for details). This is particularly useful for running 
Geonomics within an interactive Python session. Thus, :py:`Model.walk` is 
primarily designed for passively running numerous iterations of a :py:`Model`, 
to generate data for analysis, whereas :py:`Model.walk` is primarily designed
for the purposes of learning, teaching, or debugging the package, or 
developing, exploring, introspecting, or visaulizing particular :py:`Model`\s. 


-------------------------------------------------------------------------------

=================
Secondary classes
=================

The typical user will not need to access or interact with the following 
classes in any way. They will, however, parameterize them in the 
parameters file by either leaving or altering their default values. Geonomics 
sets generally sensible default parameter values wherever possible, 
but for some scenarios they may not be adequate, and for some parameters 
(e.g. the window-width used by the _DensityGridStack; see below), there is 
no "one-size-fits-most" option. Thus, it is important that the user
have a basic acquaintance with the purpose and operation of these classes.

-------------------------
:py:`_ConductanceSurface`
-------------------------

The :py:`_ConductanceSurface` class allows Geonomics
to model a :py:`Species`' 
realistic movement across a spatially varying landscape. It does this by 
creating an array of circular probability distributions (i.e. VonMises 
distributions), one for each cell on the :py:`Landscape`, from which 
:py:`Individual`\s choose their directions each time they move. To create the
:py:`_ConductanceSurface` for a :py:`Species`,
the user must indicate the :py:`Layer` 
that should be used to create it (i.e. the :py:`Layer` that represents 
landscape permeability for that :py:`Species`). 
The :py:`_ConductanceSurface`'s 
distributions can be **simple (i.e. unimodal)**, such that the 
maximum value of the distribution at each cell will point toward the
maximum value in the 8-cell neighborhood; this works best for permeability 
:py:`Layer`\s with shallow, monotonic gradients, because the differences 
between permeability values of neighboring cells can be minor (e.g. a 
gradient representing the directionality of a prevalent current). 
Alternatively, the distributions can be **mixture (i.e. multimodal)**
distributions, which are weighted sums of 8 unimodal distributions, one 
for each neighboring cell, where the weights are the relative cell 
permeabilities (i.e. the relative probabilities that an :py:`Individual` would 
move into each of the 8 neighboring cells); this works best for non-monotonic, 
complex permeability :py:`Layer`\s (e.g. a DEM of a mountainous region that is 
used as a permeability :py:`Layer`). 
(The :py:`Landscape` is surrounded by a margin of 0-permeability 
cells before the :py:`_ConductanceSurface` is calculated, such 
that :py:`Landscape` edges are treated 
as barriers to movement.) The class consists 
principally of a 3d Numpy array (y by x by z, where y and x (a.k.a i and j, 
or latitude and longitude) are the dimensions of the 
:py:`Landscape` and z is the length of the vector of values 
used to approximate the distributions in each cell.

-----------------------
:py:`_DensityGridStack`
-----------------------

The :py:`_DensityGridStack` class implements an algorithm for rapid estimating 
an array of the local density of a :py:`Species`. The resulting array has a 
spatial resolution equivalent to that of the :py:`Landscape`,
and is used in all density-dependent operations (i.e. for controlling 
population dynamics). The density is estimated 
using a sliding window approach, with the window-width determining the 
neighborhood size of the estimate (thus essentially behaving like a smoothing
parameter on the density raster that is estimated, with larger window widths
producing smoother, more homogeneous rasters). The window width can be
controlled by setting the 'density_grid_window_width' parameter in the 
'mortality' section of the :py:`Species` parameters, in a parameters file;
however, if the default value (:py:`None`) is left then the window width will
default to 1/20th of the width of the :py:`Landscape`. 
Note that setting the window width to a value less than ~1/20th of the
:py:`Landscape` width is likely to result
in dramatic increases in runtime, so this is generally advised against (but
may be necessary, depending on the user's interests). The following plot
show the estimated density rasters for a 1000x1000-cell :py:`Landscape` with
a population of 50,000 individuals, using various window widths:

.. image:: ../DensityGridStack_ww_100.jpg
   :align: center

And this plot shows how :py:`_DensityGridStack` creation (plot titled 'make')
and runtime ('calc')scale with window-width for that :py:`Landscape`:

.. image:: ../densitygridstack_compute_times.png
   :align: center

-------------
:py:`_KDTree`
-------------

The :py:`_KDTree` class is just a wrapper around :py:`scipy.spatial.cKDTree`. 
It provides an optimized algorithm (the kd-tree) for finding 
neighboring points within a given search radius.
This class is used for all neighbor-searching operations (e.g. mate-search).

-------------------------
:py:`_RecombinationPaths`
-------------------------

The :py:`_RecombinationPaths` class contains a large (and customizable) 
number of :py:`bitarray`\s, each of which indicates the genome-length 
diploid chromatid numbers (0 or 1) for a
recombinant gamete produced by an :py:`Individual` of a given :py:`Species` 
(henceforth referred to as 'recombination paths'). These recombination 
paths are generated using the genome-wide recombination rates specified by 
the :py:`Species`' :py:`GeonomicArchitecture`. They are generated during 
construction of the :py:`Model`, then drawn randomly as needed (i.e.
each time an :py:`Individual` produces a gamete). This provides a 
reasonable trade-off between realistic modelling of recombination and runtime.

-------------------------------------------------
:py:`_LandscapeChanger` and :py:`_SpeciesChanger`
-------------------------------------------------

These classes manage all of the landscape changes and demographic changes 
that were parameterized for the :py:`Landscape` and
:py:`Species` objects to which they inhere. 
The functions creating these changes are defined at the outset, 
then queued and called at their scheduled timesteps.

----------------------------------------------
:py:`_DataCollector` and :py:`_StatsCollector`
----------------------------------------------

These classes manage all of the data and statistics that should be collected 
and written to file for the :py:`Model` object to which they inhere 
(as determined by the parameters file used the create the :py:`Model`). 
The types of data to be collected, or statistics to be calculated, as 
well as the timesteps at which and methods by which they're 
collected/calculated and determined at the outset, then the 
appropriate functions called at the appropriate timesteps.
