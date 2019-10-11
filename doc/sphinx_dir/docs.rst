.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


.. code-block:: python


                    #   ::::::          :::    :: :::::::::::#
             #::::::    ::::   :::      ::    :: :: ::::::::::: ::#
          #:::::::::     ::            ::   :::::::::::::::::::::::::#
        #::::::::::                      :::::::::: :::::: ::::::::  ::#
      #  : ::::  ::                    ::::  : ::    :::::::: : ::  :    #
     # GGGGG  EEEEE OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS #
    # GG     EE    OO   OO  NNN  NN  OO   OO  MM   MM   II   CC     SS     #
    # GG     EE   OO     OO NN N NN OO     OO MMM MMM   II   CC     SSSSSS #
    # GG GGG EEEE OO     OO NN  NNN OO     OO MM M MM   II   CC         SS #
    # GG   G EE    OO   OO  NN   NN  OO   OO  MM   MM   II   CC        SSS #
     # GGGGG  EEEEE OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS #
      #     :::::::::               :::::::::: ::              ::  :   : #
        #::   :::::                    :::::: :::             :::::::  #
          #    :::                      :::::  ::              ::::: #
             #  ::                      ::::                      #
                   #                                        #
                      #  :: ::    :::             #


.. image:: ./yosemite_header_img.png
   :align: center


#########
Geonomics
#########


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
A Python package for easy construction of individual-based, spatially explicit, forward-time, and highly customizable landscape-genomics simulations
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


********
Provides
********

  1. :py:`Landscape`, :py:`GenomicArchitecture`, :py:`Trait`,
     :py:`Individual`, :py:`Species`, and :py:`Community` classes
  2. A :py:`Model` class that builds objects of the aforementioned classes 
     according to the scenario stipulated in a parameters file,
     then uses them to run numerous simulation iterations
  3. Classes to help the user to save data and basic statistics 
     from those iterations at all desired timesteps
  3. Classes to help the user model arbitrarily complex landscape- and 
     demographic-change scenarios
  4. Numerous visualization methods, to help the user design, run, explore, 
     present, and explain their models' behavior and results


*****
Intro
*****

Backward-time (i.e. coalescent) simulators abound.
But they are inadequate for simulation of many scenarios of 
interest, including: natural selection on traits with arbitrary genomic 
architectures; spatially variable natural selection; simulation of species or
populations distributed continuously and moving realistically across
complex landscapes; complex demographic change simultaneous with ongoing, 
often non-stationary environmental change; and coevolutionary interactions 
between multiple species or incipient species. Few existing forward-time 
simulators can model all of these phenomena, and the few that can often 
impose a high cost of entry (e.g. learning a new, non-standard programming
language in order to write one's desired simulations). Geonomics aims to fill 
this empty niche by combining ease of use with broad extensibility. 
If it succeeds at doing this, Geonomics should prove uniquely useful
for a wide range of purposes, from intro-level educational use to
high-quality theoretical, methodological, empirical, and
applied research.

Geonomics is written in Python, a full-fledged scripting language 
that is relatively easy to learn (and fun!). In Python, it can be pretty quick
for a new user to get up to speed and start doing useful work. For work with
Geonomics, this turnaround time should be even quicker. Geonomics aims to
require minimal Python knowledge (yet maintain high extensibility for
interested, expert users). Essentially, anyone should be able to build their
own, arbitrarily complex Geonomics models as long as they know how to install
the package, open a Python console, call Python functions, and edit some
default values in a pre-packaged script. 

Geonomics will allow you to:

  - create a :py:`Landscape` with any number of :py:`Layers` in it; 
  - create any number of :py:`Species` living on that
    :py:`Landscape`, each of which is composed of a bunch of 
    independent :py:`Individual`\s, and each of which will have a bunch of
    parameters dsecribing what it's like and how it lives;
  - optionally give the :py:`Individual`\s of any :py:`Species`
    genomes, which can optionally determine phenotypes for any number 
    of :py:`Trait`\s (all of this is controlled by the
    :py:`GenomicArchitecture` that you would create for
    the :py:`Species`);
  - simulate any number of timsesteps of the evolution of those
    :py:`Species` on that :py:`Landscape`, where each timestep can include
    movement, mating, mortality (by density-dependence and optionally also by
    natural selection), and demographic, life-history, or
    environmental changes

-------------------------------------------------------------------------------

==================
Object-orientation
==================

For a more technical understanding of the model, it may be helpful to 
understand the concept of **object-oriented programming**.  Here is a very
brief tutorial for the unacquainted:

Python is a very handy language for object-oriented programming, 
and this is the primary programming paradigm in which Geonomics is written. 
Essentially, object-orientation involves: 

  1. Defining certain types of data structures, or **classes** (e.g.
     :code:`Car`), and assigning them various behaviors, or **methods**
     (e.g. :code:`honk`, :code:`turn`, :code:`brake`);
  2. Using particular data values to create individual instances, or 
     **objects** belonging to those classes (e.g. :code:`my_1986_jeep`, or
     :code:`batmobile`);
  3. Instructing those **objects** to carry out their behaviors by 'calling' 
     their **methods** (e.g. :code:`my_1986_jeep.honk()` might return "Beepity
     beep!", wheras :code:`batmobile.honk()` might
     return "<Batman theme song>"). 
     
Geonomics defines a number of **classes**
(such as the :py:`Landscape`, :py:`Layer`,
:py:`Popualtion`, :py:`GenomicArchitecture`, and :py:`Trait` classes mentioned
above. The user will use the values they specfiy in a parameters file to
create **objects** belongining to these classes. Then the user will call
key **methods** that belong to these **objects**, to get them
to carry out certain behaviors (which are what constitute the simulation).

The subsequent documentation will present the **classes** definined in
Geonomics and their key **methods**. It will explain exactly what those methods
do, and how what they do fits into the overall structure and function of 
Geonomics models.


-------------------------------------------------------------------------------

===============
Getting started
===============

**For the beginner**, we recommend the following steps:
  1. Review the following three sections ('Model organization', 'Data
     structures and classes', and 'Operations'), to get a general
     undertsanding of the logic, components, and necessary and optional
     behaviors of a Geonomics model;
  2. Skim the subsequent section ('Parameters'), to understand the structure
     of a Geonomics parameters file;
  3. Pip-install Geonomics (:bash:`$ pip install geonomics`);
  4. Open Python and run :py:`import geonomics as gnx`;
  5. Use the :py:`gnx.make_parameters_file` function, to begin
     creating template parameters files that they can tweak as desired;
  6. Use the :py:`gnx.make_model` function and then the
     :py:`Model.walk` or :py:`Model.run` methods to instantiate and run
     the models they've parameterized;
  7. Use the various :py:`Model.plot` methods to visualize
     the behavior and results of their models.

**For the `impatient` beginner**, as soon as Geonomics has been
installed, you should be able to run the following code:

.. code-block:: python

     import geonomics as gnx

     gnx.run_default_model()

This will load the Geonomics package, create a default Geonomics
parameters file in your current working directory, 
then use that file to instantiate and run a :code:`Model` using the default
parameter values.


-------------------------------------------------------------------------------

=================
The documentation
=================

Finally, some brief notes about the structure and style of this documentation: 

  It is designed to be read from from the top down; information generally 
  becomes increasingly detailed as you scroll down). However, given the 
  interrelationships between all the components of a Geonomics 
  :py:`Model`, there are inevitably places where you'll run
  into material that relates to material from another section.
  To the extent possible, we attempt to cross-reference rather than duplicate
  information.

  We assume, throughout, that Genomics has been imported :py:`as gnx` and
  that Numpy has been imported :py:`as np`.


Merry modeling!


-------------------------------------------------------------------------------

-------------------------------------------------------------------------------

******************************
Model organization and worflow
******************************

.. image:: ./flow_diagram.png

-------------------------------------------------------------------------------

***************************
Data structures and classes
***************************

The following sections discuss the structure and function of the key
Geonomics classes. Users will interface with these classes more or less
directly when running Geonomics models, so a fundamental understanding of how 
they're organized and how they work will be useful.

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

.. image:: ./DensityGridStack_ww_100.jpg
   :align: center

And this plot shows how :py:`_DensityGridStack` creation (plot titled 'make')
and runtime ('calc')scale with window-width for that :py:`Landscape`:

.. image:: ./DensityGridStack_compute_times.pdf
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


-------------------------------------------------------------------------------

**********
Operations
**********

The following sections discuss the mechanics of core Geonomics operations. 
The material here is inevitably intertwined with some of the material in 
the "Data structures and classes" section. To the extent possible, we 
attempt to cross-reference rather than duplicate information (with 
the exception of this sentence).

======================
Movement and Dispersal
======================

Movement is optional, such that turning off movement will allow the user 
to simulate sessile organisms (which will reproduce and disperse, 
but not move after dispersal; this distinction is of course irrelevant 
for a :py:`Species` with a maximum age of 1). For :py:`Species` 
with movement, :py:`Individual`\s can
move by two distinct mechanisms. **Spatially random movement**
is the default behavior; in this case, :py:`Individual`\s 
move to next locations that are determined by a random distance drawn 
from a Wald distribution and a random direction drawn from a uniform 
circular (i.e. Von Mises) distribution.  As with most distributions used 
in Geonomics, the parameters of these distributions have sensible 
default values but can be customized in a :py:`Model`'s parameters file 
(see section 'Parameters', below). 

The alternative movement mechanism that is available is 
**movement across a permeability surface**,
using a :py:`_ConductanceSurface` object.
To parameterize a :py:`_MovemementSurface` for a :py:`Species`, the user 
must create a template parameters file that includes the 
necessary parameters section for the :py:`Species` (i.e. 
the user must set 'movement' to :py:`True` and 'movement_surface' to :py:`True` 
in the :py:`Species`' arguments to the :py:`gnx.make_parameters_file` 
function (see the docstring for that function for details and an example). 
:py:`Individual`\s move to next locations determined by a random distance drawn 
from a Wald distribution and a random direction drawn from the distribution 
at the  :py:`_ConductanceSurface` cell in which which the :py:`Individual`\s 
are currently located. For details about :py:`_ConductanceSurface` creation,
see section ':py:`_ConductanceSurface`' above, or the class' docstring.

Dispersal is currently implemeneted identically to spatially random movement 
(with the caveat that the an offspring's new location is determined 
relative its parents' midpoint). But the option to use a 
:py:`_ConductanceSurface` for dispersal will be offered soon.


-------------------------------------------------------------------------------

============
Reproduction
============

Each timestep, for each :py:`Species`, all pairs of individuals within 
a certain distance of each other (i.e. the mating radius, 
which is set in the parameters file) are identified.
These pairs are subsetted if necessary (i.e. if the :py:`Species` 
requires that :py:`Individual`\s be above a certain reproductive age, 
or that they be of opposite sexes, in order to mate; these values 
can also be changed from their defaults in the parameters file). 
Remaining pairs mate probabilistically (according to a Bernoulli 
random draw with probability equal to the :py:`Species`' birth 
rate, which is also set in the parameters file).

Pairs that are chosen to mate will produce a number of new 
offspring drawn from a Poisson distribution (with lambda set in the 
parameters file). For each offspring, sex is chosen probablistically 
(a Bernoulli random draw with probability equal to the :py:`Species`' 
sex ratio), age set to 0, and location chosen by dispersal from 
the parents' midpoint (see section 'Movement and Dispersal'). For 
:py:`Species` that have genomes, offspring genomes will be a 
fusion of two recombinant genomes from each of the two parents (where 
each recombinant is indexed out a parent's genome using a recombination 
path; see section ':py:`_RecombinationPaths`'). For :py:`Species` 
with :py:`Trait`\s in their
:py:`GenomicArchitecture`\s, offspring phenotypes are 
determined at birth. Mutations are also drawn and introduced at this 
point (see section 'Mutation for details).


-------------------------------------------------------------------------------

=========
Mortality
=========

Mortality can occur as a combination of two factors: **density dependence** 
and **natural selection**. Each :py:`Individual` has a death decision drawn 
as a Bernoulli random variable with 
:math:`P(d_{i}) = 1 - P(s_{i_{dens}})P(s_{i_{fit}})`, where :math:`P(d_{i})` 
is the probability of death of :py:`Individual` :math:`i`, and 
:math:`P(s_{i_{dens}})` and :math:`P(s_{i_{fit}})` are the probabilities of 
survival of :py:`Individual` :math:`i` given density-dependence and 
fitness. The probability of density-dependent death is contingent on an 
:py:`Individual`'s x,y location
(i.e. the cell in which they're currently located. 
And an :py:`Individual`'s probability of survival due to fitness 
is just equal to the product of their absolute fitness (:math:`\omega`) 
for each of the :py:`Individual`'s :math:`m` :py:`Trait`\s. 
Thus the equation for an :py:`Individual`'s probability of death becomes:

.. math::
   P(d_{i}) = 1 - (1 - P(d_{x,y})) \prod_{p = 1}^{m}\omega_{i,p}

The following two sections explain in detail the implementation and 
calculation of the two halves of the right side of this equation.

------------------
Density dependence
------------------

Density dependence is implemented using a spatialized form of the class 
logistic growth equation (:math:`\frac{\mathrm{d}
N_{x,y}}{\mathrm{d}t}=rN_{x,y}(1-\frac{N_{x,y}}{K_{x,y}})`, 
where the x,y subscripts refer to
values for a given cell on the :py:`Landscape`).
Each :py:`Species` has a carrying-capacity raster (a 2d Numpy array; 
attribute 'K'), which is defined in the parameters file to be 
one of the :py:`Layer`\s in the :py:`Landscape`.
The comparison between this raster and 
the population-density raster calculated at each timestep serves as the 
basis for the spatialized logistic growth equation, because both 
equations can be calculated cell-wise for the entire extent of the 
:py:`Landscape` (using the :py:`Species`'
intrinsic growth rate, the attribute 
'R', which is set in the parameters file).

The logistic equation returns an array of instantaneous population growth 
rates within each cell. We can derive from this the density-dependent 
probability of death at each cell by subtracting an array of the expected 
number of births at each cell, then dividing by the array of 
population density:

.. math::
   P(d_{x,y}) = E[N_{d;x,y}]/N_{x,y} = \frac{E[N_{b;x,y}] 
    - \frac{\mathrm{d}N_{x,y}}{\mathrm{d}t}}{N_{x,y}}

The expected number of births at each cell is calculated as a density 
raster of the number of succesful mating pairs, multiplied by the expected 
number of births per pair (i.e. the expectation of the Poisson 
distribution of the number of offspring per mating pair, which 
is just the distribution's paramater lambda). 

---------
Selection
---------

Selection on a :py:`Trait` can exhibit three regimes: **spatially divergent**, 
**universal**, and **spatially contingent**. **Spatially divergent** selection 
is the default behavior, and the most commonly used; in this form of 
selection, an :py:`Individual`'s fitness depends on the absolute difference 
between the :py:`Individual`'s phenotypic value and the environmental
value of the relevant :py:`Layer` (i.e. the :py:`Layer` that represents the 
environmental variable acting as the selective force) in the cell where 
the :py:`Individual` is located.

**Universal** selection (which can be toggled using the 'univ_adv' 
parameter with a :py:`Trait`'s section in the parameters file) occurs 
when a phenotype of 1 is optimal everywhere on the :py:`Landscape`. In other 
words, it represents directional selection on an entire :py:`Species`,
regardless of :py:`Individual`\s' spatial contexts. (Note that this can
be thought of as operating the same as spatially divergent selection,
but with the environmental variable driving natural selection being
represented by an array in which all cells are equal to 1.)

Under **spatially contingent** selection, the selection coefficient of a 
:py:`Trait` varies across space, such that the strength of selection 
is environmentally determined in some way. Importantly, this selection regime
is *not mutually exclusive* with the other two; in other words, 
selection on a certain :py:`Trait` be both spatially contingent 
and either spatially divergent or universal. Spatially contingent selection 
can be implemented by providing an array of values (equal in dimensions 
to the :py:`Landscape`) to the 'phi' value of a
:py:`Trait`, rather than a scalar 
value (which could be done within the parameters file itself, but may be 
more easily accomplished as a step between reading in a parameters file and 
instantiating a :py:`Model` object from it). (Note that non-spatailly
cotingent selection could in fact be thought of as a special case of
spatially contingent selection, but where the array of selection-coefficients
has the same value at each cell.)

All possible combinations of the three selection regimes of selection can all 
be thought of as special cases of the following equation for the fitness of 
:py:`Individual` :math:`i` for :py:`Trait` :math:`p` (:math:`\\omega_{i,p}`):

.. math::
   \omega_{i,p}= 1 - \phi_{p;x,y} (\mid e_{p;x,y} - z_{i;p} \mid)^{\gamma_{p}}

where :math:`\\phi_{p;x,y}` is the selection coefficient of trait 
:math:`p`; :math:`e_{p;x,y}` is the environmental variable of the 
relevant :py:`Layer` at :py:`Individual` :math:`i`'s x,y location
(which can also be thought of as the :py:`Individual`'s optimal 
phenotype); :math:`z_{i;p}` is :py:`Individual` :math:`i`'s (actual) 
phenotype for :py:`Trait` :math:`p`; and :math:`gamma_{p}` controls 
the curvature of the fitness function (i.e. how fitness decreases as
the absolute difference between an :py:`Individual`'s 
optimal and actual phenotypes increases; the default value of 1 causes 
fitness to decrease linearly around the optimal phenotypic value). 


-------------------------------------------------------------------------------

========
Mutation
========

Geonomics can model mutations of three different types: **neutral**, 
**deleterious**, and **trait** mutations. These terms don't map 
precisely onto the traditional population-genetic
lingo of "neutral", "deleterious", and "beneficial", but they 
are more or less analogous:

- **Neutral** mutations are the same conceptually in Geonomics as 
  they are in the field of population genetics in general: 
  They are mutations that have no effect on the fitness of
  the individuals in which they occur.
- **Deleterious** mutations in Geonomics are also conceptually the 
  same in Geonomics and in population genetics: They negatively impact 
  the fitness of the individuals in which they occur.
- **Trait** mutations are the place where the Geonomics concept and 
  the population-genetic concept diverge: In Geonomics, natural selection
  acts on the phenotype, not the genotype (although these concepts are 
  identical if a :py:`Trait` in monogenic), and it is (by default, 
  but not always; see section 'Selection', above) divergent. For this reason
  it would be a misnomer to call mutations that influence a given 
  :py:`Trait`'s phenotypes 'beneficial' -- even though that term is the closest
  population-genetic concept to this concept as it is employed in Geonomics -- 
  because the same mutant genotype in the same :py:`Individual`
  could have opposite effects on that :py:`Individual`'s fitness 
  in different environmental contexts (i.e. it could behave as
  a beneficial mutation is one region of the :py:`Landscape` 
  but as a deleterious mutation in another). 


-------------------------------------------------------------------------------

====================
Species interactions
====================

This functionality is not yet included available. But the Community class was 
created in advance recognition that this functionality could be desirable 
for future versions (e.g. to simulate coevolutionary, speciation, or 
hybridization scenarios).


-------------------------------------------------------------------------------

========================================
:py:`Landscape` and :py:`Species` change
========================================

For a given :py:`Layer`, any number of change events 
can be planned. 
In the parameters file, for each event, the user stipulates the initial
timestep; the final timestep; the end raster (i.e. the array 
of the :py:`Layer` that will exist after the event is complete, defined using
the **end_rast** parameter); and the 
interval at which intermediate changes will occur.  When the :py:`Model` is 
created, the stepped series of intermediate :py:`Layers` (and 
:py:`_ConductanceSurface` objects,
if the :py:`Layer` that is changing serves as the basis for a 
:py:`_ConductanceSurface` for any :py:`Species`) will be 
created and queued, so that they will swap out accordingly at the appropriate 
timesteps.

For a given :py:`Species`, any number of demographic change events can 
also be planned. In the parameters file, for each event, the user 
stipulates the type of the event ('monotonic', 'cyclical', 'random', or 
'custom') as well as the values of a number of associated 
parameters (precisely which parameters depdends on the type of event chosen).
As with :py:`Landscape` change events, all necessary stepwise changes will be 
planned and queued when the :py:`Model` is created, and will be 
executed at the appropriate timesteps.

It is also possible to schedule any number of instantaneous changes 
to some of the life-history parameters of a :py:`Species` (e.g. birth rate; 
the lambda parameter of the Poisson distribution determining the number of 
offspring of mating events). This functionality is currently minimalistic, 
but will be more facilitated in future versions.


-------------------------------------------------------------------------------

*************
Visualization
*************

Each :py:`Model` object has a variety of visualization methods 
(:py:`Model.plot`, :py:`Model.plot_fitness`, etc.),
which aim to help users design, run, explore, present,
and explain their :py:`Models`' behavior and results.
These methods can be called at any time (e.g. as 
soon as the :py:`Model` has been created, or after it has
run for any number of timesteps); but it is worth mentioning that some 
methods may be invalid depending on the point in model-time at 
which they're called (e.g. :py:`Model.plot_genotype`, 
:py:`Model.plot_phenotype`, and :py:`Model.plot_fitness`
cannot be run for :py:`Models` that have not yet been burned in,
as they will not yet have genomes assigned), or on
the :py:`Species` for which they're called 
(e.g. the aforementioned methods cannot create plots for a :py:`Species` 
that has no :py:`GenomicArchitecture`; and likewise, the 
:py:`Species.plot_demographic_changes` method cannot be called for a 
:py:`Species` for which demographic changes were not parameterized).

Below is a list of the visualization methods available, with example
output for each (generated from the default Geonomics :py:`Model`):


:: VIZ FUNCTION

:: PLOT
  
:: VIZ FUNCTION

:: PLOT

:: ...


-------------------------------------------------------------------------------

**********
Parameters
**********

In order to create and run a Geonomics :py:`Model`, you will need a valid
Geonomics parameters file. No worry though -- this is very easy to create!
To generate a new, template parameters file, you will simply call the
:py:`gnx.make_parameters_file` function, feeding it the appropriate
arguments (to indicate how many :py:`Species` and :py:`Layer`\s you
want to include in your :py:`Model`; which parameters sections you want
included in the file, both for those
:py:`Layer`\s and :py:`Species` and for
other components of the :py:`Model`; and the path and filename for your new
parameters file). Geonomics will then automatically create the file for you, 
arranged as you requested and saved where you requested.

When you then open that file, you will see the following:

.. code-block:: python

  #<your_filename>.py

  #This is a default parameters file generated by Geonomics
  #(by the gnx.params.make_parameters_file() function).
  
  
                        ## :: ::    :::            ##
                  ##:::   :::::    :::   ::    :: :: :::##
               ## ::::     ::           ::   ::::::::::::::##
             ##::::::                       ::::::::: :::::: :##
           ## :    :::                    :::    ::    :::::::::##
          ##ggggg eeee ooo   n   n   ooo   m   m iiiii  cccc ssss##
         ##g     e    o   o  nn  n  o   o  m   m   i   c     s    ##
         ##g     eee o     o n n n o     o mm mm   i   c     sssss##
         ##g ggg eee o     o n  nn o     o m m m   i   c         s##
         ##g   g e    o   o  n   n  o   o  m   m   i   c        ss##
          ##gggg  eeee ooo   n   n   ooo   m   m iiiii  cccc ssss##
           ##  ::::::::        ::::::::::::  :       ::  ::   : ##
             ##  ::::              :::::::  ::     ::::::::  :##
               ## :::               :::::: ::       ::::::  ##
                  ##:                ::::                ##
                        ##                         ##
  
  
  params = {
  
  ##############
  #### LAND ####
  ##############
      'land': {
  
      ##############
      #### main ####
      ##############
          'main': {
              # y,x (a.k.a. i,j) dimensions of the Landscape
              'dim':                      (20,20),

     #.
     #.
     #.

This is the beginning of a file that is really just a long but simple Python
script (hence the '.py' extension); this whole file just defines a single,
long, nested :py:`dict` (i.e. a Python 'dictionary') containing all of your
parameter values. It may look like a lot, but don't be concerned! For two
reasons:

  1. All the hard work is already done for you. You'll just need to change
     the default values where and how you want to, to set up your particular
     simulation scenario.
  2. You will probably leave a good number of the parameters defined in this
     file untouched. Geonomics does its best to set sensible default values
     for all its parameters. Though of course, you'll want to think clearly 
     nonetheless about whether the default value for each parameter 
     is satisfactory for your purposes.

Each parameter in the parameters value is preceded by a terse comment, to
remind you what the parameter does. But for detailed information about each
parameter, you'll want to refer to the following information.
What follows is a list of all of the Geonomics parameters (in the sections and
the top-to-bottom order in which they'll appear in your parameters files).
For each parameter, you will see a section with the following information:

  - a snippet of the context (i.e. lines of
    Python code) in which it appears in a parameters file; 
  - the valid Python data type(s) the parameter can take
  - the default value of the parameter
  - a ranking score, indicating how likely it is that you will want to reset
    this parameter (i.e. change it from its default value), and
    encoded as follows:

    - 'Y': almost certainly, *or* must be reset for your :py:`Model` to run
    - 'P': it is quite possible that you will want to reset this
      parameter, but this will depend on your use and scenario
    - 'N': almost certainly not, *or* no need to reset because it should be
      set intelligently anyhow (Note: this does *not* mean that you cannot
      reset the parameter! if that is the case for any value then it does not
      appear in the parameters file)

  - other relevant, detailed information about the parameter, including
    an explanation of what it defines, how its value is used, where to look
    for additioanl information about parameters related to other Python 
    packages, etcetera
   

These section will be formatted as follows:


**<param_name>**

.. code-block:: python

              #brief comment about the parameter
              '<param_name>':               <default_param_value>,

<valid Python data type(s)>

default: <default value>

reset? <ranking>

  <Explanation of what the parameter defines, how its value is used,
  and any other relevant information.>


This section should serve as your primary point of reference
if you confront any uncertainty while creating your own parameters files.
We'll start with the section of parameters that
pertains to the :py:`Landscape` object.


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


====================
Community parameters
====================

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
trait names, and effect sizes of all loci). Geonomics will create an empty
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


================
Model parameters
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
                
 

-------------------------------------------------------------------------------


*****************************
Class and function docstrings
*****************************

==============================
:py:`gnx.make_parameters_file`
==============================


Create a new parameters file.

Write to disk a new, template parameters file. The file will contain the
numbers and types of sections indicated by the parameters fed to this
function. It can often be used 'out of the box' to make a new Model
object, but typically it will be edited by the user to stipulate
the scenario being simulated, then used to instantiate a Model.

----------
Parameters
----------
filepath : str, optional
    Where to write the resulting parameters file, in /path/to/filename.py
    format. Defaults to None. If None, a file named
    "GNX_params_<datetime>.py" will be written to the working
    directory.
scapes : {int, list of dicts}, optional
    Number (and optionally, types) of Layer-parameter sections to include
    in the parameters file that is generated. Defaults to 1. Valid values
    and their associated behaviors are:

    int:
        Add sections for the stipulated number of Layers, each with default
        settings:
        
          - parameters for creating Layers of type 'random' (i.e.
            Layers that will be generated by interpolation from
            randomly valued random points)
          - no LayerChanger parameters

    [dict, ..., dict]:
        Each dict in this list should be of the form:

        {'type':    'random', 'defined', 'file', or 'nlmpy',

        'change':   bool

        }

        This will add one section of Layer parameters, with the
        contents indicated, for each dict in this list.
species : {int, list of dicts}, optional
    Number (and optionally, types) of Species-parameter sections to
    include in the parameters file that is generated. Defaults to 1. Valid
    values and their associated behaviors are:

    int:
        Add sections for the stipulated number of Species, each with
        default settings:

          - parameters for movement without a MovementSurface
          - parameters for a GenomicArchitecture with 0 Traits (i.e. with
            only neutral loci)
          - no SpeciesChanger parameters

    [dict, ..., dict]:
        Each dict should contain at least one argument from among the
        following:
        {'movement':                       bool,
        'movement_surface':                bool,
        'genomes':                         {bool, 'custom'},
        'n_traits':                        int,
        'demographic_change':              int,
        'parameter_change':                bool
        }
        This will add one section of Species parameters, customized
        as indicated, for each dict in the list. (Note that if the
        'genomes' argument is True or 'custom', a section for
        parameterization of the genomic architecture will be added,
        but if it is 'custom' then a template custom genomic architecture
        file will also be created (a CSV file), which can be filled in
        to stipulate the locus-wise values for starting allele frequency,
        recombination rate, dominance, associated traits, and effect
        sizes.)

data : bool, optional
    Whether to include a Data-parameter section in the parameters file that
    is generated. Defaults to None. Valid values and their associated
    behaviors are:

    None, False:
        Will not add a section for parameterizing data to be collected.
        No DataCollector will be created for the Model object made from
        the resulting parameters file, and no data will be collected
        during the model runs.
    True:
        Will add a section that can be used to parameterize which
        data will be collected during the model runs, when, and what
        file formats will be used to write it to disk.
        (This which will be managed by the model's DataCollector
        object.)

stats : bool, optional
    Whether to include a Stats-parameter section in the parameters file that
    is generated. Defaults to None. Valid values and their associated
    behaviors are:

    None, False:
        Will not add a section for parameterizing the statistics to be
        calculated. No StatsCollector will be created for the Model
        object made from the resulting parameters file, and no
        statistics will be calculated during the model runs.
    True:
        Will add a section that can be used to parameterize which
        statistics will be calculated during the model runs, and when.
        (This will be managed by the model's StatsCollector object.)

-------
Returns
-------
out : None
    Returns no output. Resulting parameters file will be written to the
    location and filename indicated (or by default, will be written to a
    file named "GNX_params_<datetime>.py" in the working directory).

--------
See Also
--------
sim.params.make_parameters_file

-----
Notes
-----
All parameters of this function are optional. Calling the function without
providing any parameters will always produce the parameters file for the
default model scenario. This file can be instantiated as a Model object and
run without being edited. Those three steps (create default parameters file;
create model from that parameters file; run the model) serve as a base case
to test successful package installation, and are wrapped around by the
convenience function `gnx.run_default_model`.

--------
Examples
--------
In the simplest example, we can create a parameters file for the default
model. Then (assuming it is the only Geonomics parameters file in the
current working directory, so that it can be unambiguously identified) we
can call the gnx.make_model function to create a Model object from that
file, and then call the Model.run method to run the model (setting the
'verbose' parameter to True, so that we can observe model output).

>>> gnx.make_parameters_file()
>>> mod = gnx.make_model()
>>> mod.run(verbose = True)
TODO: PUT TYPICAL MODEL OUTPUT HERE, EVEN THOUGH IT'S ONLY PRINTED?

We can use some of the function's arguments, to create a parameters
file for a model with 3 Layers and 1 Species (all with the default
components for their sections of the parameters file) and with a section
for parameterizing data collection.

>>> gnx.make_parameters_file(scapes = 3, data = True)

As a more complex example that is likely to be similar to most use cases,
we can create a parameters file for a model scenario with:

    - 2 Layers (one being an nlmpy Layer that will not change over model
      time, the other being a raster read in from a GIS file and being
      subject to change over model time);
    - 2 Species (the first having genomes, 2 Traits, and movement
      that is dictated by a MovementSurface; the second not having
      genomes but having a MovementSurface as well, and undergoing
      demographic change)
    - data-collection;
    - stats-collection;

We can save this to a file named "2-spp_2-trait_model.py" in our current
working directory.

>>> gnx.make_parameters_file(
>>>     #list of 2 dicts, each containing the values for each Layer's
>>>     #parameters section
>>>     scapes = [
>>>         {'type': 'nlmpy'},                              #scape 1 
>>>         {'type': 'gis',                                 #scape 2 
>>>          'change': True}
>>>         ],
>>>     #list of 2 dicts, each containing the values for each Species'
>>>     #parameters section
>>>     species = [
>>>         {'genomes': True,                               #spp 1
>>>          'n_traits': 2,
>>>          'movement': True,
>>>          'movement_surface': True},
>>>         {'genomes': False,                              #spp 2
>>>          'movement': True,
>>>          'movement_surface': True,
>>>          'demographic_change': True}
>>>         ],
>>>     #arguments to the data and stats parameters
>>>     data = True, stats = True, 
>>>     #destination to which to write the resulting parameter file
>>>     filepath = '2-spp_2-trait_model.py')





==========================
:py:`read_parameters_file`
==========================

Create a new ParametersDict object.

Read the Geonomics parameters file saved at the location indicated by
'filepath', check its validity (i.e. that all the Layers and Species
parameterized in that file have been given distinct names), then use the
file to instantiate a ParametersDict object.

----------
Parameters
----------
filepath : str
    String indicating the location of the Geonomics parameters file that
    should be made into a ParametersDict object.

-------
Returns
-------

An object of the ParametersDict class (a dict of nested dicts, all
of which have key-value pairs whose values can be accessed using typical
dict notation or using dot notation with the keys).

------
Raises
------
AssertionError
    If either the Layers or the Species parameterized in the parameters
    file have not all been given distinct names

--------
See Also
--------
sim.params.read
sim.params.ParametersDict

--------
Examples
--------
Read a parameters file called "null_model.py" (located in the current
working directory).

>>> gnx.read_parameters_file('null_model.py')
<class 'sim.params.ParametersDict'>
Model name:                                     GNX_params_13-10-2018_15:54:03



================
:py:`make_model`
================

Create a new Model object.

Use either a ParametersDict object or the path to a valid Geonomics
parameters file (whichever is provided to the 'parameters' argument) to
create a new Model object.

----------
Parameters
----------
parameters : {ParametersDict, str}, optional
    The parameters to be used to make the Model object.
    If `parameters` is a ParametersDict object, the object will be used to
    make the Model.
    If `parameters` is a string, Geonomics will call
    `gnx.read_parameters_file` to make a ParametersDict object, then use
    that object to make the Model.
    If `parameters` is None, or is not provided, then Geonomics will
    attempt to find a single parameters file in the current working
    directory with the filename "GNX_params_<...>.py", will use that
    file to make a ParametersDict object, then will use that object to
    make the Model.

-------
Returns
-------
out : Model
    An object of the Model class

------
Raises
------
ValueError
    If the `parameters` argument was not provided and a single, valid
    Geonomics parameters file could not be identified in the current
    working directory
ValueError
    If the `parameters` arugment was given a string that does not point
    to a valid parameters file
ValueError
    If the ParametersDict provided to the `parameters` argument, or created
    from the parameters file being used, cannot be successfully made into a
    Model

--------
See Also
--------
gnx.read_parameters_file
sim.model.Model

--------
Examples
--------
Make a Model from a single, valid "GNX_params_<...>.py" file that can
be found in the current working directory (such as a file that would be
produced by calling gnx.make_parameters_file without any arguments).

>>> gnx.make_model()
<class 'sim.model.Model'>
Model name:                                     GNX_params_13-10-2018_15:54:03
Layers:                                         0: '0'
Species:                                        0: '0'
Number of iterations:                           1
Number of burn-in timesteps (minimum):          30
Number of main timesteps:                       100
Geo-data collected:                             {}
Gen-data collected:                             {}
Stats collected:                                {}


Make a Model from a file called 'null_model.py', in the current working
directory.

>>> gnx.make_model('null_model.py')
<class 'sim.model.Model'>
Model name:                                     null_model
Layers:                                         0: 'tmp'
                                                1: 'ppt'
Species:                                        0: 'C. fasciata'
Number of iterations:                           2500
Number of burn-in timesteps (mininum):          100
Number of main timesteps:                       1000
Geo-data collected:                             {csv, geotiff}
Gen-data collected:                             {vcf, fasta}
Stats collected:                                {maf, ld, mean_fit, het, Nt}
