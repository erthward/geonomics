#########
Geonomics
#########


:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
A Python package for easy construction of individual-based, spatially explicit,
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
forward-time, and highly customizable landscape-genomics simulations
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

*****
Intro
*****

Backward-time (i.e. coalescent) simulators abound.
But they are inadequate for simulation of many scenarios of 
interest, including: natural selection on traits with arbitrary genomic 
architectures; spatially variable natural selection; simulation of populations
distributed continiously and moving realistically across
complex landscapes; complex demographic change simultaneous with ongoing, 
often non-stationary environmental change; and coevolutionary interactions 
between multiple species or incpient species. Few existing forward-time 
simulators can model all of these phenomena, and the few that can often 
impose a high cost of entry (e.g. learning a new, non-standard programming
language in order to write one's desired simulations). Geonomics aims to fill 
this empty niche by combining ease of use with , all in a popular, 
general-purpose scripting language. If it succeeds at doing this, Geonomics 
should prove uniquely useful for a wide range of purposes, from intro-level 
educational use to high-quality theoretical, methogological, empirical, and
applied research.

Geonomics is written in Python. Python is a full-fledged scripting language 
that is relatively easy to learn (and fun!). So it can be pretty quick for a
new user to get up to speed and start doing useful work. For work with
Geonomics, this turnaround time should be even quicker. Geonomics aims to
require minimal Python knowledge (yet maintain high extensibility for
interested, expert users). Essentially, anyone should be able to build their
own, arbitrarily complex Geonomics models as long as they know how to install
the package, open a Python console, call Python functions, and edit some
default values in a pre-packaged script. 

For proper acquaintance with the model, however, it will be helpful to 
understand the concept of **object-oriented programming**.  Python is a very 
handy language for object-oriented programming, and this is the primary 
programming paradigm in which Geonomics is written. Essentially, 
object-orientation involves: 

  1. Defining certain types of data structures, or **classes** (e.g. `Car`),
     and assigning them various behaviors, or **methods** (e.g. `honk`, 
     `turn`, `brake`);
  2. Using particular data values to create individual instances, or 
     **objects** belonging to those classes (e.g. `my_1986_jeep`, `batmobile`);
  3. Instructing those **objects** to carry out their behaviors by 'calling' 
     their **methods** (e.g. `my_1986_jeep.honk()` might return "Beepity
     beep!", wheras `batmobile.honk()` might return "<Batman theme song>"). 
     
Geonomics defines a number of **classes**, some with a number of **methods**
designed specifically for the user. The user will create **objects**
belongining to these classes, using their own particular data values
(i.e. with the values they define in the paramaters files they create). Then
the user will call key **methods** on these **objects**, in order to get them
to carry out their behaviors (that is to say, to get them to run a
simulation).


The subsequent documentation will present the **classes** definined in
Geonomics and their key **methods**. It will explain exactly what those methods
do, and how what they do fits into the overall structure and function of 
Geonomics models.

In brief, Geonomics provides the following components:

  1. `Land`, `GenomicArchitecture`, `Trait`, `Individual`, `Population`,
     and `Community` classes;
  2. A `Model` class that builds objects of the aforementioned classes 
     according to the scenario stipulated in a parameters file,
     uses them to run numerous simulation iterations, and saves data 
     and basic statistics from those iterations at all desired timesteps;
  3. Classes for planning arbitrarily complex landscape- and 
     demographic-change scenarios;
  4. Numerous visualization methods, to help users design, run, explore, 
     present, and explain their models' behavior and results.

**For the beginner**, we recommend the following steps:
  1. Review the following three sections ('Model organization', 'Data
     structures and classes', and 'Operations'), to get a general
     undertsanding of the logic, components, and necessary and optional
     behaviors of a Geonomics model;
  2. Skim the subsequent section ('Parameters'), to understand the structure
     of a Geonomics parameters file;
  3. Pip-install Geonomics (`$ pip install geonomics`);
  4. Open Python and run `import geonomics as gnx`;
  5. Use the `gnx.make_parameters_file` function, to begin creating template
     parameters files that they can tweak as desired;
  6. Use the `gnx.make_model` function and then the `Model.walk` or `Model.run`
     methods to instantiate and run the models they've parameterized;
  7. Use the various `Model.plot` methods to visualize the behavior and results
     of their models.

**For the *impatient* beginner**, as soon as Geonomics has been
installed, you should be able to run the following code:


.. code-block:: python

     import geonomics as gnx

     gnx.run_default_model()

This will load the Geonomics package, create a default Geonomics
parameters file in your current working directory, 
then use that file to instantiate and run a `Model` using the default
parameter values.

Lastly, please note that this documentation is designed to be read from
from the beginning (though not necessarily all the way to the end). But
given the interrelationships between all the components of the package,
you will inevitably run into material in different sections that is related.
To the extent possible, we attempt to cross-reference rather than duplicate
information.

Merry modeling!


******************
Model organization 
******************

<<<DIAGRAM HERE>>>


:: +-------+        +-------+ 
:: | stuff | -----> | stuff | 
:: +-------+        +-------+ 
::     ^                |     
::     |                v     
:: +-------+        +-------+ 
:: | stuff | <----- | stuff | 
:: +-------+        +-------+ 
         
***************************
Data structures and classes
***************************

The following sections discuss the structure and function of the key
Geonomics classes. Users will interface with these classes more or less
directly when running Geonomics models, so a fundamental understanding of how 
they're organized and how they work will be useful.

==========================
`Land` and `Scape` objects
==========================

One of the core components of a Geonomics model is the land. The land is
modeled by the `Land` class. This class is an integer-keyed `dict` composed of
numerous instances of the
class `Scape`. Each `Scape` represents a separate environmental variable (or
'layer', in GIS terminology), which is modeled a 2d Numpy array (or raster; in
attribute 'rast'), of identical dimensions to each other `Scape` in the `Land`
object, and with the values of its environmental variable 'e' constrained to
the interval [0 <= e <= 1]. Each `Scape` can be initialized from its own
parameters subsection within the 'land' parameters section of a Geonomics
parameters file. 

For each `Population` (see section '`Individuals` and `Populations`', below),
the different `Scape` layers in the `Land` can be used to model habitat 
viability, habitat connectivity, or variables imposing spatially varying
natural selection. `Land` and `Scape` objects also contain some metatdata 
(as public attributes), including the resolution (attribute `res`), upper-left
corner (`ulc`), and projection (`prj`), which default to 1, (0,0), and None but
will be set otherwise if some or all of the `Scape` layers are read in from
real-world GIS rasters.

===================================================
Genomes, `GenomicArchitecture`, and `Trait` objects
===================================================

`Individual` objects (see section '`Individuals` and `Populations`', below) can
optionally be assigned genomes. If they are, each `Individual`'s genome is 
modeled as a 2-by-L Numpy array (where 2 is the ploidy, currently fixed at
diploidy, and L is genome length) containing 0s and 1s (because
Geonomics strictly models diallelic SNPs). 

The parameter L, as well as numerous other genomic parameters (including 
locus-wise starting frequencies of the 1 alleles; locus-wise dominance effects;
locus-wise recombination rates; and genome-wide mutation rates for neutral, 
globally deleterious, and adaptive loci), are controlled by the 
`GenomicArchitecture` object pertaining to the `Population` to which an 
`Individual` belongs. (For the full and detailed list of attributes in a 
`GenomicArchitecture` object, see its class documentation, below.)
The genomes of the initial `Individual`\s in a simulation, as well as those of 
`Individual`\s in subsequent generations, are either drawn or recombined, and
are mutated, according to the values stipulated by the `GenomicArchitecture` of
their `Population`. The user can create a population with a 
`GenomicArchitecture` and with corresponding genomes by including a 'genome' 
subsection in that population's section of the Geonomics parameters file (and 
setting the section's various parameters to their desired values). 

Geonomics can model `Individual`\s' phenotypes. It does this by allowing the 
user to create an arbitrary number of distinct `Trait`\s
for each `Population`. Each trait is represented by a `Trait` object, which 
maps genomic loci onto that trait, maps effect sizes ('alpha') onto those loci,
and sets the trait's polygenic selection coefficient ('phi'). An `Individual`'s
phenotype for a given trait is calculated as the 'null phenotype' plus a 
weighted sum of the products of its 'effective genotypes' at all loci 
underlying that `Trait` and the effect sizes (i.e. 'alpha') of those loci: 
:math:`z_{i,t} = null\_genotype + \sum_{l = 0}^{n} \alpha_{t,l} g_{i,l}`, 
where :math:`z_{i,t}` is the phenotype of `Individual` i for trait t, 
:math:`g_{i, l}` is the genotype of the `Individual` at that locus, and 
:math:`\alpha_{t,l}` is the effect size of that locus for that trait.

The 'null phenotype' refers determines what would be the phenotypic value that
an `Individual` who is homozygyous for the 0 allele at all loci for a trait.
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

(For the full and detailed list of attributes in a `Trait` object, 
see its class documentation, below.)

Note that for maximal control over the `GenomicArchitecture`
of a `Population`, the user can set the value of the 'gen_arch_file' 
parameter in the parameters file to the name of a separate CSV file 
stipulating the locus numbers, starting 1-allele frequencies, dominance 
effects, traits, and inter-locus recombination rates (as columns) of 
all loci (rows) in the `GenomicArchitecture`;
these values will override any other values provided in the 'genome' 
subsection of the population's parameters.

===================================================
`Individual`, `Population`, and `Community` objects
===================================================

Being that Geonomics is an individual-based model, individuals serve as 
the fundamental units (or agents) of all simulations. They are represented by
objects of the `Individual` class. Each `Individual` has an index (saved 
as attribute 'idx'), a sex (attribute 'sex'), an age (attribute 'age'), 
an x,y position (in continuous space; attributes 'x' and 'y'), and a 
`list` of environment values (attribute 'e'), extracted from the 
`Individual`'s current cell on each `Scape` of the `Land` on which 
the `Individual` lives.

The `Population` class is an `OrderedDict` (defined by the `collections` 
package) containing all `Individaul`\s, (with 
their 'idx' attributes as keys). If a `Population` has a `GenomicArchitecture`
then the `Individual`\s in the `Population` will also each have genomes 
(attribute 'genome'),
and the `GenomicArchitecture` includes `Trait`\s then each individual will also
have a `list` of phenotype values (one per `Trait`; attribute 'z') and a 
single fitness value (attribute 'fit'). (These attributes all otherwise 
default to `None`.)

Each `Population` also has a number of other attributes of interest. Some 
of these are universal (i.e. they are created regardless of the 
parameterization of the `Model` to which a `Population` inheres). These 
include: the `Population`'s name (attribute 'name'); its current density 
raster (a Numpy array attribute called 'N'); and the number of births,
number of deaths, and terminal population size of each timestep (which are 
`list` attributes called 'n_births', 'n_deaths', and 'Nt'). If the 
`Population` was parameterized with a `GenomicArchitecture` then that will 
be created as the 'gen_arch' attribute (otherwise this attribute will be 
`None`).

All of the `Population`\s in a `Model` are collected in the `Model`'s 
`Community` object. The `Community` class is simply an integer-keyed `dict` 
of `Population`\s. For the time being, the `Community` object allows a 
Geonomics `Model` to simulate multiple `Population`\s simultaneously on 
the same `Land`, but otherwise affords no additional functionality
of interest. However, its implementation will facilitate the potential 
future development of methods for interaction between `Population`\s. 
(e.g. to simulate coevolutionary, speciation, or hybridization scenarios).


=======
`Model`
=======

Objects of the `Model` class serve as the main interface between the user 
and the Geonomics program. (While it is certainly possible for a user 
to work directly with the `Land` and `Population` or `Community` objects to 
script their own custom models, the typical user should find that the 
`Model` object allows them accomplish their goals with minimal toil.)
The main affordance of a `Model` object is the `Model.run` method, 
which, as one could guess, will run the `Model`. The typical workflow 
for creating and running a  `Model` object is as follows:

  1. Create a template paramters file containing the desired sections, 
     by calling `gnx.make_parameters_file` with all revelant arguments;
  2. Define the scenario to be simulated, by opening and editing that 
     parameters file (and optionally, creating/editing corresponding 
     files, e.g. genomic-architecture CSV files;
     or raster or numpy-array files to be used as `Scape`\s);
  3. Instantiate a `Model` object from that parameters file, by calling 
     `mod = gnx.make_model('/path/to/params_filename.py')`;
  4. Run the `Model`, by calling `mod.run()`.

For detailed information on usage of these functions, see their docstrings.
When a `Model` is run, it will:

  1. Run the burn-in (until the mininmal burn-in length stipulated in the 
     parameters file and the built-in stationarity statistics 
     determine that the burn-in is complete);
  2. Run the main model for the stipulated number of timesteps;
  3. Repeat this for the stipulated number of iterations (retaining or 
     refreshing the first run's initial `Land` and `Population` 
     objects as stipulated).

The `Model` object offers one other method, however, `Model.walk`, 
which allows the user to run a model, in either 'burn' or 'main' mode, 
for an arbitrary number of timesteps within a single iteration (see its 
docstring for details). This is particularly useful for running 
Geonomics within an interactive Python session. Thus, `Model.walk` is 
primarily designed for passively running numerous iterations of a `Model`, 
to generate data for analysis, whereas `Model.walk` is primarily designed
for the purposes of learning, teaching, or debugging the package, or 
developing, exploring, introspecting, or visaulizing particular `Model`\s. 

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

------------------
`_MovementSurface`
------------------

The `_MovementSurface` class allows Geonomics to model a `Population`'s 
realistic movement across a spatially varying landscape. It does this by 
creating an array of circular probability distributions (i.e. VonMises 
distributions), one for each cell on the `Land`, from which 
`Individual`\s choose their directions each time they move. To create the
`_MovementSurface` for a `Population`, the user must indicate the `Scape` 
that should be used to create it (i.e. the `Scape` that represents 
landscape permeability for that `Population`). The `_MovementSurface`'s 
distributions can be **simple (i.e. unimodal)**, such that the 
maximum value of the distribution at each cell will point toward the
maximum value in the 8-cell neighborhood; this works best for permeability 
`Scape`\s with shallow, monotonic gradients, because the differences 
between permeability values of neighboring cells can be minor (e.g. a 
gradient representing the directionality of a prevalent current). 
Alternatively, the distributions can be **mixture (i.e. multimodal)**
distributions, which are weighted sums of 8 unimodal distributions, one 
for each neighboring cell, where the weights are the relative cell 
permeabilities (i.e. the relative probabilities that an `Individual` would 
move into each of the 8 neighboring cells); this works best for non-monotonic, 
complex permeability `Scape`\s (e.g. a DEM of a mountainous region that is 
used as a permeability `Scape`). (The `Land` is surrounded by a margin 
of 0-permeability cells before the `_MovementSurface` is calculated, such 
that `Land` edges are treated as barriers to movement.) The class consists 
principally of a 3d Numpy array (x by y by z, where x and y are the 
dimensions of the `Land` and z is the length of the vector of values 
used to approximate the distributions in each cell.

-------------------
`_DensityGridStack`
-------------------

The `_DensityGridStack` class implements an algorithm for rapid estimating 
an array of the local density of a `Population`. The density is estimated 
using a sliding window approach, with the window-width determining the 
neighborhood size of the estimate. The resulting array has a spatial 
resolution equivalent to that of the `Land`, and is used in all
density-dependent operations.

---------
`_KDTree`
---------

The `_KDTree` class is just a wrapper around `scipy.spatial.cKDTree`. 
It provides an optimized algorithm (the kd-tree) for finding 
neighboring points within a given search radius.
This class is used for all neighbor-searching operations (e.g. mate-search).

---------------------
`_RecombinationPaths`
---------------------

The `_RecombinationPaths` class contains a large (and customizable) 
number of `bitarray`\s, each of which indicates the genome-length 
diploid chromatid numbers (0 or 1) for a
recombinant gamete produced by an `Individual` of a given `Population` 
(henceforth referred to as 'recombination paths'). These recombination 
paths are generated using the genome-wide recombination rates specified by 
the `Population`'s `GeonomicArchitecture`. They are generated during 
construction of the `Model`, then drawn randomly as needed (i.e.
each time an `Individual` produces a gamete). This provides a 
reasonable trade-off between realistic modelling of recombination and runtime.

---------------------------------------
`_LandChanger` and `_PopulationChanger`
---------------------------------------

These classes manage all of the landscape changes and demographic changes 
that were parameterized for the `Land` and `Population` objects to which 
they inhere. The functions creating these changes are defined at the outset, 
then queued and called at their scheduled timesteps.

--------------------------------------
`_DataCollector` and `_StatsCollector`
--------------------------------------

These classes manage all of the data and statistics that should be collected 
and written to file for the `Model` object to which they inhere 
(as determined by the parameters file used the create the `Model`). 
The types of data to be collected, or statistics to be calculated, as 
well as the timesteps at which and methods by which they're 
collected/calculated and determined at the outset, then the 
appropriate functions called at the appropriate timesteps.

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
for a `Population` with a maximum age of 1). For `Population`\s 
with movement, `Individual`\s can move by two distinct mechanisms. **Spatially
random movement** is the default behavior; in this case, `Individual`\s 
move to next locations that are determined by a random distance drawn 
from a Wald distribution and a random direction drawn from a uniform 
circular (i.e. Von Mises) distribution.  As with most distributions used 
in Geonomics, the parameters of these distributions have sensible 
default values but can be customized in a `Model`'s parameters file 
(see section 'Parameters', below). 

The alternative movement mechanism that is available is 
**movement across a permeability surface**, using a `_MovementSurface` object.
To parameterize a `_MovemementSurface` for a `Population`, the user 
must create a template parameters file that includes the 
necessary parameters section for the population (i.e. 
the user must set 'movement' to `True` and 'movement_surface' to `True` 
in the population's arguments to the `gnx.make_parameters_file` 
function (see the docstring for that function for details and an example). 
`Individual`\s move to next locations determined by a random distance drawn 
from a Wald distribution and a random direction drawn from the distribution 
at the  `_MovementSurface` cell in which which the `Individual`\s 
are currently located. For details about `_MovementSurface` creation, see 
section '`_MovementSurface`' above, or the class' docstring.

Dispersal is currently implemeneted identically to spatially random movement 
(with the caveat that the an offspring's new location is determined 
relative its parents' centroid). But the option to use a 
`_MovementSurface` for dispersal will be offered soon.

============
Reproduction
============

Each timestep, for each `Population`, all pairs of individuals within 
a certain distance of each other (i.e. the mating radius, 
which is set in the parameters file) are identified.
These pairs are subsetted if necessary (i.e. if the `Population` 
requires that `Individual`\s be above a certain reproductive age, 
or that they be of opposite sexes, in order to mate; these values 
can also be changed from their defaults in the parameters file). 
Remaining pairs mate probabilistically (according to a Bernoulli 
random draw with probability equal to the `Population`'s birth 
rate, which is also set in the parameters file).

Pairs that are chosen to mate will produce a number of new 
offspring drawn from a Poisson distribution (with lambda set in the 
parameters file). For each offspring, sex is chosen probablistically 
(a Bernoulli random draw with probability equal to the `Population`'s 
sex ratio), age set to 0, and location chosen by dispersal from 
the parents' centroid (see section 'Movement and Dispersal'). For 
`Population`\s that have genomes, offspring genomes will be a 
fusion of two recombinant genomes from each of the two parents (where 
each recombinant is indexed out a parent's genome using a recombination 
path; see section '`_RecombinationPaths`'). For `Population`\s 
with `Trait`\s in their `GenomicArchitecture`\s, offspring phenotypes are 
determined at birth. Mutations are also drawn and introduced at this 
point (see section 'Mutation for details).

=========
Mortality
=========

Mortality can occur as a combination of two factors: **density dependence** 
and **natural selection**. Each `Individual` has a death decision drawn 
as a Bernoulli random variable with 
:math:`P(d_{i}) = 1 - P(s_{i_{dens}})P(s_{i_{fit}})`, where :math:`P(d_{i})` 
is the probability of death of `Individual` :math:`i`, and 
:math:`P(s_{i_{dens}})` and :math:`P(s_{i_{fit}})` are the probabilities of 
survival of `Individual` :math:`i` given density-dependence and 
fitness. The probability of density-dependent death is contingent on an 
`Individual`'s x,y location (i.e. the cell in which they're currently located. 
And an `Individual`'s probability of survival due to fitness 
is just equal to the product of their absolute fitness (:math:`\omega`) 
for each of the `Individual`'s :math:`m` `Trait`\s. 
Thus the equation for an `Individual`'s probability of death becomes:

.. math::
   P(d_{i}) = 1 - (1 - P(d_{x,y})) \prod_{p = 1}^{m}\omega_{i,p}

The following two sections explain in detail the implementation and 
calculation of the two halves of the right side of this equation.

------------------
Density dependence
------------------

Density dependence is implemented using a spatialized form of the class 
logistic growth equation 
(:math:`\frac{\mathrm{d}N_{x,y}}{\mathrm{d}t}=rN_{x,y}(1-\frac{N_{x,y}}{K_{x,y}})`, 
where the x,y subscripts refer to values for a given cell on the `Land`).
Each `Population` has a carrying-capacity raster (a 2d Numpy array; 
attribute 'K'), which is defined in the parameters file to be 
one of the `Scape`\s in the `Land`. The comparison between this raster and 
the population-density raster calculated at each timestep serves as the 
basis for the spatialized logistic growth equation, because both 
equations can be calculated cell-wise for the entire extent of the 
`Land` (using the `Population`'s intrinsic growth rate, the attribute 
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

Selection on a `Trait` can exhibit three modes: **spatially divergent**, 
**universal**, and **spatially contingent**. **Spatially divergent** selection 
is the default behavior, and the most commonly used; in this form of 
selection, an `Individual`'s fitness depends on the absolute difference 
between the `Individual`'s phenotypic value and the environmental
value of the relevant `Scape` (i.e. the `Scape` that represents the 
environmental variable acting as the selective force) in the cell where 
the `Individual` is located.

**Universal** selection (which can be toggled using the 'univ_adv' 
parameter with a `Trait`'s section in the parameters file) occurs 
when a phenotype of 1 is optimal everywhere on the `Land`. In other 
words, it represents directional selection on an entire `Population`,
regardless of `Individual`\s' spatial contexts. 

Under **spatially contingent** selection, the selection coefficient of a 
`Trait` varies across space, such that the strength of selection 
is environmentally determined in some way. Importantly, this mode 
is *not mutually exclusive* with the other two; in other words, 
selection on a certain `Trait` be both spatially contingent 
and either spatially divergent or universal. Spatially contingent selection 
can be implemented by providing an array of values (equal in dimensions 
to the `Land`) to the 'phi' value of a `Trait`, rather than a scalar 
value (which could be done within the parameters file itself, but may be 
more easily accomplished as a step between reading in a parameters file and 
instantiating a `Model` object from it). 

All possible combinations of the three modes of selection can all 
be thought of as special cases of the following equation for the fitness of 
`Individual` :math:`i` for `Trait` :math:`p` (:math:`\\omega_{i,p}`):

.. math::
   \omega_{i,p}= 1 - \phi_{p} (\mid e_{p;x,y} - z_{i;p} \mid)^{\gamma_{p}}

Where :math:`\\phi_{p}` is the selection coefficient of trait 
:math:`p`, :math:`e_{p;x,y}` is the environmental variable of the 
relevant `Scape` at `Individual` :math:`i`'s x,y location
(which can also be thought of as the `Individual`'s optimal 
phenotype), :math:`z_{i;p}` is `Individual` :math:`i`'s (actual) 
phenotype for `Trait` :math:`p`, and :math:`gamma_{p}` controls 
how fitness decreases as the absolute difference between an `Individual`'s 
optimal and actual phenotypes increases (it defaults to 1, which causes 
fitness to decrease linearly around the optimal 
phenotypic value). 

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
  identical if a `Trait` in monogenic), and it is (by default, 
  but not always; see section 'Selection', above) divergent. For this reason
  it would be a misnomer to call mutations that influence a given 
  `Trait`'s phenotypes 'beneficial' -- even though that term is the closest
  population-genetic concept to this concept as it is employed in Geonomics -- 
  because the same mutant genotype in the same `Individual` could have opposite
  effects on that `Individual`'s fitness in different environmental 
  contexts (i.e. it could behave as a beneficial mutation is one region of
  the `Land` but as a deleterious mutation in another). 


======================
Population interaction
======================

This functionality is not yet included available. But the Community class was 
created in advance recognition that this functionality could be desirable 
for future versions (e.g. to simulate coevolutionary, speciation, or 
hybridization scenarios).


==========================
Land and population change
==========================

For a given `Scape`, any number of landsacpe change events can be planned. 
In the parameters file, for each event, the user stipulates the initial
timestep; the final timestep; the resulting landscape (i.e. the array 
of the `Scape` that will exist after the event is complete); and the 
interval at which intermediate changes will occur.  When the `Model` is 
created, the stepped series of intermediate landscapes (and 
`_MovementSurface` objects, if the `Scape` that is changing serves 
as the basis of a `_MovementSurface` for any `Population`) will be 
created and queued, so that they will swap out accordingly at the appropriate 
timesteps.

For a given `Population`, any number of demographic change events can 
also be planned. In the parameters file, for each event, the user 
stipulates the type of the event ('monotonic', 'cyclical', 'random', or 
'custom') as well as the values of a number of associated 
parameters (precisely which parameters depdends on the type of event chosen).
As with landscape change events, all necessary stepwise changes will be 
planned and queued when the `Model` is created, and will be 
executed at the appropriate timesteps.

It is also possible to schedule any number of instantaneous changes 
to some of the life-history parameters of a `Population` (e.g. birth rate; 
the lambda parameter of the Poisson distribution determining the number of 
offspring of mating events). This functionality is currently minimalistic, 
but will be more facilitated in future versions.



*************
Visualization
*************

Each `Population` has a wide variety of visualization methods 
(`Population.plot`, `Population.plot_fitness`, etc.), which aim to help users 
design, run, explore, present, and explain their models' behavior and results.
These methods can be called on a `Population` at any time (e.g. as 
soon as the `Population` has been created, or after the model has
run for any number of timesteps); but it is worth mentioning that some 
methods may be invalid depending on the point in model-time at 
which they're called (e.g.  `Population.plot_genotype`, 
`Population.plot_phenotype`, and `Population.plot_fitness` cannot be run 
for Populations that have not yet been burned in, as they will not yet have
genomes assigned) or the `Population` on which they're called 
(e.g. the aforementioned methods cannot create plots for a `Population` 
that has no `GenomicArchitecture`; and likewise, the 
`Population.plot_demographic_changes` method cannot be called for a 
`Population` for which demographic changes were not parameterized).

The `Land` object and its `Scape`\s also both have a `plot` method.

**********
Parameters
**********

<<<copy of the potential full parameters file, with each parameter 
cross-referenced to a subsequent section explaining its meaning/use>>>


*****************************
Class and function docstrings
*****************************

<<<COPIES OF ALL DOCSTRINGS HERE>>>


