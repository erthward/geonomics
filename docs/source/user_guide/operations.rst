.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


**********
Operations
**********

The following sections discuss the mechanics of core Geonomics operations. 

As you're reading, this conceptual diagram should serve as a useful reference.
Most of the key operations are depicted in the cycle of four boxes
in the corners.

.. image:: ../conceptual_diagram.png



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
