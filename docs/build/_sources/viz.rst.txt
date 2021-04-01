.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


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

**UNDER CONSTRUCTION!**

|construction|

.. |construction| image:: construction.gif
   :align: middle
   :width: 500
   :height: 500
