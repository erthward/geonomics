.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


********
Overview
********

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
    parameters describing what it's like and how it lives;
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

***************
Getting started
***************

**For the beginner**, we recommend the following steps:
  1. Review the diagram of a typical workflow, below.
  2. Read the following two sections ('Data structures',
     'Operations'), to get a general
     undertsanding of the logic, components, and necessary and optional
     behaviors of a Geonomics model.
  3. Skim the subsequent section ('Parameters'), to understand the structure
     of a Geonomics parameters file.
  4. Use `pip` to install Geonomics (:bash:`$ pip install geonomics`);
  5. Open Python and run :py:`import geonomics as gnx`;
  6. Use the :func:`gnx.make_parameters_file` function, to begin
     creating template parameters files that they can tweak as desired;
  7. Use the :py:`gnx.make_model` function and then the
     :py:`Model.walk` or :py:`Model.run` methods to instantiate and run
     the models they've parameterized;
  8. Use the various :py:`Model.plot` methods to visualize
     the behavior and results of their models.

**For the `impatient` beginner**, once you have installed Geonomics
you should be able to run the following code:

.. code-block:: python

     import geonomics as gnx

     gnx.run_default_model()

This will load the Geonomics package, create a default Geonomics
parameters file in your current working directory, 
then use that file to instantiate and run a :code:`Model` using the default
parameter values.

**Typical workflow:**

.. image:: ../procedural_diagram.png


----------------------------------------------------

**********
Motivation
**********

Backward-time (i.e. coalescent) simulators abound.
But they are inadequate for simulation of many scenarios of 
interest, including: natural selection on traits with arbitrary genomic 
architectures; spatially variable natural selection; simulation of species or
populations distributed continuously and moving realistically across
complex landscapes; complex demographic change simultaneous with ongoing, 
often non-stationary environmental change; and coevolutionary interactions 
between multiple species or incipient species. Few existing forward-time 
simulators can model all of these phenomena, and those that can are incredibly
powerful, but often impose a high cost of entry. Geonomics aims to fill 
this empty niche by combining ease of use with broad extensibility. 
If it succeeds at doing this, Geonomics should prove uniquely useful
for a wide range of purposes, from intro-level educational use to
high-quality theoretical, methodological, empirical, and
applied research.


