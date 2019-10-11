*********
Geonomics
*********

.. code-block:: python

                                 ##  ::::::          :::    :: ::::::::::##
                           ##:::::    ::::   :::      ::    :: :: ::::::::::: :##
                        ##::::::::     ::            ::   ::::::::::::::::::::::::##
                      ##:::::::::                      :::::::::: :::::: ::::::::  :##
                    ## : ::::  ::                    ::::  : ::    :::::::: : ::  :   ##
                   ##ggggg  eeeee ooooo   nn   nn   ooooo   mm   mm iiiiii  ccccc sssss##
                  ##gg     ee    oo   oo  nnn  nn  oo   oo  mm   mm   ii   cc     ss    ##
                  ##gg     eeee oo     oo nn n nn oo     oo mmm mmm   ii   cc     ssssss##
                  ##gg ggg eeee oo     oo nn  nnn oo     oo mm m mm   ii   cc         ss##
                  ##gg   g ee    oo   oo  nn   nn  oo   oo  mm   mm   ii   cc        sss##
                   ##ggggg  eeeee ooooo   nn   nn   ooooo   mm   mm iiiiii  ccccc sssss##
                    ##    :::::::::             :::::::::::: ::              ::  :   :##
                      ##:   :::::                   ::::::: :::             ::::::: ##
                        ##   :::                     ::::::  ::              :::::##
                           ## ::                      ::::                     ##
                                 ##                                      ##
              
 

**Geonomics** allows users to build and run arbitrarily complex, forward-time,
agent-based, and spatially explicit simulations for landscape genomics. It is
designed to allow novice Python users to create sophisticated simulations with
minimal code, while also allowing advanced users a high level of extensibility
and customizability.

We will continue to expand and add functionality in future versions. Please
contact us with questions, suggestions, or requests!

Main Features
-------------
The following is a short list of highlights. For the full monty, please see the
`original methods paper <PAPER_URL_HERE>`_ and the `homepage <https://github.com/drewhart/geonomics>`_.

    - a model object, which serves as the primary user interface and which
      contains all other model components
    - a landscape object consisting of an arbitrary number of environmental
      raster layers
    - a community object consisting of an arbitrary number of species objects,
      each consisting of an arbitrary number of individuals
    - an optional genomic-architecture object, upon which individuals' genomes
      are based
    - spatialized logistic growth regulating local population densities
    - the capability to model realistic movement and offspring dispersal
      across conductance surfaces
    - neutral and non-neutral evolution capabilities, with spatially contingent
      selection
    - demographic- and environmental-change capabilities
    - the capability to run an arbitrary number of iterations of a model
    - the capability to sample data and a variety of statistics at any desired
      timesteps during a model run
    - numerous visualization methods to aid in model design, exploration,
      analysis, and presentation

Quickstart
----------
The impatient beginner can try the following::

  >>> import geonomics as gnx
  >>> mod = gnx.run_default_model()

This will build and run **geonomics**' default model, return its `Model` object
as `mod`, and leave its parameters file in your current working directory under
the name 'GEONOMICS_default_model_params.py'.


.. note::
Disclaimer
----------
**Geonomics** claims no affiliation with the philosophy and economic ideology
`Georgism <https://en.wikipedia.org/wiki/Georgism>`_, sometimes referred to as
'geonomics'. It is a portmanteau of **geo**\graphy and ge\ **nomics**.
I just thought it sounded neat, and found it delightfully confusing.
