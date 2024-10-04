.. role:: py(code)
   :language: python

.. role:: bash(code)
   :language: bash

.. role:: underline
   :class: underline


.. code-block:: python
         
       .                  .   ::::::          :::    :: :::::::::::.           .
                    .::::::    ::::   :::      ::    :: :: ::::::::::: ::.                 .
    .      .     .:::::::::     ::            ::   :::::::::::::::::::::::::.        .   
               .::::::::::                      :::::::::: :::::: ::::::::  ::.
         .   .  : ::::  ::                    ::::  : ::    :::::::: : ::  :    .      .
            . 66666 :6666: 00000   66   66   66666   00   00 666666  66666 00000 .
           . 66     66    00   00  666  66  66   66  00   00   66   66     00     .         .
 .         . 66     66   00     00 66 6 66 66     66 000 000   66   66     000000 .
           . 66 666 6666 00     00 66  666 66     66 00 0 00   66   66         00 .
           . 66   6 66    00   00  66   66  66   66  00   00   66   66        000 .     .
            . 66666 :6666: 00000   66   66   66666   00   00 666666  66666 00000 .         .
      .      .    : ::::::::               :::::::::: ::              ::  :   : .
     .         .:    :::::                    :::::: :::             :::::::  .      .
  .              .    :::                      :::::  ::              ::::: .              . 
           .        .  ::                      ::::                      .
                          . ::                                     .                .
          

Geonomics |version|
===================

         

Description
***********

A Python package for simulation of genomic evolution
on complex and dynamic landscapes.
Provides for easy construction of individual-based, spatially explicit,
forward-time simulation models under arbitrarily complex scenarios.

Using minimal code, build models with arbitrarily complex scenarios,
including spatially varying selection, selection on multiple, monogenic or
polygenic traits, and non-stationary demographic and environmental change.


Key Features
************

1. An object-oriented scripting framework, allowing for easy model construction, customization, and extension

2. Model set-up from a single, well annotated parameters file

3. Tools for customizable collection of data throughout a simulation

4. Ability to model complex evolutionary scenarios, including custom demographic change, spatially varying selection, and multiple polygenic traits 

5. Ability to model complex spatial scenarios, including multi-layer simulated or real-world landscapes, resistance-based movement, and non-stationary environmental change

6. Optional reportng of the ancestral recombination graph (ARG) and spatial pedigree using `tskit <https://tskit.dev/tskit/docs/stable/introduction.html>`

7. Numerous visualization tools, to help the user design models and explore results


About the Documentation
***********************

This documentation is designed to be read from from the top down,
as information becomes increasingly detailed.

To jump right in, check out :ref:`Getting started` and :ref:`Examples`.

For more information, see :ref:`Motivation`, :ref:`Data structures`,
and :ref:`Operations`.

For fine detail about particular :ref:`Parameters`,
:ref:`Data structures`, or :ref:`Operations`, see those sections.

**Merry modeling!**


.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   Installation <install>
   Quick Start <quickstart>
   Intro Video <intro_video>
   Workflow Diagram <procedural_diagram>
   Motivation <motivation>
   Examples <examples>


.. toctree::
   :maxdepth: 2
   :caption: User Guide

   Conceptual Diagram <conceptual_diagram>
   Overview <overview>
   Visualization <viz>


.. toctree::
   :maxdepth: 1
   :caption: Reference Guide

   Parameters <params>
   API <modules>


------------------------------------------------------------------

Attribution
...........

This package was written by Drew Ellison Hart, as part of his PhD work.
It is available to freely distribute and modify, with proper
attribution, under the MIT License.

Should you use Geonomics for research, education, or other purposes,
please cite it as:

    Terasaki Hart, D.E., Bishop, A.P., Wang, I.J. 2021. Geonomics: forward-time, spatially explicit, and arbitrarily complex landscape genomic simulations. Manuscript submitted for publication.

Should you have any questons or
concerns, please feel free to get in touch at drew.hart<at>berkeley<dot>edu !

------------------------------------------------------------------

Disclaimer
..........

**Geonomics** claims no affiliation with the philosophy and economic ideology
`Georgism <https://en.wikipedia.org/wiki/Georgism>`_, sometimes referred to as
'geonomics'. It is a portmanteau of
:underline:`geo`\graphy and ge\ :underline:`nomics`.

We just thought it sounded neat, and found it delightfully confusing.


Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

