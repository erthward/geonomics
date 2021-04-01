.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash

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

Geonomics is written in Python, a full-fledged scripting language 
that is relatively easy to learn (and fun!). In Python, it can be pretty quick
for a new user to get up to speed and start doing useful work. For work with
Geonomics, this turnaround time should be even quicker. Geonomics aims to
require minimal Python knowledge (yet maintain high extensibility for
interested, expert users). Essentially, anyone should be able to build their
own, arbitrarily complex Geonomics models as long as they know how to install
the package, open a Python console, call Python functions, and edit some
default values in a pre-packaged script. 

