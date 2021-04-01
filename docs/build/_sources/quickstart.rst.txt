.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash


Getting started
***************

Quick start
-----------

**For the beginner**, we recommend the following steps:
  1. Review the diagram of a typical workflow, below.
  2. Read the 'Overview' subsection (below), and the following two sections
     ('Data structures' and 'Operations'), to get a general
     undertsanding of the logic, components, and necessary and optional
     behaviors of a Geonomics model.
  3. Skim the 'Parameters Guide', to understand the structure
     and use of a Geonomics parameters file.
  4. Use `pip` to install Geonomics (:bash:`$ pip install geonomics`);
  5. Open Python and run :py:`import geonomics as gnx`;
  6. Use the `gnx.make_parameters_file` function, to begin
     creating template parameters files that they can tweak as desired;
  7. Use the :py:`gnx.make_model` function and then the
     :py:`Model.walk` or :py:`Model.run` methods to instantiate and run
     the models they've parameterized;
  8. Use the various :py:`Model.plot` methods to visualize
     the behavior and results of their models.

As you work through these steps, it may be helpful to refer to the
:ref:`Workflow Diagram` for guidance.


*Really* Quick Start
--------------------

**For the impatient beginner**:

1. Install Geonomics
2. Launch Python
3. Import Geonomics

.. code-block:: python

     >>> import geonomics as gnx

4. Run the default model, and start mussin' around!

.. code-block:: python

     >>> mod = gnx.run_default_model()

This will load the Geonomics package as `gnx`, create a default Geonomics
parameters file in your current working directory, then use that file
to instantiate and run a :code:`Model` using the default parameter values.
