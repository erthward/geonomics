--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------

# CODE CHANGES (ranked in priority from 1 (least) to 5):

        5 should I allow r_distr_alpha to have a value and r_distr_beta to be None, thus fixing all recomb rates at the alpha's value?? this is of course getting far from the alpha/beta meaning coming from the beta distr parameterization, so might want to change their names in that case and then just explain in the docs?

        5 add a 'mating_interval' argument, to define the number of timesteps between mating events? or is this meaningless really?

        5 change 'iterations' to 'runs' across the package and docs

        5 add argument to the mod.plot functions that allows it to take an axes instance and plot there!

        5 for packaging, based on picking through numpy and pandas, it looks as though I need to import all my top-level stuff that I want users to access publicly/easily in my main __init__.py, but then import specific functions/classes where needed within the submodules using "from geonomics.<module> import submodule". See, for example, a comparison of numpy/__init__.py's from . import linalg statement and then the specific import statements in numpy/linalg/linalg.py
                - BUT STILL NOT SURE IF I NEED THE __init__.py FILES IN ALL MY MODULES' SUBDIRECTORIES...
                - looks like now I have a bunch of duplicated, unnecessary, and undesired stuff that is accessible from gnx.<...> once I import (e.g. I don't want gnx.C, gnx.f, gnx.structs; and I think gnx.sys is okay, by comparison with np.sys and pd.np.sys, but I don't want both gnx.sys and gnx.species.sys, nor gnx.np and gnx.species.np and gnx.individual.np, etc...) ---> how import things like numpy in one place at top of package, then get access everywhere else??? have to replace import sys with import genomics.sys in species.py, for example??

        5 only implement the recombination paths thing if recombination rates < 0.5 are present??

        5 fix recombination-rate draws, which throw an 'a <= 0' (i.e. alpha <= 0) error if the user sets r_distr_alpha = 0 in the params

        5 consider changing the ParametersDict class to be something more formally structured, rather than the currently implemented dyanmic-attribute dict, because if users (including myself, having already made this mistake...) try to use gnx.read_parameters_file() to create a params object, then edit it manually in a script before running mod = gnx.make_model(), if they set the dict's attributes using dot notation (i.e. dict.key) then it actually creates a separate attribute than if they set the dict's attributes using classic dict['key'] notation, and this could wind up being a HUGE problem...

        5 MAKE FINAL DECISION ABOUT HOW IMPLEMENT MONO vs POLYGENIC TRAITS! Because if I made 0 the base phenotype for mono and poly, then poly individuals could overshoot z = 1.0 but
          could never undershoot z = 0.0; however, if I leave mono baseline at 0 and poly baseline at 0.5, and then a mono trait undergoes mutation, then from one timestep to the next
          and individ with genotype 1|1 at the originally monogenic trait-locus would go from having phenotype 1.0 before the mutation to 1.5 after the timestep (because ops.selection
          uses traits[n_trait].n_loci to determine how to calculate the phenotype for each trait!
        
        3 maybe add a plot_recombinants method, to get a nice simple visual of recombination (with vertical lines on top of any arbitrary list of loci, if provided as an extra argument)

        5 change polygenic baseline phenotype to 0 also

        5 decide whether to add a data output parameter/option, and a format, for saving the genomic architecture for a trait (to be able to afterward check non-neutral loci, etc)

        5 check whether the Layers' rasters are in fact always checked/coerced to 0 <= e <= 1

        5 change the data-writing function so that each z and e value gets dumped into its own column, by trait and layer name

        5 use the PyTricks email from 18/12/2018 to make all my check-elements-in-a-list-equal calls faster

        5 add parameters to 'defined'-type Layer, so that just a straight-up np.array can be fed in and used as the Layer's raster
                - And add to the docs the explanation that folks, if they want, can import packages and add additional prep code/etc prior to the `params = {` line in the params file!

        5 double-check what the T/F 'sex' param in species.mating is actually doing!

        6 decide whether or not to add a switch determining whether movement surafces are calculated up front, before model is run, or are calculated during the model run (b/c of the memory/compute time ttrade-off)

        5 better burn-in functions for determining stable spatial distribution of individuals across timesteps?

        5 explain that trait mutations are normal distributon (and why) but that deleterious mtuations are gamma (based on the lit)

        5 consider adding a generalized function for making a distribution's draw be a function of phenotype, then apply that to movement, dispersal, fecundity, etc?

        5 consider adding ability for life-history parameters to be determined by a trait's phenotype

        5 touch up all comments (only keep necessary, non-obvious comments; and maybe separate each comment-topped section by a single line break?)

        5 add type annotation to all functions

        5 fix/debug all unit tests

        5 write lots of checks on the types of all the data in the params file, with informative errors

        2 finalize nice string representations for all classes


--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------

# DOCS CHANGES:
- make clearer that "dispersal" means "offspring" and is differentiated from "movement"
- make clear that a Species' K will change because of landscape changes in the Layer it's based on, or because of parameterized changes, OR BOTH (BUT BOTH WILL BE A MESS,
  so user should either include a Layer in the Landscape that won't change, just to start off a Pop's K raster, then tweak that K raster with demo changes, OR plan Landscape
  changes knowing that the Species' K raster will be affected, but NOT BOTH!)
- write docstrings for remaining classes
- add images
- explain dist-parameters' name format somewhere
- put together some examples (from simplest upward), in a Cookbook section of the docs
      - and use one of them as an off-the-bat example that can serve as a vehicle for continuing to discuss all the additional possible parts of a model?
- add to the Intro (or somewhere near there) in the docs a list of the sorts of questions the package could answer???
- rework Intro based on feedback from lab meeting, and on Nick's notes (see his email)
- add mention of "geonomics" the philosophical movememnt (and disambiguation)
- add a graphical representation of the model components to the docs (Species composed of Individuals, optionally having genomes, optionally subtending phenotypes, and living on a Landscape with multiple Layers)


--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------

# FINAL STEPS:

- write popgen-based theoretical validations-tests

- add docstrings to all functions (look up current standards first)

- add an ./ex directory, which will contain a Yosemite simulation

- finalize pkg structure

- prep for PyPi release

- script and run sims for the first study
        - genomic conflict and adaptation to climate change (2 traits, 2 environmental gradients, 1 static, 1 dynamic)
        - LOOK INTO YOSEMITE CHIPMUNKS

- put together analyses and manuscript

- send out for publication
