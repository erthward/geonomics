import geonomics as gnx
import msprime
import tskit
import warnings



# TODO:
    # 0. double-check all aspects of TableCollection in gnx alreay and in
    #    add_individs fn

    # 1. finish drafting add_individs fn

    # 2. debug and test

    # 3. work out function and params-handling for instantiating population
    #    using one or more msprime pops, then debug and test that

    # 4. integrate everything into pkg

    # 5. produce a clean, simple example nb or two (AND TRY TO PUT ON SITE!)

    # 6. clean up repo and publish new version

    # 7. circle back to group with news and script


# OTHER TODO:
    # : DOUBLE CHEKC IM UNDERSTANDING AND USING POPULATIONS TABLE AND
    # nodes['population'] FIELD CORRECTLY!

    # : UPDATE GNX GENERALLY SO THAT IT IS PROPERLY FILLING IN THE
    # provenances INFO IN ANY NEW TABLECOLLECTION!

    # how does the 'location' column in individuals table work again?
    # AND HOW OFTEN IS IT BEING UPDATED?

    # need to do mutations manually? will not work for target allele freqs >0.5

    # how to take allele freq args/parameters?
    # SEEMS PRETTY REASONABLE THAT WE SHOULD ALLOW DIFF PARAMS FOR EACH MSPRIME
    # POP!

    # figure out how to alter node and edge info so that there are no problems
    # when fusing into main TableCollection

    # review species._set_genomes_and_tables and genome.make_starting_mutations
    # once more, and figure if/how I should best integrate this into them
    # (e.g., _set_genomes_and_tables can check for t=0 components in the 'init'
    # and call that general-purpose function within its body if so

    # add spot in 'init' params to allow quick creation of a model using 
    # arbitrary msprime demog! then add to params docco too

    # add the extra functions and start-up code needed to handle those params


def _handle_msprime_rate(rate, L):
    """
    Either return a simple float, if the given rate
    is homogeneous, or else return an msprime.RateMap object.
    Used for msprime recombination and mutation rates.
    """
    assert type(rate) in [float, list, numpy.ndarray]
    if (isinstance(rate, float):
        return rate
    elif len(np.unique(recomb_rate)) == 1):
        return rate[0]
    else:
        assert len(rate) == L
        # NOTE: position array has to be one longer than the sequence length,
        #       so that rightmost position brackets a stretch of genomic
        #       intervals equal to the sequence_length
        positions = np.linspace(0, L, L+1)
        ratemap = msprime.RateMap(position=positions, rate=rates)
        return ratemap


def sim_msprime_individs(n,
                         L,
                         recomb_rate,
                         mut_rate,
                         demography=None,
                         population_size=None,
                         gnx_spp_use_tskit=True,
                         ancestry_model=None,
                         random_seed=None,
                        ):
    '''
    Use an msprime simulation to generate new Geonomics Individuals.

    Parameters
    ----------

    n: int
        The number of Geonomics Individuals to be simulated

    L: int
        The length of the simulated genome.

    recomb_rate: {float, list, numpy.ndarray}
        Either a scalar float, indicating the genome-wide recombination rate;
        or a 1d numpy.ndarray or 1d-array-like of length L-1, indicating the
        inter-locus recombination rates along the length-L simulated genome.
        (Note: This parameter is expressed according to the convention
        used by Geonomics, as explained above. This differs from the convention
        used by msprime, the translation to which will be handled internally by
        Geonomics.)

    mut_rate: {float, list, numpy.ndarray}
        Either a scalar float, indicating the genome-wide neutral mutation rate;
        or a 1d numpy.ndarray or 1d-array-like of length L, indicating the
        mutation rates of all loci in the the length-L simulated genome.
        (Note: This parameter is expressed according to the convention
        used by Geonomics, as explained above. This differs from the convention
        used by msprime, the translation to which will be handled internally by
        Geonomics.)

    demography: msprime.Demography, optional, default: None
        An `msprime.Demography` object to be used by `msprime.sim_ancestry`.
        `msprime` allows specification of complex demographic histories,
        including multiple extant and ancestral populations with fixed
        or variable sizes, growth rates, or migration rates, as well as
        timed demographic events. Defaults to None, in which case
        `msprime.sim_ancestry` will simulate a single population
        whose size will be determined by the `population_size` parameter
        (or will default to 1, if `population_size` is also None).
        Cannot be specified along with `population_size`.
        (See the `msprime.Demography` documentation and references therein
        for more details.)

    population_size: int, optional, default: None
        The population size of the default single population to be simulated
        by `msprime.sim_ancestry`. That default single population will only be
        simulated if a more complex demographic scenario is not specified using
        the `demography` parameter (i.e., if `demography` is None).
        Thus, `population_size` cannot be specified along with `demography`.
        If both are None then `msprime` simulates a single population of size 1.
        (See the `msprime.sim_ancestry` documentation for more details.)

    gnx_spp_use_tskit: bool, optional, default: True
        Whether or not the Geonomics `Species` for which the `Individual`s are
        intended uses tskit (i.e., has Species.gen_arch.use_tskit as True).
        If not, `Individual`s will each store a copy of their Lx2 genotype
        array in the 'g' attribute; otherwise, this attribute will be None.

    ancestry_model: {`msprime.AncestryModel`, list}, optional, default: None
        The `msprime.AncestryModel instance (or list of them) to be used by
        `msprime.sim_ancestry`. See `msprime.sim_ancestry` documentation and
        references therein for more details.

    random_seed: int, optional, default: None
        The random seed value to be provided to `msprime.sim_ancestry` and
        `msprime.sim_mutations`. Use if reproducibility is desirable.
        See `msprime.sim_ancestry` and `msprime.sim_mutations` documentation
        for detail.

    Returns
    -------
    :class:`list`
        A list of `geonomics.structs.individual.Individual` objects.
        Individuals are assigned sequential index integers starting from 0,
        placed at x,y coordinates (0, 0), and otherwise given default attribute
        values. Attributes should be updated as necessary in downstream code.

    '''
    # handle the recombination and mutation rates
    if not isinstance(recomb_rate, float):
        # NOTE: tack a 0.5 onto the righthand side of the recomb_rate vector,
        #       to satisfy the difference between msprime (for which L
        #       expresses length in intervals, each of which needs a recomb
        #       rate specified, hence the recomb_rate vector needing to have a
        #       number of rates equal to L) and gnx (for which L expressed
        #       length in loci, between which recomb rates are needed, hence
        #       needing L-1 recomb rates);
        #       TODO: DOUBLE-CHECK THE FOLLOWING LOGIC
        #       this 0.5 value will not affect gnx data because it effectively
        #       sits between sites L-1 and L, but a simulated gnx genome of
        #       length L only includes sites [0, L-1]
        recomb_rate = [*recomb_rate] + [0.5]
    recomb_rate = _handle_msprime_rate(recomb_rate, L)
    mut_rate = _handle_msprime_rate(mut_rate, L)

    # simulate ancestry and mutations
    # TODO: look back at how gnx uses genome positions, because 
    #       it looks like we need to provide discrete_genome=false
    #       to sim_ancestry and true to sim_mutations;
    #       or perhaps instead just fine to use false for both
    #       and then subtract 0.5 from all non-zero 'left' values
    #       from all non-L 'right' values when handles the edges table?
    ts = msprime.sim_ancestry(samples=n,
                              sequence_length=L,
                              demography=demography,
                              population_size=population_size,
                              recombination_rate=recomb_rate,
                              ploidy=2,
                              model=ancestry_model,
                              random_seed=random_seed,
                             )
    # NOTE: mutations use a simple binary 0|1 model, with 0 always the
    #       ancestral state, and with mutations being state-independent
    #       (i.e., always placing a 1, even if this results in a silent
    #       mutation of 1→1);
    #       this model is the best aligned to gnx data
    mut_model = msprime.BinaryMutationModel(state_independent=True)
    mts = msprime.sim_mutations(tree_sequence=ts,
                                rate=mut_rate,
                                model=mut_model,
                               )
    # get the mutated tskit.TableCollection object
    mtc = mts.dump_tables()

    # double-check the number of Individuals
    assert np.sum(mtc.nodes.asdict()['individual']>-1) == 2*n

    # create the Individuals, assigning placeholder values for the attributes
    # (to be later updated by downstream code)
    individs = {i: gnx.structs.individual.Individual(idx=i,
                                                     x=0,
                                                     y=0) for i in range(n)}

    # get the information linking Individuals to rows in TableCollection tables
    nodes_tab_individs_col = mtc.nodes.asdict()['individual']
    for i, ind in individs.items():
        # get the nodes table rows corresponding to this Individual
        nodes_tab_ids = np.where(nodes_tab_individs_col==i)[0]
        assert len(nods_tab_ids) == 2
        # set hidden attributes tracking the individs' rows in the 'nodes' table
        # and id in the 'individuals' table
        ind._nodes_tab_ids = dict(zip(range(2), nodes_tab_ids))
        ind._individuals_tab_id = i
        # if the gnx Model these Individuals are intended for is not using tskit
        # then walk the variants generator and compile the individ's
        # genotype array, then store it in the attribute 'g'
        if not gnx_spp_use_tskit:
            var_gen = mts.variants(samples=nodes_tab_ids)
            gt = np.stack([next(var_gen).genotypes for _ in range(L)])
            ind.g = gt

    return individs, mtc



##############
# Model method
def add_individs(n,
                 coords,
                 recip_spp=0,
                 source_spp=None,
                 source_msprime_params=None,
                 individs=None,
                ):
    """
    Add individuals to a recipient Species object, either using a second
    Species object as the source population or feeding a dict of parameters
    into msprime to simulate the source population.

    Parameters
    ----------
    n: int
        The number of Individuals to be added.

    coords: {tuple, list, numpy.ndarray}
        A tuple, list of numpy.ndarray indicating the x,y coordinates where the
        new Individuals should be added to the recipient Species.
        If `coords` contains 2 values then all added Individuals will be placed
        at the x,y coordinate pair indicated by those 2 values.
        Otherwise, `coords` may be of size `n` x 2, indicating the x,y
        coordinate pair at which each of the `n` Individuals will be placed.

    recip_spp: {gnx.structs.species.Species, int, str}, optional, default: None
        An int or str referring to the Species in this Model to which
        Individuals should be added.

    source_spp: {gnx.structs.species.Species, int, str}, optional, default: None
        A Geonomics.structs.species.Species object to act as the source
        population from which individuals will be taken
        (or an int or str referring to that Species within this Model).

    source_msprime_params: dict, optional, default: None
        A dict of keyword args to be fed to `gnx.sim_msprime_individs()`,
        to parameterize the msprime source population model from which
        new Individuals will be drawn. These include:
            recomb_rate         (required)
            mut_rate            (required)
            demography          (optional; default: None)
            population_size     (optional; default: None)
            ancestry_model      (optional; default: None)
            random_seed         (optional; default: None)
        See `gnx.sim_msprime_individs` for details.

    individs: {tuple, list, numpy.ndarray}, optional, default: None
        A 1d numpy.ndarray or equivalent array-like containing the int indices
        of the Individuals to be taken from `source_spp` and added to
        `recip_spp`. If None

    Returns
    -------
    None
        Alters `recip_spp` in place by adding the provided
        (if `source_spp` is not None) or simulated (if `source_msprime_params`
        is not None) Individuals.

    """

    warnings.warn(("This function is new and has been checked but not "
                   "extensively tested. User beware. Please immediately report "
                   "any suspected bugs at "
                   "https://github.com/erthward/geonomics/issues "))

    # get the recipient species
    recip_spp = self.comm[self._get_spp_num(recip_spp)]

##############
#Species method
    # make sure either a Species object is given or a set of valid msprime args
    assert ((source_spp is not None and source_msprime_params is None) or
            (source_spp is None and source_msprime_params is not None)), (
                    "Source population must be provided either as "
                    "a gnx Species object ('source_spp') "
                    "or as a set of valid args for gnx.run_msprime_sim() "
                    "('source_msprime_params'), but not both.")
    if source_spp is not None:
        assert (isinstance(source_spp, gnx.structs.species.Species) or
                isinstance(source_spp, int) or
                isinstnace(source_spp, str)), ("The "
                    "value given to 'source_spp' identifies a gnx "
                    "Species object that will serve as the source "
                    "population, so must be either a gnx Species "
                    "instance (i.e., an object of the class "
                    "gnx.structs.species.Species) whose Individuals will "
                    "be added to the recipient Species, or a valid int or str "
                    "identifier of another Species in this model whose"
                    "Individuals will be added to the recipeint Species.")

    # either get Individuals and their tskit.TableCollection
    # from the given gnx Species object...
    if source_spp is not None:
        # raise warning about potential unforeseen effects of introducing
        # Individuals from an incorrectly parameterized source Species
        warnings.warn(("Introducing Individuals from one Species object "
                       "into another could cause unflagged issues. "
                       "Geonomics will check and prevent conflict between strictly "
                       "incompatible attributes of the Species, "
                       "GenomicArchitectures, and Traits. However, attributes "
                       "whose compatibility depends on the scientific "
                       "reasoning behind the scenario being simulated "
                       "(e.g., mutational parameters, "
                       "carrying capacity parameters, etc.) "
                       "will be allowed to differ between the Species and, "
                       "if these differences were unintended but go "
                       "undetected, could lead to incorrect or misleading data "
                       "and results. Please use care when planning, "
                       "implementing, and checking your code."))

        # get the source Species from within this model, if necessary
        if isinstance(source_spp, int) or isinstance(source_spp, str):
            source_spp = self.comm[self._get_spp_num(source_spp)]

        # assert that either n or individs is provided, not both,
        # and that size of n or size and contents of individs are valid
        assert ((n is not None and individs is None) or
                (n is None and individs is not None)), ("If using another "
                                                        "Geonomics Species "
                                                        "object as the source "
                                                        "population then "
                                                        "either `n` or "
                                                        "`individs` must be "
                                                        "provided, "
                                                        "but not both.")
        if individs is None:
            assert n <= len(source_spp), ("`n` must not be greater than "
                                          "the population size of "
                                          "`source_spp`.")
        else:
            assert (len(np.unique(individs)) <= len(source_spp) and
                    np.all([ind in source_spp for ind in individs]), (
                        "`individs` contains Individual indices that do "
                        "not exist in `source_spp`.")

        # check that all Species params that need to match are the same
        # in the source population's Species and the recipient Species
        spp_attrs_to_check = ['K_layer',
                              'selection',
                              'sex_ratio',
                              'move',
                             ]
        for attr in spp_attrs_to_check:
            assert np.all(getattr(source_spp, attr) ==
                          getattr(recip_spp, attr)), (
                              "both the source and recipient Species "
                              "must have identical values for "
                              f"the attribute '{attr}'")

        # check that all gen_arch params that need to match are the same
        # in the source population's Species and the recipient Species
        source_gen_arch = source_spp.gen_arch
        recip_gen_arch = recip_spp.gen_arch
        assert source_gen_arch.L == recip_gen_arch.L, ("source "
                                           "Species must have a "
                                           "simulated genome of the same "
                                           "length as the recpient Species.")
        gen_arch_attrs_to_check = ['sex',
                                   'use_tskit',
                                   'x',
                                  ]
        # NOTE: not requiring p to be the same in both populations facilitates
        #       creating loci that differ in their standing genetic variation in
        #       the two populations, or even loci that are fixed differently
        #       in one or more populations
        # NOTE: making it so that attrs like mutation rate, dominance params,
        #       and the arrays of neutral and deleterious loci needn't be the
        #       same (mutated loci almost certainly wouldn't be
        #       if non-neutral mutation has been happening independently in
        #       both populations). it isn't entirely realistic to model this
        #       this way, but also not entirely unrealistic for the model to
        #       behave such that loci that were neutral in a population's prior
        #       environment are suddenly not neutral in the new environment
        for attr in gen_arch_attrs_to_check:
            assert np.all(getattr(source_gen_arch, attr) ==
                          getattr(recip_gen_arch, attr)), ("both the source "
                                                     "and recipient Species' "
                                                     "GenomicArchitectures "
                                                     "objects must have "
                                                     "identical values for the "
                                                     f"attribute '{attr}'")
        # check recombination attributes
        recomb_attrs_to_check = ['_rates',
                                 '_r_distr_alpha',
                                 '_r_distr_beta',
                                 '_jitter_breakpoints',
                                ]
        for attr in recomb_attrs_to_check:
            assert np.all(getattr(source_gen_arch.recombinations, attr) ==
                          getattr(recip_gen_arch.recombinations, attr)), (
                                                     "both the source "
                                                     "and recipient Species' "
                                                     "Recombinations objects "
                                                     "must have identical "
                                                     "values for the "
                                                     f"attribute '{attr}'")

        # check traits and their params
        if recip_gen_arch.traits is None:
            assert source_gen_arch.traits is None, ("Individuals "
                                              "from a Species with "
                                              "Traits cannot be added to a "
                                              "Species without them.")
        else:
            assert len(source_gen_arch.traits) == len(recip_gen_arch.traits), (
                "The source population's Species object must have the "
                "same number of Traits as the recipient Species.")
            trt_attrs_to_check = ['name',
                                  'phi',
                                  'lyr_num',
                                  'max_alpha_mag',
                                  'gamma',
                                  'univ_adv',
                                 ]
            # NOTE: not requiring 'n_loci' to be the same because they're
            # unlikely to be the same if non-neutral mutation has been happening
            # independently in both populations
            for n, trt in recip_gen_arch.traits.items():
                for attr in trt_atts_to_check:
                    assert np.all(getattr(trt, attr) ==
                                  getattr(source_gen_arch.traits[n], attr)), (
                                        f"Attribute '{attr}' of trait {n} "
                                         "in the source population's Species "
                                         "object must be the same as "
                                         "that attribute in the recipient "
                                          "population's Species.")

        # get individuals from the given Species object
        if individs is not None:
            individs = [ind for idx, ind in source_spp.items() if idx in individs]
        else:
            individs = [*source_spp.values()]

        # get their tskit.TableCollection
        source_tc = source_spp._tc

    # ... or simulate them w/ the given msprime args
    # and the recipient Species' genomic architecture
    elif source_spp is None:
        assert individs is None, ("The 'individs' argument indicates the "
                                  "indices of Individuals to be added from "
                                  "a gnx Species object serving as a source "
                                  "population, so can only be used if "
                                  "'source_msprime_params' is None and "
                                  "'source_spp' is provided.")

        # make sure that Species has no selection before allowing use of msprime
        assert not recip_spp.selection, ("msprime is a coalescent simulator, "
                                         "so cannot be used to simulate "
                                         "selection. Consider using a second "
                                         "Geonomics model to simulate a source "
                                         "population with selection.")

        # check that required msprime params are provided and any others
        # provided are among those allowed
        required_msprime_params = ['recomb_rate', 'mut_rate']
        allowed_msprime_params = ['demography', 'population_size',
                                  'ancestry_model', 'random_seed',
                                 ]
        all_msprime_params = required_msprime_params + allowed_msprime_params
        assert np.all([k in all_msprime_params for
                        k in source_msprime_params]), ("The only parameters "
                                "accepted in 'source_msprime_params' are: "
                               f"{', '.join(all_msprime_params)}.")
        assert np.all([k in source_msprime_params for
                        k in required_msprime_params]), ("Both "
                               f"{' and '.join(required_msprime_params)} "
                                "must be provided within "
                                "'source_msprime_params'.")

        # run the simulation and return template Individuals
        use_tskit = recip_spp.gen_arch.use_tskit
        individs, source_tc = sim_msprime_individs(n=n,
                                                   L=recip_spp.gen_arch.L,
                                                   gnx_spp_use_tskit=use_tskit,
                                                   **source_msprime_params,
                                                   )

    # reset/set all source Individuals' x,y coordinates
    assert np.atleast_2d(coords).shape in [(1,2), (len(individs),2)], (
        "The 'coords' arg must be provided either "
        "a single x,y coordinate pair (at which all Individuals "
        "will be introduced) or an nx2 array of x,y coordinate pairs "
        "(indicating the introduction locations for "
        "each of the n Individuals being introduced.")
    if np.atleast_2d(coords).shape == (1, 2):
        coords = np.stack([coords for _ in range(len(individs))])
        assert coords.shape == (len(individs), 2)
    [ind._set_pos(*coords[i,:]) for i, ind in enumerate(individs)]


    # reset source Individuals' indices
    # and increment the recipient Species's max idx concordantly
    max_idx_b4 = recip_spp.max_ind_idx
    for ind in individs:
        recip_spp.max_ind_idx += 1
        ind.idx = recip_spp.max_ind_idx
    max_idx_af = recip_spp.max_ind_idx
    assert max_idx_af - max_idx_b4 == len(individs)

    # TODO: recalculate Individuals' node IDs, individual IDs, anything else
    #       necessary, by adding spp._tc.nodes.num_rows,
    #       and update this information within the Individuals' attributes
    #       LOOK AT SIMPLE EXAMPLE IN MY SCRIPT UNIONING 2 SMALL TCs!

    # sort and simplify the recipient Species' tskit.TableCollection
    recip_spp._sort_and_simplify_table_collection()

    # TODO: make sure each source Individual's introduction time is flagged
    #       somewhere within the TableCollection

    # union source individuals' TableCollection into that of recipient Species
    recip_spp._tc.union(source_tc,
                        # NOTE: (for now at least) assuming that the introduced
                        #       individuals come from totally independent pop,
                        #       do not share a MRCA within our data, and thus
                        #       have no common nodes that need to be mapped,
                        #       hence we provide a list of tskit.NULLs
                        #       to union()'s node_mapping argument
                        node_mapping=[tskit.NULL] * len(individs),
                        # NOTE: populations of the newly added nodes
                        #       will be flagged by an integer
                        #       that is 1 greater than the largest
                        #       population-flag integer in the pre-introduction
                        #       recipient Species' TableCollection.nodes table
                        add_populations=True,
                        # NOTE: best to just allow it to do this,
                        #       in case we generalize things later on,
                        #       but for now there will actually be nothing to
                        #       check because we assume no shared nodes between
                        #       source and recipient TableCollections
                        check_shared_equality=True,
                        # NOTE: should copy provenance info over from the source
                              # Individuals' TableCollection to that of the
                              # receipient Species
                        record_provenance=True,
                             )

    pop_size_b4 = len(recip_spp)
    # TODO: add the individs to the recipient Species, keyed by their idx values

    pop_size_af = len(recip_spp)
    assert (pop_size_af - pop_size_b4) == len(individs), ("Incorrect number "
                            "of Individuals added to recipient Species.")

    # TODO: reset/set source Individuals' environmental and fitness values

    # TODO: CHECK FOR ANY OTHER INDIVID- OR SPP-LEVEL ATTRS NEEDING RESET!

    # TODO: reset the Species' coords and cells attrs

    # TODO: reset the KDTree after adding the new individuals

    # print informative output
    print(f"\n{len(individs)} Individuals successfully removed.\n")






source_pops = {'A': {'initial_size':10000,
                     'growth_rate':0,
                     'description':'first'},
               'B': {'initial_size':1000,
                     'growth_rate':0.1,
                     'description':'second'},
               'C': {'initial_size':100,
                     'growth_rate': -0.1,
                     'description':'last',
                    },
              }

source_pop_intros = {'A': [25, (4, 6)],
                     'B': [50, (1, 3)],
                     'C': [100, (5,6)],
                    }

source_pop_allele_freqs = {'A': [0.1]*100,
                           'B': [0.9]*100,
                           'C': [0]*100,
                          }


def _make_msprime_demographies(source_pops):
    """
    use all of the distinct source populations identified in the provided
    parameters dict to create msprime.Demography objects that
    can be used to draw new gnx individuals
    """
    demogs = []
    for name, params in source_pops.items():
        disallowed_params = ['default_sampling_time',
                             'initially_active',
                            ]
        # NOTE: need to disallow certain parameters that might break gnx models
        #       render output inaccurate
        for dp in disallowed_params:
            assert dp not in params, ((f"parameter '{dp}' is not allowed in "
                                        "msprime simulations used to generate "
                                        "geonomics models; please reformulate "
                                        "population parameters, then rerun."))
        demog = msprime.Demography()
        demog.add_population(name=name, **params)
        demogs.append(demog)
    return demogs


def _get_msprime_demographies(source_pops):
    """
    either grab the msprime.Demography object
    (or list of msprime.Demography objects) provided in the gnx parameters file
    (for more complicated models involving e.g., migrations, splits, etc.),
    or grab the demographic parameters dict provided in the gnx parameters file
    and use it to create a list of msprime.Demography objects
    """
    if (isinstance(source_pops, msprime.Demography) and
        source_pops.num_populations > 0):
        return [source_pops]
    elif (isinstance(source_pops, list) and
          set([type(sp) for sp in source_pops]) == set([msprime.Demography])):
        return source_pops
    elif isinstance(source_pops, dict):
        return _make_msprime_demographies(source_pops)
    else:
        raise ValueError((f"'source_pops' component of species 'init' "
                           "parameters must be either an msprime.Demography "
                           "with at least one population or a dictionary of "
                           "name-keyed demographic parameters that can be used "
                           "to instantiate msprime.Population objects for a "
                           "Demography object."))


def _add_individs_from_msprime_source_pops(n_individs,
                                           demogs,
                                           spp,
                                           mut_rates=None
                                          ):
    """
    simulate the numbers of individuals indicated by n_individs
    (either an int, if demogs has only 1 population,
    or else a population-name-keyed dict of ints)
    using the given msprime.Demography and gnx.GenomicArchitecture objects,
    and adding then to the given tskit.TableCollection
    """
    assert (len(spp.gen_arch.nonneut_loci) == 0 and
            spp.gen_arch.traits is None), ("msprime can only currently "
                                           "be used to seed a model with "
                                           "individuals from source "
                                           "populations if the model "
                                           "has no selection. "
                                           "This constraint may be "
                                           "relaxed in the future.")
    assert spp.gen_arch.use_tskit, ("to use msprime to seed a population with "
                                    "individuals from source populations, "
                                    "please set the 'use_tskit' parameter to "
                                    "True".)
    assert ((isinstance(n_individs, int) and len(demogs) == 1) or
            (isinstance(n_individs, dict) and
                (len(demogs) == 1 and
                 isinstance(demogs[0], msprime.Demography) and
                 demogs[0].num_populations == len(n_individs)) or
                (len(n_individs) == len(demogs)))):
    if mut_rates is not None:
        assert len(mut_rates) == len(demogs)
    # get the tree-seq object for the simulated ancestry
    if len(np.unique(spp.gen_arch.recombinations._rates[1:])) == 1:
        recomb_rate = spp.gen_arch.recombinations._rates[1]
    else:
        # NOTE: position array has to be one longer than the sequence length,
        #       so that rightmost position brackets a stretch of genomic
        #       intervals equal to the sequence_length
        recomb_rate= msprime.RateMap(position=np.linspace(0,
                                                          spp.gen_arch.L,
                                                          spp.gen_arch.L+1,
                                                         ),
                              rate=[*spp.gen_arch.recombinations._rates[1:]]+[0.5])
    if len(demogs) == 1 and isinstance(demogs[0], msprime.Demography):
        ts = msprime.sim_ancestry(samples=n_individs,
                                  demography=demogs[0],
                                  sequence_length=spp.gen_arch.L,
                                  recombination_rate=recomb_rate,
                                  ploidy=2,
                                  # NOTE: DiscreteTimeWrightFisher ('dtwf') may be
                                  #       closest match to the discrete operations
                                  #       of gnx's forward-time simulations, but
                                  #       probably not necessary to match that and
                                  #       probably simpler, faster, and most
                                  #       reasonably expected to just use the
                                  #       standard coalescent
                                  model='hudson',
                                 )
        tcs = [ts.dump_tables()]
    else:
        tcs = []
        i = 0
        for n, demog in zip(n_individs, demogs):
            ts = msprime.sim_ancestry(samples=n,
                                      demography=demog,
                                      sequence_length=spp.gen_arch.L,
                                      recombination_rate=recomb_rate,
                                      ploidy=2,
                                      model='hudson',
                                     )
            # add mutations according to the gen_arch
            # (using starting allele frequencies as the mutation rate,
            # so that mutations will wind up approximating them)
           if mut_rates is None and len(np.unique(spp.gen_arch.p)) == 1:
                mut_rate = spp.gen_arch.p[0]
                # don't bother with mutations if mut_rate == 0
                if mut_rate == 0:
                    break
           else:
               if mut_rates is not None:
                   mut_rate_vals = mut_rates[i]
                else:
                    mut_rate_vals = spp.gen_arch.p
                # NOTE: position array has to be one longer than the sequence length,
                #       so that rightmost position brackets a stretch of genomic
                #       intervals equal to the sequence_length
                mut_rate= msprime.RateMap(position=np.linspace(0,
                                                          spp.gen_arch.L,
                                                          spp.gen_arch.L+1,
                                                         ),
                                          rate=mut_rate_vals,
                                         )
            mts = msprime.sim_mutations(tree_sequence=ts,
                                        rate=mut_rate,
                                        # state-independent binary mutation
                                        # means that all initial mutations
                                        # will create 1s, exactly as we want
                                        # when we want to create
                                        # specific initial allele frequencies
                                        model=msprime.BinaryMutationModel(state_independent=True),
                                       )
            tcs.append(mts.dump_tables())
            i += 1

    # TODO: need some assert statements that can be run based on comparing
    #       calculations using spp._tc before and after the additions!

    # for each TableCollection
    # (i.e., each introduction event from a source population)
    for tc in tcs:
        # increment all node and edge IDs to be > max node and edge IDs in the
        # species' current TableCollection

        # add new individuals to the TableCollection, matching parent nodes, etc

        # add new individuals to the species, using correct TableCollection
        # ids, correction geolocation, etc

        # sort and simplify spp._tc

    return



