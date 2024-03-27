import geonomics as gnx
import msprime
import tskit


"""
TODO:
    1. finish up working through all aspects of the add_individs code, except
       tskit pieces of code
    2. clarify and firm up all tskit bits (USE CODE AT BOTT!)
    3. integrate into package
    4. test and debug (including all arg variants!)
    5. produce simple example script/notebook
    6. clean up repo
    7. push to new package release and version
    8. circle back to group!
"""

##############
# Model method
def add_individs(coords,
                 recip_spp=0,
                 source_spp=None,
                 source_msprime_params=None,
                 individs=None,
                 # TODO: ARE THESE LAST PARAMS USEFUL?
                 pop_info=None,
                 tskit_provenance_info=None,
                ):
    # TODO:
    """
    DOCCO!
    """

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
        # get the source Species from within this model, if necessary
        if isinstance(source_spp, int) or isinstance(source_spp, str):
            source_spp = self.comm[self._get_spp_num(source_spp)]

        # TODO: check that all genomic and other parameters in the given
        #       Species object match up to those in the Species to which
        #       they're being added

        # get individuals from the given Species object
        if individs is not None:
            individs = [ind for idx, ind in source_spp.items() if idx in individs]
        else:
            individs = [*source_spp.values()]

        # get their tskit.TableCollection
        source_tc = source_spp._tc

    # ... or simulate them w/ the given msprime args
    # and the recipient Species' genomic architecture
    else:
        assert individs is None, ("The 'individs' argument indicates the "
                                  "indices of Individuals to be added from "
                                  "a gnx Species object serving as a source "
                                  "population, so can only be used if "
                                  "'source_msprime_params' is None and "
                                  "'source_spp' is provided.")
        # TODO: CHECK MODEL INVOLVES NO SELECTION, RETURN INFORMATIVE WARNING!

        # TODO: check that all genomic and other parameters in the msprime
        #       match/jive with the parameters of the recipient Species

        # TODO: SET UP AND VALIDATE/2x-CHECK MSPRIME ANCESTRY SIMS!

        # TODO: GENERALIZE/RECYCLE CODE ALREADY WRITTEN TO SIMULATE MUTS?

        # TODO: GENERALIZE/RECYCLE CODE ALREADY WRITTEN TO CREATE gnx Individs
        #       FROM THOSE RESULTS (BUT IGNORE x,y, ETC, TO BE SET BELOW)

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

    # TODO: reset source Individuals' indices,
    #       and increment the recipient Species's max idx in lock step

    # TODO: reset/set source Individuals' environmental and fitness values

    # TODO: CHECK FOR ANY OTHER INDIVID- OR SPP-LEVEL ATTRS NEEDING RESET!

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





# OTHER TODO:

    # : DOUBLE CHEKC IM UNDERSTANDING AND USING POPULATIONS TABLE AND
    # nodes['population'] FIELD CORRECTLY!

    # : UPDATE GNX GENERALLY SO THAT IT IS PROPERLY FILLING IN THE
    # provenances INFO IN ANY NEW TABLECOLLECTION!


    # how does the 'location' column in individuals table work again?

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


def _cull_individuals(spp, ids):
    """
    remove the Individuals in the given list of ids from the given Species
    """


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



