import numpy as np
import geonomics as gnx
import msprime
import tskit
import warnings
from collections import Counter as C



# TODO PICKUP:

    # debug functionality to add individs from another gnx model (getting
    # duplications of the same individual idx, which means I just to rset source
    # individuals' idxs on them and in their tc

    # why check_shared_equality failing?

    # move code to proper methods

    # make sure model plots correctly after adding individs,
    # once the fns are proper methods
    # (because the spp was becoming disassociated from the model when I was
    # playing around)

    # keep debugging with all various combos of args
    # TODO: TEST WITH BOTH DEBUG ARGS==TRUE IN SORT_AND_SIMPLIFY FN!

    # 3. work out function and params-handling for instantiating population
    #    using one or more msprime pops, then debug and test that

    # 4. integrate everything into pkg

    # 5. produce a clean, simple example nb or two (AND TRY TO PUT ON SITE!)

    # 6. clean up repo and publish new version

    # 7. circle back to group with news and script



def _format_msprime_rate(rate, L):
    """
    Either return a simple float, if the given rate
    is homogeneous, or else return an msprime.RateMap object.
    Used for msprime recombination and mutation rates.
    """
    assert (type(rate) in [float, list, np.ndarray] or rate == 0 or rate == 1)
    if (isinstance(rate, float) or rate == 0 or rate == 1):
        return rate
    elif len(np.unique(rate)) == 1:
        return rate[0]
    else:
        assert len(rate) == L
        # NOTE: position array has to be one longer than the sequence length,
        #       so that rightmost position brackets a stretch of genomic
        #       intervals equal to the sequence_length
        positions = np.linspace(0, L, L+1)
        ratemap = msprime.RateMap(position=positions, rate=rate)
        return ratemap


def sim_msprime_individs(n,
                         L,
                         recomb_rate,
                         mut_rate,
                         demography=None,
                         population_size=None,
                         ancestry_model=None,
                         random_seed=None,
                         gnx_spp_use_tskit=True,
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

    ancestry_model: {`msprime.AncestryModel`, list}, optional, default: None
        The `msprime.AncestryModel instance (or list of them) to be used by
        `msprime.sim_ancestry`. See `msprime.sim_ancestry` documentation and
        references therein for more details.

    random_seed: int, optional, default: None
        The random seed value to be provided to `msprime.sim_ancestry` and
        `msprime.sim_mutations`. Use if reproducibility is desirable.
        See `msprime.sim_ancestry` and `msprime.sim_mutations` documentation
        for detail.

    gnx_spp_use_tskit: bool, optional, default: True
        Whether or not the Geonomics `Species` for which the `Individual`s are
        intended uses tskit (i.e., has Species.gen_arch.use_tskit as True).
        If not, `Individual`s will each store a copy of their Lx2 genotype
        array in the 'g' attribute; otherwise, this attribute will be None.

    Returns
    -------
    :class:`list`
        A list of `geonomics.structs.individual.Individual` objects.
        Individuals are assigned sequential index integers starting from 0,
        placed at x,y coordinates (0, 0), and otherwise given default attribute
        values. Attributes should be updated as necessary in downstream code.

    '''
    # handle the recombination and mutation rates
    if (not isinstance(recomb_rate, float) and
        not (recomb_rate == 0)):
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
    recomb_rate = _format_msprime_rate(recomb_rate, L)
    mut_rate = _format_msprime_rate(mut_rate, L)

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
                              discrete_genome=True,
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
    individs = [gnx.structs.individual.Individual(idx=0,
                                                  x=0,
                                                  y=0) for i in range(n)]

    # get the information linking Individuals to rows in TableCollection tables
    nodes_tab_individs_col = mtc.nodes.asdict()['individual']
    for i, ind in enumerate(individs):
        # get the nodes table rows corresponding to this Individual
        nodes_tab_ids = np.where(nodes_tab_individs_col==i)[0]
        assert len(nodes_tab_ids) == 2
        # set hidden attributes tracking the individs' rows in the 'nodes' table
        # and id in the 'individuals' table
        ind._set_nodes_tab_ids(nodes_tab_ids)
        ind._individuals_tab_id = i
        # if the gnx Model these Individuals are intended for is not using tskit
        # then walk the variants generator and compile the individ's
        # genotype array, then store it in the attribute 'g'
        if not gnx_spp_use_tskit:
            var_gen = mts.variants(samples=nodes_tab_ids)
            gt = np.stack([next(var_gen).genotypes for _ in range(L)])
            ind.g = gt

    return individs, mtc



############
# Genome fn?

def _prep_tablecollection_for_gnx_spp_union(tc,
                                            recip_spp,
                                            coords,
                                           ):
    '''
    Reformat a tskit.TableCollection produced by an msprime simulation
    so that it conforms to Geonomics' expected format, and thus can be fed
    into the .union() method of the given Geonomics Species' current TableCollection.

    The code and comments walk through each table, column by column,
    and specify if and how the contents need to be reformatted or updated
    in order to match Geonomics conventions.
    '''

    # reformat nodes table content:
    #-----------------------------

    # TODO: DOUBLE CHECK (REALLY 0 = STARTING GNX TIME STEP?)
    # gnx uses the time column as follows:
        # - smallest values are the most recent individuals
        # - nodes with val 0 originated at very beginning of sim's main
        # - after that, new nodes are assigned a time of origination = -1*spp.t
    # msprime times are expressed as time before present (i.e. samples are 0s).
    # thus, convert msprime to gnx convention by subtracting from the whole
    # column the value of spp.t at the time of their introduction
    # (i.e., the time of their simulation and sampling)
    time = tc.nodes.time - recip_spp.t
    tc.nodes.time = time

    # msprime uses the flags column to indicate extant nodes
    # (i.e., docs say "NODE_IS_SAMPLE = 1"); gnx uses the flags to distinguish
    # between nodes in 'real' gnx Individuals and msprime nodes;
    # thus, we actually just want to keep this column as is, because the
    # samples are the only nodes that will be assigned to new gnx Individuals,
    # and thus are the only ones we would otherwise coerce to 1s anyhow!

    # nothing to do the population column; TableCollection.union()
    # will be asked to use the next available int in the recipient Species'
    # populations table to automatically record the population these nodes came
    # from

    # the individual column needn't be touched; TableCollection.union()
    # will automatically increment this column by using the next individual id
    # available in the recipient Species' individuals table;
    # we will just need to later track what the next id was and how many
    # individuals were added, so that TableCollection.individuals ids can
    # be matched to gnx Individuals via the Individuals._individuals_tab_id attr

    # metadata, metadata_offset, and metadata_schema columns are currently
    # unused by both msprime and gnx


    # reformat individuals table content:
    #-----------------------------------

    # all individuals in the individuals table are extant and will be turned
    # into gnx Individual objects, and gnx uses a 1 in the 'flags' column to
    # indicated 'real' Individuals (i.e., individuals that were actual
    # gnx.Individuals at some point), so set the entire 'flags' column to 1s
    flags = tc.individuals.flags*0+1

    # gnx stores individs' starting x,y coordinates in the 'location' column
    # (NOTE: also the phenotype and fitness values, but those are not relevant
    # in the neutral-loci-only msprime simulations permitted by gnx);
    # TODO: think about if/how sweeps would be allowed and how represented/matched up to genarch
    # location_offset tells the TableCollection to jump 2 values (x, y) per row
    # NOTE: 'ragged columns' actually consist of the column's contents in a
    #       flattened array format and a column of n+1 index values that can be
    #       offset and used to index the start and end positions of each row's
    #       contents within that array
    location_offset = [int(n) for n in np.linspace(0,
                                                   2*tc.individuals.num_rows,
                                                   tc.individuals.num_rows+1)]

    # parents column is unused by geonomics
    # (wasn't implemented until 21/01/21, after gnx was developed;
    # tskit-dev commit 8ebfd8b2900217dab2b82d6f530c830238037fc8).
    # useful but not necessary to implement;
    # this column is also empty in msprime, so leaving empty for now

    # gnx uses the metadata column to store each individual's
    # gnx Individual idx, for later matching to rows in the individuals table
    # when the table is updated after tskit's simplify algo filters individuals;
    # determine the next n indx values using the Species.max_ind_idx attribute
    # and the length of the individuals table
    gnx_idxs = [*range(recip_spp.max_ind_idx+1,
                       recip_spp.max_ind_idx + tc.individuals.num_rows + 1)]
    metadata=[idx.to_bytes(length=4, byteorder='little') for idx in gnx_idxs]
    metadata_offset = [int(n) for n in np.linspace(0,
                                                   4*tc.individuals.num_rows,
                                                   tc.individuals.num_rows+1)]

    # NOTE the set_columns() method fails because the metadata column, although
    #      supposed to be a binary data type according to the docs, is actually
    #      coded as np.int8 and so fails to convert the 4-byte binary data we
    #      are trying to provide; never ran into this problem in tha past
    #      because the .add_row() method allows the binary data for some reason
    #      and that's the only method I've used before (just adding 1
    #      individual at a time as they're born);
    #      however, somewhat conveniently, we were already planning to change
    #      every column already containing data in this table anyhow, so for
    #      now we are just:
    #           1.) double-checking the length of the original table
    #           2.) clearing the table
    #           3.) manually looping and adding new rows
    #           4.) confirming identical resulting length
    assert (len(flags) ==
            coords.shape[0] ==
            (len(location_offset)-1) ==
            len(metadata) ==
            (len(metadata_offset)-1)), ("Tne 'flags', 'location', and "
                                        "'metadata' columns for the "
                                        "individuals table must all "
                                        "be the same length as the coords "
                                        "array's 0th dimension and the "
                                        "'location_offset' and "
                                        "'metadata_offset' columns must be "
                                        "one longer.")
    original_inds_tab_len = tc.individuals.num_rows
    assert original_inds_tab_len == coords.shape[0], ("Length "
                                             "of individuals table "
                                             "generated by msprime "
                                             "differs from length of coords "
                                             "attempting to be added to it.")
    tc.individuals.clear()
    for i, flag in enumerate(flags):
        tc.individuals.add_row(flags=flag,
                               location=coords[i, :],
                               parents=None,
                               metadata=metadata[i],
                              )
    assert tc.individuals.num_rows == original_inds_tab_len, ("Reconstructed "
                                            "individuals table is not the same "
                                            "length as the original table "
                                            "generated by msprime.")


    # reformat edges table content:
    #-----------------------------

    # geonomics models all edges as being delineated by positions belonging
    # to the set [0, 0.5, 1.5, ..., L-2.5, L-1.5, L];
    # this is more complicated than the saner msprime
    # convention of using the set [0, 1, 2, ..., L-2, L-1, L],
    # but this was a decision I made based on modeling recombination
    # breakpoints halfway between loci and it is equally functional, so...
    # c'est la vie.
    L = recip_spp.gen_arch.L
    left = tc.edges.left
    right = tc.edges.right
    left = np.clip(left - 0.5, a_min=0, a_max=L)
    right = np.clip([r-0.5 if r != L else r for r in right], a_min=0, a_max=L)
    tc.edges.left = left
    tc.edges.right = right

    # parent and child columns are both fine as is
    # (have double-checked that node ids are all automatically and correctly
    # updated by the tskit.TableCollection.union() method)

    # metadata column is currently unused by both msprime and gnx


    # reformat sites table content:
    #-----------------------------

    # double-check that position and ancestral state columns are identical in
    # the recipient Species' TableCollection;
    # if so, then just duplicate that table here
    # (because the only difference between the two in the metadata column,
    # which in gnx stores whether a site started as a neutral ('n')
    # or non-neutral ('t', for trait-influencing) site
    # TODO: ACTUALLY DON'T EVEN RUN THESE CHECKS BECAUSE SITES TABLE
    #       CAN LACK SOME SITES DEPENDING ON RECOMB RATE, N_E, ETC
    #print(tc.sites)
    #print(recip_spp._tc.sites)
    #for col in ['position', 'ancestral_state']:
    #    assert np.all(getattr(tc.sites, col) ==
    #                  getattr(recip_spp._tc.sites, col)), (
    #                        f"The {col} column in the sites table "
    #                        f"generated by msprime disagrees with the {col} "
    #                         "column in the recipient Species' sites table.")
    #tc.sites.replace_with(recip_spp._tc.sites)


    # reformat mutations table content:
    #---------------------------------

    # site, node, and parent columns are the core pieces needed to recover the
    # genotypes produced by the mutations, and we want whatever genotypes
    # msprime has simulated, so these columns should all remain untouched
    # site vals are all in interval [0, L-1]
    assert (np.all(0 <= tc.mutations.site) and
            np.all(tc.mutations.site <= (recip_spp.gen_arch.L-1))), ("The "
                        "site column in the mutations table contains values "
                        "not in the interval [0, recip_spp.gen_arch.L-1].")

    # Mutation times in msprime are expressed in generations before present
    # (i.e., larger values are further toward the root of the coalescent tree).
    # gnx also expresses time this way, but in gnx 0 demarcates the start of the
    # forward-time simulation, and thus negative values represent -1 times the
    # gnx time step. Thus, to convert to gnx time we just need to subtract the
    # current gnx time step from msprime's time column (same as how we handled
    # time for the nodes table, above; and in fact, if we didn't do this we
    # would set ourselves up to fail to satisfy the tskit requirements
    # for the comparison between the time of a mutation and the times of the
    # nodes surrounding it in the tree
    # Times for mutations created at the start of a gnx simulation, to
    # match the starting allele frequency array ('p') provided to a gnx
    # model, are artificial as they are all forcibly mutated at at once.
    # Thus, those times are NaN. Unfortunately, tskit does not all
    # mutation times to be mixed known and unknwon, so that means we must
    # set all mutation times to NaN
    time = [tskit.UNKNOWN_TIME] * len(tc.mutations.time)
    tc.mutations.time = time

    # double-check that the derived_state_offset column is just serial
    # integers from 0 to the table length, inclusive
    assert np.all(tc.mutations.derived_state_offset ==
                  np.linspace(0,
                              tc.mutations.num_rows,
                              tc.mutations.num_rows + 1)), ("The "
                "derived_state_offset column of the mutations table should "
                "contain serial integers from 0 to the table length, "
                "inclusive, but does not.")

    # metadata, metadata_offset, and metadata_schema are all currently unused
    # in both msprime and gnx


    # reformat migrations table content:
    #-----------------------------------

    # this table is not used in gnx, but if it was used in a complex msprime
    # Demography that simulated individuals for introduction into gnx then it
    # will not be empty in this TableCollection; in that case, its contents
    # will just be unioned into the otherwise empty migrations table of the 
    # recipient Species and then left as is


    # reformat populations table content:
    #-----------------------------------

    # This table is unused by gnx, but the msprime table has
    # basic population names (e.g., 'pop_0') and descriptions ('' by default,
    # but a user could specify).
    # We won't touch or do anything with this data, and the
    # TableCollection.union() method will pull it along without editing,
    # but the rows in this table will be appended to existing rows in the gnx
    # model and the nodes table will be constructed so as to match those indices

    # reformat provenances table content:
    #-----------------------------------

    # no need to do anything here; all provenances info will be retained and
    # combined when the union() method is called on the recipient Species'
    # TableCollection

    # return the gnx idxs that were set in the individuals table, for
    # downstream validation
    return gnx_idxs





##############
# Model method
def add_individs(n,
                 coords,
                 land,
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
    if isinstance(recip_spp, int) or isinstance(recip_spp, str):
        recip_spp = self.comm[self._get_spp_num(recip_spp)]

    # and make sure it has already been burned in
    assert recip_spp.burned, ("Individuals cannot be manually added to a "
                              "Species that hasn't been burned in yet.")

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

    # get the next n gnx Individuals' ids
    next_n_idxs = [*range(recip_spp.max_ind_idx + 1,
                        recip_spp.max_ind_idx + n + 1)]

    # check and handle the coords param
    assert np.atleast_2d(coords).shape in [(1,2), (n,2)], (
        "The 'coords' arg must be provided either "
        "a single x,y coordinate pair (at which all Individuals "
        "will be introduced) or an nx2 array of x,y coordinate pairs "
        "(indicating the introduction locations for "
        "each of the n Individuals being introduced.")
    if np.atleast_2d(coords).shape == (1, 2):
        coords = np.stack([coords for _ in range(n)])
        assert coords.shape == (n, 2)
    # make sure coords are all within the landscape
    for coord in coords:
        assert (0 <= coord[0] <= recip_spp._land_dim[0] and
                0 <= coord[1] <= recip_spp._land_dim[1]), ("Coords must be "
                                                     "valid for the recipient "
                                                     "Species' Landscape "
                                                     "(i.e., x and y must be "
                                                     "in the inclusive "
                                                     " intervals "
                                                     "[0, land_dim_x] and"
                                                     "[0, land_dim_y], "
                                                     "respectively.")

    # either get Individuals and their tskit.TableCollection
    # from the given gnx Species object...
    if source_spp is not None:
        # raise warning about potential unforeseen effects of introducing
        # Individuals from an incorrectly parameterized source Species
        warnings.warn(("Introducing Individuals from one Species object "
                       "into another could cause unnoticed issues. "
                       "Geonomics will check and prevent conflict between strictly "
                       "incompatible attributes of the Species, "
                       "GenomicArchitectures, and Traits. However, attributes "
                       "whose compatibility depends on the scientific "
                       "reasoning behind the scenario being simulated "
                       "(e.g., mutational parameters, "
                       "carrying capacity parameters, etc.) "
                       "will be allowed to differ between the source and "
                       "recipient Species objects and, "
                       "if these differences were unintended but go "
                       "undetected, could lead to incorrect or misleading data "
                       "and results. Please use care when planning, "
                       "implementing, and checking your code."))

        # get the source Species from within this model, if necessary
        if isinstance(source_spp, int) or isinstance(source_spp, str):
            source_spp = self.comm[self._get_spp_num(source_spp)]

        # assert that exactly one of the 'n' and 'individs' params is provided,
        # and that the size of 'n' or size and contents of 'individs' are valid
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
                    np.all([ind in source_spp for ind in individs])), (
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
            individs = [*source_spp.values()][:n]

        # TODO UPDATE INDIVIDUALS TABLE TO REFLECT THESE NEW IDX VALS, EITHER
        # EHERE OR BELOW
        # get the list of indexes they should wind up having
        next_gnx_idxs = [*range(recip_spp.max_ind_idx + 1,
                                recip_spp.max_ind_idx + 1 + len(individs))]

        # get their tskit.TableCollection
        # TODO SUBSET TO ONLY THE INDIVIDUALS/NODES NEEDED!
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

        # run all TableCollection edits that are necessary to align msprime
        # output with the TableCollection conventions used in gnx
        # NOTE: edits the TableCollection in place;
        #       returns the next n gnx Individual idxs assigned to those
        #       Individuals for incorporation into the given Species
        next_gnx_idxs = _prep_tablecollection_for_gnx_spp_union(tc=source_tc,
                                                        recip_spp=recip_spp,
                                                        coords=coords,
                                                               )

    # set the source Individuals' coordinates
    [ind._set_pos(*coords[i,:]) for i, ind in enumerate(individs)]

    # update the index values identifying Individuals in gnx and tskit data
    max_idx_b4 = recip_spp.max_ind_idx
    new_idxs = []
    for i, ind in enumerate(individs):
        # increment the recipient Species' max_ind_idx
        recip_spp.max_ind_idx += 1
        # set the Individual's gnx index
        ind.idx = recip_spp.max_ind_idx
        new_idxs.append(recip_spp.max_ind_idx)
        # set the Individuals' ids for the tskit nodes table
        ind._nodes_tab_ids = {k: v + recip_spp._tc.nodes.num_rows for
                                            k, v in ind._nodes_tab_ids.items()}
        # set the Individuals' ids for the tskit individuals table
        # NOTE: the inverse (the Individuals' gnx idx values, as recored in the
        #       tskit individuals table's metadata) should already be correct
        #       because _prep_tablecollection_for_gnx_spp_union
        #       generated the index values and set them in the TableCollection
        ind._individuals_tab_id += recip_spp._tc.individuals.num_rows
    max_idx_af = recip_spp.max_ind_idx
    assert max_idx_af - max_idx_b4 == len(individs)
    # compare the idx values that were set on the Individuals to those that
    # were set in the TableCollection.individuals table
    assert np.all(np.array(next_gnx_idxs) == np.array(new_idxs)), ("Individaul "
                            "idx values set in the TableCollection.individuals "
                            "table do not agree with those set on the "
                            "Individuals themselves.")

    # sort and simplify the recipient Species' tskit.TableCollection
    recip_spp._sort_and_simplify_table_collection()

    # union source individuals' TableCollection into that of recipient Species
    nodes_tab_len_b4 = recip_spp._tc.nodes.num_rows
    individs_tab_len_b4 = recip_spp._tc.individuals.num_rows
    recip_spp._tc.union(source_tc,
                        # NOTE: (for now at least) assuming that the introduced
                        #       individuals come from totally independent pop,
                        #       do not share a MRCA within our data, and thus
                        #       have no common nodes that need to be mapped,
                        #       hence we provide a list of tskit.NULLs
                        #       to union()'s node_mapping argument
                        node_mapping=[tskit.NULL] * source_tc.nodes.num_rows,
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
                        # TODO: FIGURE OUT WHY FAILING WHEN MSPRIME SIM AND GNX
                        #       SIM CLEARLY HAVE NO SHARED PORTIONS!
                        check_shared_equality=False,
                        # NOTE: should copy provenance info over from the source
                        #       Individuals' TableCollection to that of the
                        #       receipient Species
                        record_provenance=True,
                             )
    nodes_tab_len_af = recip_spp._tc.nodes.num_rows
    individs_tab_len_af = recip_spp._tc.individuals.num_rows

    # check unioned tables have correct lengths
    assert (nodes_tab_len_af - nodes_tab_len_b4) == source_tc.nodes.num_rows
    assert ((individs_tab_len_af - individs_tab_len_b4) ==
                                                source_tc.individuals.num_rows)

    # check that all gnx idxs stored in the individuals table's metadata column
    # occur only once (aside from 0, which is a null value for individuals that
    # never existed as 'real' gnx Individual objects because they were in the
    # past of a coalescent simulation)
    recip_tc_idxs = [int.from_bytes(recip_spp._tc.individuals[i].metadata,
                  'little') for i in range(recip_spp._tc.individuals.num_rows)]
    assert len(recip_tc_idxs) == recip_spp._tc.individuals.num_rows
    counts = C(recip_tc_idxs)
    for k, v in counts.items():
        if k != 0:
            assert v == 1, (f"gnx Individual idx {i} occurs more than once "
                             "in the metadata column of the union of "
                             "the recipient and source Species' "
                             "TableCollections individuals tables.")

    # add the Individuals to the recipient Species's dict, keyed by idx
    pop_size_b4 = len(recip_spp)
    for ind in individs:
        recip_spp[ind.idx] = ind
    pop_size_af = len(recip_spp)
    assert (pop_size_af - pop_size_b4) == len(individs), ("Incorrect number "
                            "of Individuals added to recipient Species.")

    # update each Individual's _individuals_tab_id and _nodes_tab_ids
    # to match the new, unioned TableCollection
    nodes_tab_individuals = recip_spp._tc.nodes.individual
    for idx, ind in recip_spp.items():
        print(np.argwhere(np.array(recip_tc_idxs) == idx))
        individuals_tab_id = np.argwhere(np.array(recip_tc_idxs) == idx).ravel()
        assert len(individuals_tab_id) == 1
        individuals_tab_id = individuals_tab_id[0]
        ind._individuals_tab_id = individuals_tab_id
        nodes_tab_ids = np.argwhere(nodes_tab_individuals ==
                                    individuals_tab_id).ravel()
        assert len(nodes_tab_ids) == 2
        ind._nodes_tab_ids = {i: id for i, id in enumerate(nodes_tab_ids)}

    # final check that indices check out in the other direction
    # (i.e., from nodes tab to individuals tab to gnx Individuals
    # and that nodes ids check out too
    curr_node_ct = 0
    for node_id, node in enumerate(recip_spp._tc.nodes):
        # NOTE: ignore -1s; they are msprime flags for non-sample nodes but
        #       beahve as a Python indexer of the last Individual in the
        #       Species and thus throw an error
        if node.individual != -1:
            gnx_individual_id = int.from_bytes(recip_spp._tc.individuals[
                                            node.individual].metadata, 'little')
            if gnx_individual_id in [*recip_spp]:
                curr_node_ct += 1
                assert (recip_spp[gnx_individual_id]._individuals_tab_id ==
                        node.individual and node_id in
                        recip_spp[gnx_individual_id]._nodes_tab_ids.values())
    assert curr_node_ct == 2*len(recip_spp), ("Number of current Individuals "
                                           "in the recipient Species does not "
                                           "match the number expected based on "
                                           "TableCollection contents.")

######################
# back to model method
    # sort and simplify the recipient Species' tskit.TableCollection again
    recip_spp._sort_and_simplify_table_collection()

    # reset the Species' coords and cells attrs
    recip_spp._set_coords_and_cells()

    # reset the Species' KDTree after adding the new individuals
    recip_spp._set_kd_tree()

    # reset the density grids
    recip_spp._set_dens_grids(self.land)

    # reset/set source Individuals' environmental and fitness values
    recip_spp._set_e(self.land)

    #reset phenotypes and fitnesses
    if recip_spp.gen_arch is not None and recip_spp.gen_arch.traits is not None:
        [ind._set_z(recip_spp.gen_arch) for ind in recip_spp.values()]
        recip_spp._calc_fitness(set_fit=True)

    # print informative output
    print(f"\n{len(individs)} Individuals successfully "
          f"added to Species {recip_spp.idx} ('{recip_spp.name}').\n")




#
##
##source_pops = {'A': {'initial_size':10000,
##                     'growth_rate':0,
##                     'description':'first'},
##               'B': {'initial_size':1000,
##                     'growth_rate':0.1,
##                     'description':'second'},
##               'C': {'initial_size':100,
##                     'growth_rate': -0.1,
##                     'description':'last',
##                    },
##              }
##
##source_pop_intros = {'A': [25, (4, 6)],
##                     'B': [50, (1, 3)],
##                     'C': [100, (5,6)],
##                    }
##
##source_pop_allele_freqs = {'A': [0.1]*100,
##                           'B': [0.9]*100,
##                           'C': [0]*100,
##                          }
##
##
##def _make_msprime_demographies(source_pops):
##    """
##    use all of the distinct source populations identified in the provided
##    parameters dict to create msprime.Demography objects that
##    can be used to draw new gnx individuals
##    """
##    demogs = []
##    for name, params in source_pops.items():
##        disallowed_params = ['default_sampling_time',
##                             'initially_active',
##                            ]
##        # NOTE: need to disallow certain parameters that might break gnx models
##        #       render output inaccurate
##        for dp in disallowed_params:
##            assert dp not in params, ((f"parameter '{dp}' is not allowed in "
##                                        "msprime simulations used to generate "
##                                        "geonomics models; please reformulate "
##                                        "population parameters, then rerun."))
##        demog = msprime.Demography()
##        demog.add_population(name=name, **params)
##        demogs.append(demog)
##    return demogs
##
##
##def _get_msprime_demographies(source_pops):
##    """
##    either grab the msprime.Demography object
##    (or list of msprime.Demography objects) provided in the gnx parameters file
##    (for more complicated models involving e.g., migrations, splits, etc.),
##    or grab the demographic parameters dict provided in the gnx parameters file
##    and use it to create a list of msprime.Demography objects
##    """
##    if (isinstance(source_pops, msprime.Demography) and
##        source_pops.num_populations > 0):
##        return [source_pops]
##    elif (isinstance(source_pops, list) and
##          set([type(sp) for sp in source_pops]) == set([msprime.Demography])):
##        return source_pops
##    elif isinstance(source_pops, dict):
##        return _make_msprime_demographies(source_pops)
##    else:
##        raise ValueError((f"'source_pops' component of species 'init' "
##                           "parameters must be either an msprime.Demography "
##                           "with at least one population or a dictionary of "
##                           "name-keyed demographic parameters that can be used "
##                           "to instantiate msprime.Population objects for a "
##                           "Demography object."))
##
##
##def _add_individs_from_msprime_source_pops(n_individs,
##                                           demogs,
##                                           spp,
##                                           mut_rates=None
##                                          ):
##    """
##    simulate the numbers of individuals indicated by n_individs
##    (either an int, if demogs has only 1 population,
##    or else a population-name-keyed dict of ints)
##    using the given msprime.Demography and gnx.GenomicArchitecture objects,
##    and adding then to the given tskit.TableCollection
##    """
##    assert (len(spp.gen_arch.nonneut_loci) == 0 and
##            spp.gen_arch.traits is None), ("msprime can only currently "
##                                           "be used to seed a model with "
##                                           "individuals from source "
##                                           "populations if the model "
##                                           "has no selection. "
##                                           "This constraint may be "
##                                           "relaxed in the future.")
##    assert spp.gen_arch.use_tskit, ("to use msprime to seed a population with "
##                                    "individuals from source populations, "
##                                    "please set the 'use_tskit' parameter to "
##                                    "True.")
##    assert ((isinstance(n_individs, int) and len(demogs) == 1) or
##            (isinstance(n_individs, dict) and
##                (len(demogs) == 1 and
##                 isinstance(demogs[0], msprime.Demography) and
##                 demogs[0].num_populations == len(n_individs)) or
##                (len(n_individs) == len(demogs)))):
##    if mut_rates is not None:
##        assert len(mut_rates) == len(demogs)
##    # get the tree-seq object for the simulated ancestry
##    if len(np.unique(spp.gen_arch.recombinations._rates[1:])) == 1:
##        recomb_rate = spp.gen_arch.recombinations._rates[1]
##    else:
##        # NOTE: position array has to be one longer than the sequence length,
##        #       so that rightmost position brackets a stretch of genomic
##        #       intervals equal to the sequence_length
##        recomb_rate= msprime.RateMap(position=np.linspace(0,
##                                                          spp.gen_arch.L,
##                                                          spp.gen_arch.L+1,
##                                                         ),
##                              rate=[*spp.gen_arch.recombinations._rates[1:]]+[0.5])
##    if len(demogs) == 1 and isinstance(demogs[0], msprime.Demography):
##        ts = msprime.sim_ancestry(samples=n_individs,
##                                  demography=demogs[0],
##                                  sequence_length=spp.gen_arch.L,
##                                  recombination_rate=recomb_rate,
##                                  ploidy=2,
##                                  # NOTE: DiscreteTimeWrightFisher ('dtwf') may be
##                                  #       closest match to the discrete operations
##                                  #       of gnx's forward-time simulations, but
##                                  #       probably not necessary to match that and
##                                  #       probably simpler, faster, and most
##                                  #       reasonably expected to just use the
##                                  #       standard coalescent
##                                  model='hudson',
##                                 )
##        tcs = [ts.dump_tables()]
##    else:
##        tcs = []
##        i = 0
##        for n, demog in zip(n_individs, demogs):
##            ts = msprime.sim_ancestry(samples=n,
##                                      demography=demog,
##                                      sequence_length=spp.gen_arch.L,
##                                      recombination_rate=recomb_rate,
##                                      ploidy=2,
##                                      model='hudson',
##                                     )
##            # add mutations according to the gen_arch
##            # (using starting allele frequencies as the mutation rate,
##            # so that mutations will wind up approximating them)
##           if mut_rates is None and len(np.unique(spp.gen_arch.p)) == 1:
##                mut_rate = spp.gen_arch.p[0]
##                # don't bother with mutations if mut_rate == 0
##                if mut_rate == 0:
##                    break
##           else:
##               if mut_rates is not None:
##                   mut_rate_vals = mut_rates[i]
##                else:
##                    mut_rate_vals = spp.gen_arch.p
##                # NOTE: position array has to be one longer than the sequence length,
##                #       so that rightmost position brackets a stretch of genomic
##                #       intervals equal to the sequence_length
##                mut_rate= msprime.RateMap(position=np.linspace(0,
##                                                          spp.gen_arch.L,
##                                                          spp.gen_arch.L+1,
##                                                         ),
##                                          rate=mut_rate_vals,
##                                         )
##            mts = msprime.sim_mutations(tree_sequence=ts,
##                                        rate=mut_rate,
##                                        # state-independent binary mutation
##                                        # means that all initial mutations
##                                        # will create 1s, exactly as we want
##                                        # when we want to create
##                                        # specific initial allele frequencies
##                                        model=msprime.BinaryMutationModel(state_independent=True),
##                                       )
##            tcs.append(mts.dump_tables())
##            i += 1
##
##    # TODO: need some assert statements that can be run based on comparing
##    #       calculations using spp._tc before and after the additions!
##
##    # for each TableCollection
##    # (i.e., each introduction event from a source population)
##    for tc in tcs:
##        # increment all node and edge IDs to be > max node and edge IDs in the
##        # species' current TableCollection
##
##        # add new individuals to the TableCollection, matching parent nodes, etc
##
##        # add new individuals to the species, using correct TableCollection
##        # ids, correction geolocation, etc
##
##        # sort and simplify spp._tc
##
##    return
##
##
##
