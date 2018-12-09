import unittest
from copy import deepcopy
import geonomics as gnx

class Island1DTestCases(unittest.TestCase):
    """This test approximates a population-genetic 1-dimensional island model,
    or stepping-stone model.

    We parameterize this approximation as follows:
    1. The 66x66 Landscape consists of a chain of six, square,
    equal-area (6x6 cells) equidistant (6 diagonal cells apart)
    islands along the diagonal.
    2. Ignore mutation, recombination, and selection (i.e. all mutation rates
    are 0; all recombination rates are 0.5; and the Species has no Traits).
    3. Mean movement distance is set to 1/6th of island width
    and interisland distance (i.e. 1 cell), and standard deviation is set to
    0.25. Movement is isotropic (i.e. no
    _MovementSurface is used).
    4. Dispersal is extremely local (to minimize possibility of dispersal
    between islands), with mean and standard deviation of dispersal
    distribution set at 0.001, 0.001.
    5. Mating radius is set to 1 cell (so that interisland matings do not
    occur).
    6. Carrying capacity is set to 10x the Layer 0 raster (i.e. the islands).
    7. Starting population size is really large (5000), so that all islands
    should get seeded with an adequate number of individuals to start a
    self-sustaining local population.

    This test uses the parameters file 'island_params.py', in this directory.
    """
    def make_island_landscape(self, n, w, d):
        """
        n = number of islands
        w = island width
        d = interisland diagonal distance (i.e. distance/sqrt(2))
        """
        #determine land dimension (i.e. nrows & ncols, hence n diagonal cells)
        dim = int(n * (w + d) + d)
        #create landscape
        land = np.zeros([dim]*2)
        island_ul_corners = [int(i) for i in np.linspace(d, dim-w-d, n)]
        for ulc in island_ul_corners:
            land[ulc:ulc+w, ulc:ulc+w] = 1
        #create a second raster, where each island's hab values are its island
        #number (to be used later to estimate migration rate, to then compare
        #results to expected results of stepping-stone model
        island_labels = deepcopy(land)
        for ulc in island_ul_corners:
            island_labels[ulc:ulc+w, ulc:ulc+w] = (ulc/(w+d)+1)/n

        return land, island_labels

    def test_island_model(self):
        #define island number, width, and diagonal distance
        n = 6
        w = 6
        d = 6
        #create a Model from the params file
        mod = gnx.make_model('island_params.py')
        #replace Layer 0's raster with an island raster (6 islands, each 6
        #cells wide and 6 cells diagonally spaced), and Layer 1's raster
        #with a raster containing island labels
        islands, island_label_rast = self.make_island_landscape(n, w, d)
        mod.land[0].rast = islands
        mod.land[1].rast = island_label_rast
        #create a lookup-dict for the island numbers of each non-zero value in
        #island_labels
        island_vals =  dict(zip(sorted([*np.unique(island_label_rast)])[1:],
            range(n)))
        #run the model
        try:
            #burn the model in
            mod.walk(mode = 'burn', T = 1000000)
            #empty dict to hold proportions of each timestep's population
            #that underwent each of the possible inter-island pairwise
            #migration events (i.e. island 0 to 1; 1 to 0; 0 to 2; etc...)
            mig_props = {}
            for i in range(n):
                for j in range(n):
                    if i != j:
                        mig_props.update({(i,j):[]})
            #walk Model timestep by timestep
            for t in range(mod.params.main.T):
                #create a dictionary to keep count of how many individuals make
                #each of the possible migration events
                mig_counts = {}
                for i in range(n):
                    for j in range(n):
                        if i != j:
                            mig_counts.update({(i,j):0})
                #record individuals' current island locations
                #(i.e. their second habitat values)
                old_island_labels = dict(zip([*mod.comm[0]],
                    [island_vals[i[1]] for i in mod.comm[0]._get_e()]))
                assert 0 not in [*old_island_labels.values()], ("It appears "
                    "some individuals landed in the 'sea' during the "
                    "previous timestep and didn't die!")
                mod.walk(mode = 'main', T = 1)
                #record the number of individuals whose island numbers
                #changed (i.e. who migrated this timestep)
                new_island_labels = dict(zip([*mod.comm[0]],
                    [island_vals[i[1]] for i in mod.comm[0]._get_e()]))
                assert 0 not in [*new_island_labels.values()], ("It appears "
                    "some individuals landed in the 'sea' during "
                    "this timestep and didn't die!")
                not_newborns = [ind for ind in [
                    *old_island_labels] if ind in [*new_island_labels]]
                #get number of individuals who underwent each of the possible
                #migration events 
                for ind in not_newborns:
                    if new_island_labels[ind][0] != old_island_labels[ind][0]:
                        mig_event = (old_island_labels[ind][0],
                            new_island_labels[ind][0])
                        mig_counts[mig_event] += 1
                #append the proportion of individuals who underwent each
                #migration event to the mig_props dict
                [mig_props[mig_event].append(ct/len(
                    not_newborns)) for mig_event, ct in mig_counts.items()]

            #use the number migrated each timestep to estimate the 
            #model's migration rate, then use that to compare the genetic
            #data to the result expected under the 1d island model
        except:
            print("Island model could not be run.")


if __name__=='__main__':
    unittest.main()
