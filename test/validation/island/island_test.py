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
    3. Mean movement and dispersal distances are set to 1/6th of island width
    and interisland distance (i.e. 1 cell). Movement is isotropic (i.e. no
    _MovementSurface is used).
    4. Carrying capacity is set to 10x the Layer 0 raster (i.e. the islands)

    This test uses the parameters file 'island_params.py', in this directory.
    """
    def make_island_landscape(self, n, w, d):
        """
        n = number of islands
        w = island width
        d = interisland diagonal distance (i.e. distance/sqrt(2))
        """
        #determine land dimension (i.e. nrows & ncols, hence n diagonal cells)
        dim = int(n * (w + d) - d)
        #create landscape
        land = np.zeros([dim]*2)
        island_ul_corners = [int(i) for i in np.linspace(0, dim-w, n)]
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
        #create a Model from the params file
        mod = gnx.make_model('island_params.py')
        #replace Layer 0's raster with an island raster (6 islands, each 6
        #cells wide and 6 cells diagonally spaced), and Layer 1's raster
        #with a raster containing island labels
        islands, island_labels = self.make_island_landscape(6,6,6)
        mod.land[0].rast = islands
        #run the model
        try:
            #burn the model in
            mod.walk(mode = 'burn', T = 1000000)
            #empty list to hold all numbers of migration events
            mig_ratios = []
            #walk it one timestep at a time
            for t in range(mod.params.main.T):
                #record individuals' current island locations
                #(i.e. their second habitat values)
                island_labels = [ind[1] for ind in mod.comm[0]._get_e()]
                island_labels = dict(zip([*mod.comm[0]],
                    [i[1] for i in mod.comm[0]._get_e()]))
                mod.walk(mode = 'main', T = 1)
                #record the number of individuals whose island numbers
                #changed (i.e. who migrated this timestep)
                new_island_labels = dict(zip([*mod.comm[0]],
                    [i[1] for i in mod.comm[0]._get_e()]))
                not_newborns = [ind for ind in [
                    *island_labels] if ind in [*new_island_labels]]
                migrated = [new_island_labels[ind][0] != island_labels[ind][
                    0] for ind in not_newborns] 
                mig_ratio = sum(migrated)/len(not_newborns)
                mig_ratios[t] = mig_ratio

            #use the number migrated each timestep to estimate the 
            #model's migration rate, then use that to compare the genetic
            #data to the result expected under the 1d island model
        except:
            print("Island model could not be run.")


if __name__=='__main__':
    unittest.main()
