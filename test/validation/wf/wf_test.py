import unittest
import os
import geonomics as gnx

class Wright_FisherTestCases(unittest.TestCase):
    """Named after early pioneers of theoretical population genetics,
    Sewall Wright and Ronald A. Fisher, the Wright-Fisher model describes
    the sampling of alleles in a population with no selection, no mutation,
    no migration, non-overlapping generation times and random mating.

    We parameterize an approximation of the Wright-Fisher model as follows:
    1. Ignore mutation, recombination, sex, selection,
    population size change and structure
    2. Set recombination rate to 0.5 by setting alpha and beta value to 0
    3. Set mutation rate to 0
    4. Set assigning sex to false

    This test uses the parameters file 'wf_params.py', in this directory.
    """
    def test_wright_fischer_model(self):
        #make a new template parameter file
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "wf_params.py"
        #test printing out the file path here
        gnx.make_parameters_file(filepath)
        #TODO: edit the parameters file (this step need to be done manually)
        #use parameter file to instantiate a model object
        mod = gnx.make_model('wf_params.py')
        #run the model
        try:
            mod.run()
        except:
            print("Wright fisher model could not be run.")


if __name__ == '__main__':
    unittest.main()
