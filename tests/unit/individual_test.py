import unittest
import os
import geonomics as gnx


class IndividualTestCases(unittest.TestCase):
    """
    Unit tests for Individual.py.
    """
    def test_make_individual_simple_case(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GNX_test_params.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = gnx.make_landscape(params)
        genome_arch = gnx.structs.genome._make_genomic_architecture(params,
                                                                    land)
        ind = gnx.structs.individual._make_individual(1, True, None,
                                                      genome_arch)
        self.assertIsInstance(ind, gnx.structs.individual.Individual)
        print("It's an individual!")

    def test_make_individual_largeNumber(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GNX_test_params.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = gnx.make_landscape(params)
        genome_arch = gnx.structs.genome._make_genomic_architecture(params,
                                                                    land)
        ind = gnx.structs.individual._make_individual(1000000, True, None,
                                                      genome_arch)
        self.assertIsInstance(ind, gnx.structs.individual.Individual)
        print("It's an individual with a large index number!")


if __name__ == '__main__':
    # from structs import genome
    unittest.main()
