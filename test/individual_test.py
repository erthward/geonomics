import unittest
import os
import geonomics as gnx
import structs
from structs import genome
from sim import model
from structs import individual



class IndividualTestCases(unittest.TestCase):
    """
    Unit tests for Individual.py. 
    """
    def test_make_individual_simple_case(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = gnx.make_landscape(params)
        genome_arch = genome._make_genomic_architecture(params, land)
        ind = individual._make_individual(1, True, None, genome_arch)
        self.assertIsInstance(ind, individual.Individual)

    def test_make_individual_largeNumber(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = gnx.make_landscape(params)
        genome_arch = genome._make_genomic_architecture(params, land)
        ind = individual._make_individual(1000000, True, None, genome_arch)
        self.assertIsInstance(ind, individual.Individual)

if __name__ == '__main__':
    from structs import genome
    unittest.main()
