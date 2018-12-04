import unittest
from sim import model
from structs import landscape
from structs import population
import geonomics as gnx
import os

class PopulationTestCases(unittest.TestCase):
    """
    Unit tests for population.py. 
    """
    def test_make_population(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = gnx.make_landscape(params)
        #TODO: Need to update here
        #in order to make population we need to have land, params, pop_params
        pop = population._make_population(land, "aribitrary name", params)
        self.assertIsInstance(pop, population.Population)


if __name__ == '__main__':
    unittest.main()
