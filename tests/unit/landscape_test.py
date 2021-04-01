import unittest
from sim import burnin
from sim import model
import os
import geonomics as gnx
from structs import landscape
import copy



class LandscapeTestCases(unittest.TestCase):
    """
    Unit tests for Landscape.py. 
    """
    def test_make_random_scape(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        dim = params.land.main.dim
        try:
            random_scape = landscape._make_random_lyr(dim, 10)
        except:
            print("Can not make the random scape")

    def test_make_defined_scape(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        dim = params.land.main.dim
        res = params.land.main.res
        ulc = params.land.main.ulc
        #TODO: The make defined scape has wrong number of parameters as called in the function
        #defined_scape = landscape.make_defined_scape(dim,)

    def test_make_land(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        dim = params.land.main.dim
        res = params.land.main.res
        ulc = params.land.main.ulc
        land = landscape._make_landscape(params, 2)
        self.assertIsInstance(land, landscape.Landscape)



if __name__ == '__main__':
    unittest.main()
